#include "stage1_range_resizer.h"
#include "stage1_config.h"

#include <fftw3.h>

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <limits>
#include <stdexcept>

namespace gmti {
namespace sim_stage1 {

RangeFftZeroPadResizer::RangeFftZeroPadResizer(int src_len, int dst_len)
    : src_len_(src_len),
      dst_len_(dst_len),
      in_(nullptr),
      spec_src_(nullptr),
      inv_in_(nullptr),
      time_full_(nullptr),
      forward_plan_(nullptr),
      inverse_plan_(nullptr),
      spec_shift_(static_cast<size_t>(src_len)),
      spec_pad_shift_(static_cast<size_t>(dst_len))
{
    if (src_len_ <= 0 || dst_len_ <= src_len_ || (src_len_ % 2) != 0 || (dst_len_ % 2) != 0) {
        throw std::invalid_argument("invalid range resize lengths");
    }
    in_ = fftwf_malloc(sizeof(fftwf_complex) * static_cast<size_t>(src_len_));
    spec_src_ = fftwf_malloc(sizeof(fftwf_complex) * static_cast<size_t>(src_len_));
    inv_in_ = fftwf_malloc(sizeof(fftwf_complex) * static_cast<size_t>(dst_len_));
    time_full_ = fftwf_malloc(sizeof(fftwf_complex) * static_cast<size_t>(dst_len_));
    if (!in_ || !spec_src_ || !inv_in_ || !time_full_) {
        cleanup();
        throw std::bad_alloc();
    }
    forward_plan_ = fftwf_plan_dft_1d(src_len_,
                                      static_cast<fftwf_complex *>(in_),
                                      static_cast<fftwf_complex *>(spec_src_),
                                      FFTW_FORWARD,
                                      FFTW_ESTIMATE);
    inverse_plan_ = fftwf_plan_dft_1d(dst_len_,
                                      static_cast<fftwf_complex *>(inv_in_),
                                      static_cast<fftwf_complex *>(time_full_),
                                      FFTW_BACKWARD,
                                      FFTW_ESTIMATE);
    if (!forward_plan_ || !inverse_plan_) {
        cleanup();
        throw std::runtime_error("failed to create FFTW plans");
    }
}

RangeFftZeroPadResizer::~RangeFftZeroPadResizer()
{
    cleanup();
}

void RangeFftZeroPadResizer::cleanup()
{
    if (forward_plan_) fftwf_destroy_plan(static_cast<fftwf_plan>(forward_plan_));
    if (inverse_plan_) fftwf_destroy_plan(static_cast<fftwf_plan>(inverse_plan_));
    fftwf_free(in_);
    fftwf_free(spec_src_);
    fftwf_free(inv_in_);
    fftwf_free(time_full_);
    forward_plan_ = nullptr;
    inverse_plan_ = nullptr;
    in_ = nullptr;
    spec_src_ = nullptr;
    inv_in_ = nullptr;
    time_full_ = nullptr;
}

bool RangeFftZeroPadResizer::resize(const std::complex<float> *src,
                                    std::vector<std::complex<float> > &dst)
{
    if (!src || !forward_plan_ || !inverse_plan_) return false;
    fftwf_complex *in = static_cast<fftwf_complex *>(in_);
    fftwf_complex *spec_src = static_cast<fftwf_complex *>(spec_src_);
    fftwf_complex *inv_in = static_cast<fftwf_complex *>(inv_in_);
    fftwf_complex *time_full = static_cast<fftwf_complex *>(time_full_);
    for (int n = 0; n < src_len_; ++n) {
        in[n][0] = src[n].real();
        in[n][1] = src[n].imag();
    }
    fftwf_execute(static_cast<fftwf_plan>(forward_plan_));
    for (int k = 0; k < src_len_; ++k) {
        const int src_index = (k + src_len_ / 2) % src_len_;
        spec_shift_[static_cast<size_t>(k)] =
            std::complex<float>(spec_src[src_index][0], spec_src[src_index][1]);
    }
    std::fill(spec_pad_shift_.begin(), spec_pad_shift_.end(),
              std::complex<float>(0.0f, 0.0f));
    const int insert = (dst_len_ - src_len_) / 2;
    std::copy(spec_shift_.begin(), spec_shift_.end(), spec_pad_shift_.begin() + insert);
    for (int k = 0; k < dst_len_; ++k) {
        const std::complex<float> v =
            spec_pad_shift_[static_cast<size_t>((k + dst_len_ / 2) % dst_len_)];
        inv_in[k][0] = v.real();
        inv_in[k][1] = v.imag();
    }
    fftwf_execute(static_cast<fftwf_plan>(inverse_plan_));
    dst.resize(static_cast<size_t>(dst_len_));
    const float scale = 1.0f / static_cast<float>(src_len_);
    for (int n = 0; n < dst_len_; ++n) {
        dst[static_cast<size_t>(n)] =
            std::complex<float>(time_full[n][0] * scale, time_full[n][1] * scale);
    }
    return true;
}

SignalStats computeStats(const std::vector<std::complex<float> > &x)
{
    SignalStats s;
    if (x.empty()) return s;
    double sum_abs = 0.0;
    double sum_abs2 = 0.0;
    for (size_t i = 0; i < x.size(); ++i) {
        const float re = x[i].real();
        const float im = x[i].imag();
        if (std::isnan(re) || std::isnan(im)) s.has_nan = true;
        if (std::isinf(re) || std::isinf(im)) s.has_inf = true;
        const double a = std::abs(x[i]);
        sum_abs += a;
        sum_abs2 += a * a;
        if (a > s.max_abs) s.max_abs = a;
    }
    s.mean_abs = sum_abs / static_cast<double>(x.size());
    s.power = sum_abs2;
    s.rms = std::sqrt(sum_abs2 / static_cast<double>(x.size()));
    double var = 0.0;
    for (size_t i = 0; i < x.size(); ++i) {
        const double d = std::abs(x[i]) - s.mean_abs;
        var += d * d;
    }
    s.std_abs = std::sqrt(var / static_cast<double>(x.size()));
    return s;
}

bool runRangeResizeSelfTests(const std::string &output_dir,
                             RangeFftZeroPadResizer &resizer,
                             std::string &err)
{
    const int n0 = 4096;
    std::ofstream md(pathJoin(pathJoin(output_dir, "reports"), "range_resize_validation.md").c_str());
    std::ofstream csv(pathJoin(pathJoin(output_dir, "debug"), "range_resize_stats.csv").c_str());
    if (!md || !csv) {
        err = "failed to open range resize validation outputs";
        return false;
    }
    csv << "case,input_mean_abs,input_rms,input_max_abs,output_mean_abs,output_rms,output_max_abs,has_nan,has_inf\n";
    md << "# range resize validation\n\n";
    md << "4096 -> 11820 uses FFT, fftshift, centered zero padding, ifftshift, and FFTW backward IFFT scaled by 1/src_len.\n\n";

    std::vector<std::complex<float> > in(static_cast<size_t>(n0));
    std::vector<std::complex<float> > out;
    const char *names[3] = {"constant", "single_tone", "random"};
    for (int t = 0; t < 3; ++t) {
        for (int i = 0; i < n0; ++i) {
            if (t == 0) in[static_cast<size_t>(i)] = std::complex<float>(1.0f, 0.0f);
            else if (t == 1) {
                const double ph = 2.0 * M_PI * 17.0 * i / static_cast<double>(n0);
                in[static_cast<size_t>(i)] =
                    std::complex<float>(static_cast<float>(std::cos(ph)), static_cast<float>(std::sin(ph)));
            } else {
                const float re = static_cast<float>((std::rand() % 2001 - 1000) / 1000.0);
                const float im = static_cast<float>((std::rand() % 2001 - 1000) / 1000.0);
                in[static_cast<size_t>(i)] = std::complex<float>(re, im);
            }
        }
        if (!resizer.resize(in.data(), out)) {
            err = "range resize selftest resize failed";
            return false;
        }
        const SignalStats a = computeStats(in);
        const SignalStats b = computeStats(out);
        csv << names[t] << "," << a.mean_abs << "," << a.rms << "," << a.max_abs << ","
            << b.mean_abs << "," << b.rms << "," << b.max_abs << ","
            << boolText(b.has_nan) << "," << boolText(b.has_inf) << "\n";
        md << "## " << names[t] << "\n\n"
           << "- input rms: " << a.rms << "\n"
           << "- output rms: " << b.rms << "\n"
           << "- output max_abs: " << b.max_abs << "\n"
           << "- nan/inf: " << boolText(b.has_nan) << "/" << boolText(b.has_inf) << "\n\n";
    }
    md << "This is an engineering interpolation for stage 1. It does not convert old 38 us LFM physics into strict new 130 us LFM echoes.\n";
    return true;
}

} // namespace sim_stage1
} // namespace gmti

