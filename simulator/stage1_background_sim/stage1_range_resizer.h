#pragma once

#include <complex>
#include <string>
#include <vector>

namespace gmti {
namespace sim_stage1 {

struct SignalStats {
    double mean_abs = 0.0;
    double std_abs = 0.0;
    double rms = 0.0;
    double max_abs = 0.0;
    double power = 0.0;
    bool has_nan = false;
    bool has_inf = false;
};

class RangeFftZeroPadResizer {
public:
    RangeFftZeroPadResizer(int src_len, int dst_len);
    ~RangeFftZeroPadResizer();
    RangeFftZeroPadResizer(const RangeFftZeroPadResizer &) = delete;
    RangeFftZeroPadResizer &operator=(const RangeFftZeroPadResizer &) = delete;
    bool resize(const std::complex<float> *src, std::vector<std::complex<float> > &dst);

private:
    void cleanup();
    int src_len_;
    int dst_len_;
    void *in_;
    void *spec_src_;
    void *inv_in_;
    void *time_full_;
    void *forward_plan_;
    void *inverse_plan_;
    std::vector<std::complex<float> > spec_shift_;
    std::vector<std::complex<float> > spec_pad_shift_;
};

SignalStats computeStats(const std::vector<std::complex<float> > &x);
bool runRangeResizeSelfTests(const std::string &output_dir,
                             RangeFftZeroPadResizer &resizer,
                             std::string &err);

} // namespace sim_stage1
} // namespace gmti

