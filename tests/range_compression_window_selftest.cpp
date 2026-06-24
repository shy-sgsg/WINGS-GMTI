#include "config_structs.hpp"
#include "rangeCompress.hpp"

#include <cmath>
#include <complex>
#include <iostream>
#include <vector>

int main()
{
    Config cfg{};
    cfg.pulse_len = 11820;
    cfg.rg_len = 4096;
    cfg.range_fft_len = 12288;
    cfg.range_crop_start = 3864;
    cfg.range_compress_len = 4096;
    cfg.pulse_num = 1;
    cfg.process_pulse_num = 1;
    cfg.hasRefFunc = 0;
    cfg.fs = 225.0e6;
    cfg.Br = 50.0e6;
    cfg.Tr = 38.0e-6;

    std::vector<std::complex<double>> input(
        static_cast<size_t>(cfg.pulse_len),
        std::complex<double>(0.0, 0.0));
    input[0] = std::complex<double>(1.0, 0.0);

    std::vector<std::complex<double>> output;
    if (!rangeCompressFFT(cfg, input, output)) {
        std::cerr << "rangeCompressFFT returned false" << std::endl;
        return 1;
    }
    if (output.size() != static_cast<size_t>(cfg.range_compress_len)) {
        std::cerr << "output size mismatch: " << output.size() << std::endl;
        return 2;
    }
    for (size_t i = 0; i < output.size(); ++i) {
        if (!std::isfinite(output[i].real()) ||
            !std::isfinite(output[i].imag())) {
            std::cerr << "non-finite output at index " << i << std::endl;
            return 3;
        }
    }

    std::cout << "[SELFTEST] PASS: 11820 -> zero-pad 12288 -> crop [3864, 7960)"
              << " -> 4096 samples" << std::endl;
    return 0;
}
