#ifndef DBS_FUSION_TYPES_HPP
#define DBS_FUSION_TYPES_HPP

#include "config_structs.hpp"
#include "DbsTypes.hpp"

#include <vector>

struct FusionBeamMeta {
    int beam_index = 0;
    int slot = -1;

    double theta_sq = 0.0;
    double theta_true = 0.0;
    double beam_pointing_bias_deg = 0.0;

    double fd_ctr_wrapped = 0.0;
    double fd_ctr_unwrapped = 0.0;
    int prf_ambiguity_k = 0;

    double phase_slope = 0.0;
    double phase_intercept = 0.0;
    int az_center = 0;

    double utc_mid = 0.0;
    GMTIOutput::Plane plane;

    double PRF = 0.0;
    double fc_hz = 0.0;
    double lambda = 0.0;
};

struct DetectionRaw {
    int beam_index = 0;
    int slot = -1;
    int prow = 0;
    int pcol = 0;

    double range_m = 0.0;
    double af_wrapped = 0.0;
    double phase = 0.0;
    double amplitude = 0.0;
    double utc_mid = 0.0;
};

struct FusionGroupContext {
    std::vector<int> periodList;
    RDData rd;
    MetaPack meta;
    std::vector<FusionBeamMeta> beam_meta;
    std::vector<std::vector<DetectionRaw>> detections;
    std::vector<char> done;

    void reset(const std::vector<int> &periods) {
        periodList = periods;
        const size_t n = periods.size();
        rd = RDData();
        rd.amp.resize(n);
        rd.fd_axis.resize(n);
        rd.rg_axis.resize(n);
        meta = MetaPack();
        meta.beams.resize(n);
        beam_meta.assign(n, FusionBeamMeta());
        detections.assign(n, std::vector<DetectionRaw>());
        done.assign(n, 0);
        for (size_t i = 0; i < n; ++i) {
            beam_meta[i].beam_index = periods[i];
            beam_meta[i].slot = static_cast<int>(i);
        }
    }
};

#endif // DBS_FUSION_TYPES_HPP
