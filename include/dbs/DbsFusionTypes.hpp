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
    double p38_pre_k = std::numeric_limits<double>::quiet_NaN();
    double p38_pre_b = std::numeric_limits<double>::quiet_NaN();
    double p38_pre_rmse = std::numeric_limits<double>::quiet_NaN();
    double p38_refit_k = std::numeric_limits<double>::quiet_NaN();
    double p38_refit_b = std::numeric_limits<double>::quiet_NaN();
    double p38_refit_rmse = std::numeric_limits<double>::quiet_NaN();
    int p38_refit_sample_count = 0;
    double p38_refit_inlier_ratio = std::numeric_limits<double>::quiet_NaN();
    int p38_refit_valid = 0;
    double p38_used_k = std::numeric_limits<double>::quiet_NaN();
    double p38_used_b = std::numeric_limits<double>::quiet_NaN();
    std::string p38_used_source;

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
    double af_row = 0.0;
    double af_phase = 0.0;
    double af_total = 0.0;
    double af_motion = 0.0;
    double af_geometry = 0.0;
    double phase = 0.0;
    double phi_motion = 0.0;
    double v_radial = 0.0;
    int motion_comp_valid = 0;
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
