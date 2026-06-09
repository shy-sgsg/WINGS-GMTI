#include "dbs/DbsFusion.hpp"
#include "dbs/DbsStitcher.hpp"
#include "rotation_xy.hpp"
#include "unwrap_fd.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

namespace {

static bool estimate_angle_from_fd(const FusionBeamMeta &m, double &angleDeg)
{
    const double v = m.plane.V;
    const double lambda = m.lambda;
    if (!(v > 0.0) || !(lambda > 0.0)) {
        return false;
    }

    double ratio = -m.fd_ctr_wrapped * lambda / (2.0 * v);
    ratio = std::max(-1.0, std::min(1.0, ratio));
    angleDeg = std::asin(ratio) * 180.0 / M_PI;
    return std::isfinite(angleDeg);
}

static std::vector<size_t> center_indices(size_t n)
{
    std::vector<size_t> out;
    if (n == 0) {
        return out;
    }
    if (n == 1) {
        out.push_back(0);
        return out;
    }

    const size_t mid = n / 2;
    if ((n % 2) == 0) {
        out.push_back(mid - 1);
        out.push_back(mid);
    } else {
        out.push_back(mid);
        if (mid + 1 < n) {
            out.push_back(mid + 1);
        }
    }
    return out;
}

static inline int sgn(double x, double eps = 1e-6)
{
    return (x > eps) - (x < -eps);
}

static int flight_flag_by_sign_local(double vE, double vN, double eps = 1e-6)
{
    const int sE = sgn(vE, eps);
    const int sN = sgn(vN, eps);
    if (sE == 0 && sN == 0) {
        return 0;
    }
    if (sN < 0 && sE < 0) {
        return 1;
    }
    if (sN > 0 && sE > 0) {
        return 2;
    }
    if (sN > 0 && sE < 0) {
        return 3;
    }
    if (sN < 0 && sE > 0) {
        return 4;
    }
    if (sE == 0) {
        return (sN > 0) ? 3 : 4;
    }
    if (sN == 0) {
        return (sE > 0) ? 2 : 1;
    }
    return 0;
}

static inline double normalize_azimuth_deg(double angleDeg)
{
    angleDeg = std::fmod(angleDeg, 360.0);
    if (angleDeg < 0.0) {
        angleDeg += 360.0;
    }
    return angleDeg;
}

static inline double wrap180_deg(double angleDeg)
{
    angleDeg = std::fmod(angleDeg + 180.0, 360.0);
    if (angleDeg < 0.0) {
        angleDeg += 360.0;
    }
    return angleDeg - 180.0;
}

static inline double beam_center_relative_dir_deg(int squintSide, double thetaTrueDeg)
{
    const double sideDir = (squintSide == 1) ? -90.0 : 90.0;
    return wrap180_deg(sideDir - thetaTrueDeg);
}

static inline double location_beam_gate_deg(const Config &cfg)
{
    if (cfg.loc_beam_gate_deg > 0.0) {
        return cfg.loc_beam_gate_deg;
    }
    return std::max(0.0, cfg.beamwidth_deg * 0.5) + 0.5;
}

} // namespace

bool estimateBeamPointingBiasByCenterBeams(const std::vector<int> &periodList,
                                           const std::vector<FusionBeamMeta> &beamMeta,
                                           double &biasDeg)
{
    biasDeg = 0.0;
    if (periodList.empty() || beamMeta.size() != periodList.size()) {
        return false;
    }

    const std::vector<size_t> indices = center_indices(periodList.size());
    double sumBias = 0.0;
    int count = 0;

    for (size_t idx : indices) {
        const FusionBeamMeta &m = beamMeta[idx];
        double angleDeg = 0.0;
        if (!estimate_angle_from_fd(m, angleDeg)) {
            continue;
        }
        sumBias += (angleDeg - m.theta_sq);
        ++count;
    }

    if (count == 0) {
        return false;
    }

    biasDeg = sumBias / static_cast<double>(count);
    if (std::abs(biasDeg) > 5.0) {
        std::cerr << "[fusion][fd][warn] beam_pointing_bias exceeds 5 deg: "
                  << biasDeg << std::endl;
    }
    return std::isfinite(biasDeg);
}

bool applyBeamPointingBiasToFusionContext(double biasDeg,
                                          FusionGroupContext &ctx)
{
    const size_t n = ctx.periodList.size();
    if (ctx.beam_meta.size() != n || ctx.meta.beams.size() != n ||
        ctx.rd.fd_axis.size() != n) {
        return false;
    }

    for (size_t i = 0; i < n; ++i) {
        FusionBeamMeta &m = ctx.beam_meta[i];
        m.beam_pointing_bias_deg = biasDeg;
        m.theta_true = m.theta_sq + biasDeg;
        m.fd_ctr_unwrapped = unwrap_prf_to_model(m.fd_ctr_wrapped,
                                                 m.PRF,
                                                 m.theta_true,
                                                 m.plane.V,
                                                 m.fc_hz);
        if (m.PRF > 0.0) {
            m.prf_ambiguity_k = nearest_int((m.fd_ctr_unwrapped - m.fd_ctr_wrapped) / m.PRF);
        }

        MetaPerBeam &dbsMeta = ctx.meta.beams[i];
        dbsMeta.fd_ctr = static_cast<float>(m.fd_ctr_unwrapped);
        dbsMeta.angle_deg = static_cast<float>(m.theta_true);

        for (float &fd : ctx.rd.fd_axis[i]) {
            fd += static_cast<float>(m.fd_ctr_unwrapped);
        }

        // std::cout << "[fusion][fd] beam=" << m.beam_index
        //           << " theta_sq=" << m.theta_sq
        //           << " theta_true=" << m.theta_true
        //           << " fd_wrapped=" << m.fd_ctr_wrapped
        //           << " fd_unwrapped=" << m.fd_ctr_unwrapped
        //           << " fa_shift=" << (m.fd_ctr_unwrapped - m.fd_ctr_wrapped)
        //           << " k_prf=" << m.prf_ambiguity_k << std::endl;
    }

    return true;
}

bool relocateFusionDetections(const FusionGroupContext &ctx,
                              const Config &cfg,
                              std::vector<GMTIOutput> &results)
{
    const size_t n = ctx.periodList.size();
    if (ctx.beam_meta.size() != n || ctx.detections.size() != n) {
        return false;
    }

    results.clear();
    results.resize(n);

    auto deg2rad = [](double d) { return d * M_PI / 180.0; };

    size_t totalRaw = 0;
    size_t totalKept = 0;
    size_t totalDropSin = 0;
    size_t totalDropPx = 0;
    size_t totalDropGate = 0;

    for (size_t slot = 0; slot < n; ++slot) {
        const FusionBeamMeta &m = ctx.beam_meta[slot];
        const auto &rawList = ctx.detections[slot];
        GMTIOutput &out = results[slot];

        out.utcMid = m.utc_mid;
        out.plane = m.plane;
        out.detect.prow.reserve(rawList.size());
        out.detect.pcol.reserve(rawList.size());
        out.detect.row_af.reserve(rawList.size());
        out.MT.clear();
        out.MT.reserve(rawList.size() * 8);

        const double thetaRot = m.plane.V_angle;
        const double cosT = std::abs(std::cos(deg2rad(thetaRot)));
        const double sinT = std::abs(std::sin(deg2rad(thetaRot)));
        const double vE = m.plane.V * std::cos(deg2rad(m.plane.V_angle));
        const double vN = m.plane.V * std::sin(deg2rad(m.plane.V_angle));
        int flag = flight_flag_by_sign_local(vE, vN);
        flag += cfg.squint_side * 4;

        const double lambda = (m.lambda > 0.0) ? m.lambda : cfg.lambda;
        if (!(lambda > 0.0) || !(m.plane.V > 0.0)) {
            continue;
        }

        const double faShift = m.fd_ctr_unwrapped - m.fd_ctr_wrapped;
        size_t kept = 0;
        size_t dropSin = 0;
        size_t dropPx = 0;
        size_t dropGate = 0;
        for (const DetectionRaw &d : rawList) {
            double afRansac = d.af_wrapped;
            if (std::abs(m.phase_slope) >= 1e-12) {
                afRansac = (d.phase - m.phase_intercept) / m.phase_slope;
            }
            const double af = afRansac + faShift;
            const double sinA = af * lambda / (2.0 * m.plane.V);
            if (std::abs(sinA) > 1.0) {
                ++dropSin;
                continue;
            }

            const double Rg = d.range_m;
            const double py = Rg * sinA;
            const double dz = m.plane.H - cfg.MT_nowz;
            const double px2 = Rg * Rg - py * py - dz * dz;
            if (px2 < 0.0) {
                ++dropPx;
                continue;
            }
            const double px = (px2 > 0.0) ? std::sqrt(px2) : 0.0;

            double xP = 0.0;
            double yP = 0.0;
            rotation_xy(py, px, flag, cosT, sinT, xP, yP);
            xP += m.plane.E;
            yP += m.plane.N;

            double lat = 0.0;
            double lng = 0.0;
            (void)Gaussp3RV(xP, yP, cfg.L0, lat, lng);

            const double dE = xP - m.plane.E;
            const double dN = yP - m.plane.N;
            const double targetAzimuthDeg = normalize_azimuth_deg(std::atan2(dN, dE) * 180.0 / M_PI);
            const double direction = wrap180_deg(m.plane.V_angle - targetAzimuthDeg);
            const double beamCenterDir = beam_center_relative_dir_deg(cfg.squint_side, m.theta_true);
            const double beamHalfWidth = location_beam_gate_deg(cfg);
            const double beamDirErr = wrap180_deg(direction - beamCenterDir);
            if (std::abs(beamDirErr) > beamHalfWidth) {
                // std::cout << "[fusion][loc][gate] drop target beam=" << m.beam_index
                //           << " dir=" << direction
                //           << " beam_center_dir=" << beamCenterDir
                //           << " theta_true=" << m.theta_true
                //           << " err=" << beamDirErr
                //           << " gate=" << beamHalfWidth << std::endl;
                ++dropGate;
                continue;
            }
            const double range = std::sqrt(dE * dE + dN * dN);

            out.detect.prow.push_back(d.prow);
            out.detect.pcol.push_back(d.pcol);
            out.detect.row_af.push_back(af);

            out.MT.push_back(lat);
            out.MT.push_back(lng);
            out.MT.push_back(cfg.MT_nowz);
            out.MT.push_back(xP);
            out.MT.push_back(yP);
            out.MT.push_back(m.utc_mid);
            out.MT.push_back(direction);
            out.MT.push_back(range);
            ++kept;
        }

        totalRaw += rawList.size();
        totalKept += kept;
        totalDropSin += dropSin;
        totalDropPx += dropPx;
        totalDropGate += dropGate;
        // std::cout << "[fusion][loc][summary] beam=" << m.beam_index
        //           << " raw=" << rawList.size()
        //           << " kept=" << kept
        //           << " drop_sin=" << dropSin
        //           << " drop_px=" << dropPx
        //           << " drop_gate=" << dropGate
        //           << " theta_true=" << m.theta_true
        //           << " fd_shift=" << faShift << std::endl;
    }

    std::cout << "[fusion][loc][total] raw=" << totalRaw
              << " kept=" << totalKept
              << " drop_sin=" << totalDropSin
              << " drop_px=" << totalDropPx
              << " drop_gate=" << totalDropGate << std::endl;

    return true;
}

bool runDbsFusionImaging(const FusionGroupContext &ctx,
                         const Config &cfg,
                         bool useGpu)
{
    if (ctx.rd.amp.empty() || ctx.meta.beams.empty()) {
        std::cerr << "[fusion][dbs] empty RD/meta cache" << std::endl;
        return false;
    }

    const size_t beamCount = ctx.rd.amp.size();
    const size_t rdRows = ctx.rd.amp.empty() ? 0U : static_cast<size_t>(ctx.rd.amp.front().rows);
    const size_t rdCols = ctx.rd.amp.empty() ? 0U : static_cast<size_t>(ctx.rd.amp.front().cols);
    const double rdCacheGiB = static_cast<double>(beamCount) * static_cast<double>(rdRows) *
                              static_cast<double>(rdCols) * static_cast<double>(sizeof(float)) /
                              (1024.0 * 1024.0 * 1024.0);
    std::cout << "[fusion][dbs] imaging start: beams=" << beamCount
              << " rd=" << rdRows << "x" << rdCols
              << " rd_amp_cache~" << rdCacheGiB << " GiB"
              << " PRF=" << (ctx.beam_meta.empty() ? cfg.PRF : ctx.beam_meta.front().PRF)
              << " Rmin=" << cfg.R_min
              << " fs=" << ((cfg.R_bin > 0.0) ? (C / (2.0 * cfg.R_bin)) : cfg.fs)
              << " out_res=" << cfg.dbs_out_res_m
              << " range_skip=" << cfg.dbs_range_skip
              << " beam_skip=" << cfg.dbs_beam_skip
              << " useGpu=" << (useGpu ? 1 : 0) << std::endl;

    DbsStitcher stitcher;
    stitcher.setLonRef(cfg.L0);

    Bounds bounds;
    Grid grid;
    if (!stitcher.estimateMosaicExtent(cfg, ctx.rd, ctx.meta, bounds, grid)) {
        std::cerr << "[fusion][dbs] estimateMosaicExtent failed" << std::endl;
        return false;
    }

    Mosaic mosaic;
    const bool built = useGpu
        ? stitcher.buildMosaicGPU(cfg, ctx.rd, ctx.meta, grid, mosaic, false)
        : stitcher.buildMosaic(cfg, ctx.rd, ctx.meta, grid, mosaic);
    if (!built) {
        std::cerr << "[fusion][dbs] buildMosaic failed" << std::endl;
        return false;
    }

    if (!stitcher.writeProducts(cfg, grid, mosaic, bounds)) {
        std::cerr << "[fusion][dbs] writeProducts failed" << std::endl;
        return false;
    }
    std::cout << "[fusion][dbs] writeProducts success:" << cfg.result_add << std::endl;

    return true;
}
