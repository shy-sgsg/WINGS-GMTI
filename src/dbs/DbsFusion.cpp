#include "dbs/DbsFusion.hpp"
#include "dbs/DbsStitcher.hpp"
#include "rotation_xy.hpp"
#include "unwrap_fd.hpp"
#include "trig_lut.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <iostream>
#include <sstream>
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
    angleDeg = gmti::trig_lut::asin(ratio) * 180.0 / M_PI;
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

static std::vector<int> parse_debug_beam_env()
{
    std::vector<int> out;
    const char *env = std::getenv("DBS_DEBUG_SINGLE_BEAMS");
    if (!env || !*env) {
        return out;
    }
    std::stringstream ss(env);
    std::string token;
    while (std::getline(ss, token, ',')) {
        try {
            const int v = std::stoi(token);
            if (v > 0) {
                out.push_back(v);
            }
        } catch (const std::exception &) {
        }
    }
    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end()), out.end());
    return out;
}

static bool debug_beam_matches(int requested,
                               size_t slot,
                               const FusionGroupContext &ctx)
{
    if (requested == static_cast<int>(slot) + 1) {
        return true;
    }
    if (slot < ctx.periodList.size() && requested == ctx.periodList[slot]) {
        return true;
    }
    if (slot < ctx.beam_meta.size() && requested == ctx.beam_meta[slot].beam_index) {
        return true;
    }
    return false;
}

static void amplitude_stats(const Image2D<float> &img,
                            size_t &valid,
                            double &mean,
                            double &p95,
                            double &p99,
                            double &maxv)
{
    valid = 0;
    mean = 0.0;
    p95 = 0.0;
    p99 = 0.0;
    maxv = 0.0;
    std::vector<float> vals;
    vals.reserve(img.buf.size());
    for (float v : img.buf) {
        const float a = std::fabs(v);
        if (!(a > 0.0f) || !std::isfinite(a)) {
            continue;
        }
        vals.push_back(a);
        mean += a;
        if (a > maxv) {
            maxv = a;
        }
    }
    valid = vals.size();
    if (valid == 0) {
        return;
    }
    mean /= static_cast<double>(valid);
    auto percentile = [&](double q) -> double {
        size_t idx = static_cast<size_t>(
            std::min<double>(valid - 1, std::floor(q * static_cast<double>(valid - 1))));
        std::nth_element(vals.begin(), vals.begin() + static_cast<std::ptrdiff_t>(idx),
                         vals.end());
        return vals[idx];
    };
    p95 = percentile(0.95);
    p99 = percentile(0.99);
}

static void log_dbs_diagnostics(const FusionGroupContext &ctx,
                                const Config &cfg,
                                const Grid *grid,
                                const Bounds *bounds)
{
    if (bounds && grid) {
        const size_t nx = grid->x.size();
        const size_t ny = grid->y.size();
        std::cout << "[fusion][dbs][diag][grid]"
                  << " boundsX=[" << bounds->minX << "," << bounds->maxX << "]"
                  << " boundsY=[" << bounds->minY << "," << bounds->maxY << "]"
                  << " nx=" << nx
                  << " ny=" << ny
                  << " dx=" << grid->dx
                  << " dy=" << grid->dy
                  << " pixels=" << (nx * ny)
                  << " max_pixels=" << cfg.dbs_max_mosaic_pixels
                  << std::endl;
    }

    const size_t n = ctx.meta.beams.size();
    std::cout << "[fusion][dbs][diag] beams=" << n
              << " loc_beam_gate_deg=" << location_beam_gate_deg(cfg)
              << " dbs_interp_mode=" << cfg.dbs_interp_mode
              << " dbs_beam_skip=" << cfg.dbs_beam_skip
              << " dbs_range_skip=" << cfg.dbs_range_skip
              << std::endl;
    for (size_t i = 0; i < n; ++i) {
        const MetaPerBeam &bm = ctx.meta.beams[i];
        const FusionBeamMeta *fm =
            (i < ctx.beam_meta.size()) ? &ctx.beam_meta[i] : nullptr;
        const std::vector<float> *fdv =
            (i < ctx.rd.fd_axis.size()) ? &ctx.rd.fd_axis[i] : nullptr;
        const float fd_first = (fdv && !fdv->empty()) ? fdv->front() : 0.0f;
        const float fd_mid = (fdv && !fdv->empty()) ? (*fdv)[fdv->size() / 2] : 0.0f;
        const float fd_last = (fdv && !fdv->empty()) ? fdv->back() : 0.0f;
        size_t valid = 0;
        double mean = 0.0, p95 = 0.0, p99 = 0.0, maxv = 0.0;
        if (i < ctx.rd.amp.size()) {
            amplitude_stats(ctx.rd.amp[i], valid, mean, p95, p99, maxv);
        }
        std::cout << "[fusion][dbs][diag][beam]"
                  << " slot=" << (i + 1)
                  << " period=" << (i < ctx.periodList.size() ? ctx.periodList[i] : -1)
                  << " beam_index=" << (fm ? fm->beam_index : -1)
                  << " angle_deg=" << bm.angle_deg
                  << " theta_sq=" << (fm ? fm->theta_sq : 0.0)
                  << " theta_true=" << (fm ? fm->theta_true : bm.angle_deg)
                  << " beam_pointing_bias=" << (fm ? fm->beam_pointing_bias_deg : 0.0)
                  << " fd_ctr_wrapped=" << (fm ? fm->fd_ctr_wrapped : 0.0)
                  << " fd_ctr_unwrapped=" << (fm ? fm->fd_ctr_unwrapped : bm.fd_ctr)
                  << " unwrap_k=" << (fm ? fm->prf_ambiguity_k : 0)
                  << " meta_fd_ctr=" << bm.fd_ctr
                  << " fd_axis[first,mid,last]=" << fd_first << ","
                  << fd_mid << "," << fd_last
                  << " amp_valid=" << valid
                  << " amp_mean=" << mean
                  << " amp_p95=" << p95
                  << " amp_p99=" << p99
                  << " amp_max=" << maxv
                  << std::endl;
    }
}

} // namespace

bool fusionSlotHasSignal(const FusionGroupContext &ctx, size_t slot)
{
    if (slot >= ctx.rd.amp.size() || ctx.rd.amp[slot].empty()) {
        return false;
    }
    const Image2D<float> &amp = ctx.rd.amp[slot];
    for (float v : amp.buf) {
        if (std::isfinite(v) && std::fabs(v) > 0.0f) {
            return true;
        }
    }
    return false;
}

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
        const double cosT = std::abs(gmti::trig_lut::cos(deg2rad(thetaRot)));
        const double sinT = std::abs(gmti::trig_lut::sin(deg2rad(thetaRot)));
        const double vE = m.plane.V * gmti::trig_lut::cos(deg2rad(m.plane.V_angle));
        const double vN = m.plane.V * gmti::trig_lut::sin(deg2rad(m.plane.V_angle));
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
            const double targetAzimuthDeg = normalize_azimuth_deg(gmti::trig_lut::atan2(dN, dE) * 180.0 / M_PI);
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

    FusionGroupContext activeCtx;
    for (size_t slot = 0; slot < ctx.meta.beams.size(); ++slot) {
        if (!fusionSlotHasSignal(ctx, slot)) {
            continue;
        }
        activeCtx.periodList.push_back(slot < ctx.periodList.size()
                                           ? ctx.periodList[slot]
                                           : static_cast<int>(slot + 1));
        activeCtx.rd.amp.push_back(ctx.rd.amp[slot]);
        activeCtx.rd.fd_axis.push_back(ctx.rd.fd_axis[slot]);
        activeCtx.rd.rg_axis.push_back(ctx.rd.rg_axis[slot]);
        activeCtx.meta.beams.push_back(ctx.meta.beams[slot]);
        if (slot < ctx.beam_meta.size()) {
            activeCtx.beam_meta.push_back(ctx.beam_meta[slot]);
        }
        if (slot < ctx.detections.size()) {
            activeCtx.detections.push_back(ctx.detections[slot]);
        }
        if (slot < ctx.done.size()) {
            activeCtx.done.push_back(ctx.done[slot]);
        }
    }
    activeCtx.rd.nEff = ctx.rd.nEff;
    if (activeCtx.meta.beams.empty()) {
        std::cerr << "[fusion][dbs] no active nonzero beam after zero-fill filtering" << std::endl;
        return false;
    }
    std::cout << "[fusion][dbs] active nonzero beams for mosaic:";
    for (int p : activeCtx.periodList) {
        std::cout << ' ' << p;
    }
    std::cout << " (" << activeCtx.meta.beams.size() << " / "
              << ctx.meta.beams.size() << ")" << std::endl;

    Bounds bounds;
    Grid grid;
    if (!stitcher.estimateMosaicExtent(cfg, activeCtx.rd, activeCtx.meta, bounds, grid)) {
        std::cerr << "[fusion][dbs] estimateMosaicExtent failed" << std::endl;
        return false;
    }
    log_dbs_diagnostics(ctx, cfg, nullptr, nullptr);
    log_dbs_diagnostics(activeCtx, cfg, &grid, &bounds);

    Mosaic mosaic;
    const bool built = useGpu
        ? stitcher.buildMosaicGPU(cfg, activeCtx.rd, activeCtx.meta, grid, mosaic, false)
        : stitcher.buildMosaic(cfg, activeCtx.rd, activeCtx.meta, grid, mosaic);
    if (!built) {
        std::cerr << "[fusion][dbs] buildMosaic failed" << std::endl;
        return false;
    }

    const std::vector<int> debugSingleBeams = parse_debug_beam_env();
    for (int requestedBeam : debugSingleBeams) {
        for (size_t slot = 0; slot < activeCtx.meta.beams.size(); ++slot) {
            if (!debug_beam_matches(requestedBeam, slot, activeCtx)) {
                continue;
            }
            RDData oneRd;
            oneRd.nEff = activeCtx.rd.nEff;
            oneRd.amp.push_back(activeCtx.rd.amp[slot]);
            oneRd.fd_axis.push_back(activeCtx.rd.fd_axis[slot]);
            oneRd.rg_axis.push_back(activeCtx.rd.rg_axis[slot]);
            MetaPack oneMeta;
            oneMeta.beams.push_back(activeCtx.meta.beams[slot]);
            Mosaic singleMosaic;
            if (!stitcher.buildMosaic(cfg, oneRd, oneMeta, grid, singleMosaic)) {
                std::cerr << "[fusion][dbs][debug] single-beam mosaic failed for requested beam "
                          << requestedBeam << " slot=" << (slot + 1) << std::endl;
                continue;
            }
            const int period = slot < activeCtx.periodList.size()
                                   ? activeCtx.periodList[slot]
                                   : requestedBeam;
            char name[128];
            std::snprintf(name, sizeof(name),
                          "DBS_debug_beam%d_slot%zu_period%d_on_current_grid.png",
                          requestedBeam, slot + 1, period);
            if (!stitcher.writeDebugMosaicImage(cfg, singleMosaic, name)) {
                std::cerr << "[fusion][dbs][debug] write single-beam debug image failed: "
                          << name << std::endl;
            }
        }
    }

    if (!stitcher.writeProducts(cfg, grid, mosaic, bounds, &activeCtx.meta)) {
        std::cerr << "[fusion][dbs] writeProducts failed" << std::endl;
        return false;
    }
    std::cout << "[fusion][dbs] writeProducts success:" << cfg.result_add << std::endl;

    return true;
}
