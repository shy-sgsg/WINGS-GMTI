#include "GMTIProcessor.hpp"
#include "geo/geoProj.hpp"
#include <vector>
#include <array>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>

static void print_xy_preview(
    const std::vector<std::array<double,2>>& dets,
    int maxItems)
{
    const int n = std::min<int>(static_cast<int>(dets.size()), std::max(1, maxItems));
    for (int i = 0; i < n; ++i) {
        std::cout << " (" << dets[i][0] << ", " << dets[i][1] << ")";
    }
    if (static_cast<int>(dets.size()) > n) std::cout << " ...";
}

// 中位数（忽略 NaN），仅用于小型数组
static double median_of(std::vector<double> v) {
    v.erase(std::remove_if(v.begin(), v.end(),
                           [](double x){ return !(x==x); }), v.end()); // 去 NaN
    if (v.empty()) return std::nan(""); 
    std::sort(v.begin(), v.end());
    const size_t n = v.size();
    return (n % 2) ? v[n/2] : 0.5*(v[n/2 - 1] + v[n/2]);
}

// 解析单帧 res.MT（stride=6/8）→ 本帧的 {x,y} 集合，并取出该帧 utc
static bool parse_one_frame_MT(
    const std::vector<double>& MT,               // 扁平：6*N
    std::vector<std::array<double,2>>& xy_frame, // 输出：{x,y}
    std::vector<double>& dir_frame,
    std::vector<double>& range_frame,
    double& utc                                   // 输出：该帧 utcMid
){
    xy_frame.clear();
    dir_frame.clear();
    range_frame.clear();
    utc = std::numeric_limits<double>::quiet_NaN();

    if (MT.empty()) {
        // 允许空帧：无点 && utc=NaN（后续估计 T 时会自动忽略）
        return true;
    }
    const size_t stride = (MT.size() % 8 == 0) ? 8 : 6;
    if (MT.size() % stride != 0) {
        std::cerr << "parse_one_frame_MT: MT.size() 非 6/8 的倍数，size=" << MT.size() << "\n";
        return false;
    }
    const size_t N = MT.size() / stride;
    xy_frame.reserve(N);
    dir_frame.reserve(N);
    range_frame.reserve(N);

    for (size_t i=0; i<N; ++i) {
        const size_t base = i * stride;
        const double utc_i = MT[base + 5];
        const double xP = MT[base + 3];
        const double yP = MT[base + 4];

        // 同一周期内 utc 相同，取第一个有效值作为该帧时间
        if (!(utc == utc) && std::isfinite(utc_i)) {
            utc = utc_i;
        }

        // 允许用 NaN 哨兵保存“仅时间无检测”的帧，避免生成伪检测点
        if (!std::isfinite(xP) || !std::isfinite(yP)) {
            continue;
        }

        xy_frame.push_back(std::array<double,2>{xP, yP});
        if (stride == 8) {
            dir_frame.push_back(MT[base + 6]);
            range_frame.push_back(MT[base + 7]);
        } else {
            dir_frame.push_back(0.0);
            range_frame.push_back(0.0);
        }
        // 其他维度 lat/lng/z 暂不用于跟踪
    }
    return true;
}

bool GMTIProcessor::buildMTposAndT(
    const std::vector<std::vector<double>>& MT_frames,
    const Config& cfg,
    std::vector<std::vector<std::array<double,2>>>& MT_pos,
    std::vector<std::vector<double>>& MT_dir,
    std::vector<std::vector<double>>& MT_range,
    double& T_out)
{
    MT_pos.clear();
    MT_dir.clear();
    MT_range.clear();
    T_out = 0.0;

    const size_t K = MT_frames.size();
    if (K == 0) {
        // 无帧：给个兜底 T
        T_out = (cfg.PRF > 0 && cfg.pulse_num > 0)
              ? (double(cfg.pulse_num) / cfg.PRF)
              : 1.0;
        return true;
    }

    // 先解析每帧，并收集 utc
    std::vector<double> utcs; utcs.reserve(K);
    MT_pos.resize(K);
    MT_dir.resize(K);
    MT_range.resize(K);
    for (size_t k=0; k<K; ++k) {
        double utc_k = 0.0;
        if (!parse_one_frame_MT(MT_frames[k], MT_pos[k], MT_dir[k], MT_range[k], utc_k)) {
            return false;
        }
        utcs.push_back(utc_k);
    }

    const bool dbg = (cfg.track_debug_level > 0);
    const int dbgMaxFrames = std::max(1, cfg.track_debug_frames);
    const int dbgMaxItems = std::max(1, cfg.track_debug_points);

    // 若 utc 非单调，按 utc 排序（保持帧同步）
    std::vector<size_t> idx(K);
    for (size_t i=0;i<K;++i) idx[i] = i;
    if (cfg.mt_sort_by_utc) {
        std::stable_sort(idx.begin(), idx.end(),
            [&](size_t a, size_t b){
                const bool va = std::isfinite(utcs[a]);
                const bool vb = std::isfinite(utcs[b]);
                if (va != vb) return va; // 有效 utc 排在前面
                if (!va) return a < b;   // 都无效时保持稳定
                return utcs[a] < utcs[b];
            });
    }

    bool sorted = true;
    for (size_t i=1;i<K;++i) if (idx[i] < idx[i-1]) { sorted=false; break; }
    if (cfg.mt_sort_by_utc && !sorted) {
        // 生成排序后的副本
        std::vector<std::vector<std::array<double,2>>> MT_pos_sorted; MT_pos_sorted.reserve(K);
        std::vector<std::vector<double>> MT_dir_sorted; MT_dir_sorted.reserve(K);
        std::vector<std::vector<double>> MT_range_sorted; MT_range_sorted.reserve(K);
        std::vector<double> utcs_sorted; utcs_sorted.reserve(K);
        for (size_t i=0;i<K;++i){
            MT_pos_sorted.push_back(MT_pos[idx[i]]);
            MT_dir_sorted.push_back(MT_dir[idx[i]]);
            MT_range_sorted.push_back(MT_range[idx[i]]);
            utcs_sorted.push_back(utcs[idx[i]]);
        }
        MT_pos.swap(MT_pos_sorted);
        MT_dir.swap(MT_dir_sorted);
        MT_range.swap(MT_range_sorted);
        utcs.swap(utcs_sorted);
    }

    if (dbg) {
        std::cout << "[TRACK-INPUT] frames=" << K
                  << " mt_sort_by_utc=" << (cfg.mt_sort_by_utc ? "true" : "false")
                  << " reordered=" << ((cfg.mt_sort_by_utc && !sorted) ? "true" : "false")
                  << std::endl;

        const size_t nShow = std::min<size_t>(K, static_cast<size_t>(dbgMaxFrames));
        for (size_t i = 0; i < nShow; ++i) {
            const size_t src = (cfg.mt_sort_by_utc ? idx[i] : i);
            std::cout << "[TRACK-INPUT] procFrame=" << i
                      << " <- srcFrame=" << src
                      << " utc=" << utcs[i]
                      << " detN=" << MT_pos[i].size()
                      << " firstXY:";
            print_xy_preview(MT_pos[i], dbgMaxItems);
            std::cout << std::endl;
        }
    }

    // 用 utc 的相邻差的**中位数**估计 T（忽略 <=0 的差值）
    std::vector<double> diffs; diffs.reserve(K>1?K-1:0);
    for (size_t i=1;i<K;++i) {
        if (!std::isfinite(utcs[i]) || !std::isfinite(utcs[i-1])) {
            continue;
        }
        const double d = utcs[i] - utcs[i-1];
        if (d > 0) diffs.push_back(d);
    }
    double T_est = median_of(diffs);

    // 兜底：若估不出来，用 cfg 的脉冲参数推算
    if (!(T_est == T_est) || T_est <= 0.0) {
        if (cfg.PRF > 0 && cfg.pulse_num > 0) {
            T_est = double(cfg.pulse_num) / cfg.PRF;
        } else {
            T_est = 1.0; // 最后兜底
        }
    }

    T_out = T_est;
    if (dbg) {
        std::cout << "[TRACK-INPUT] T_est=" << T_out
                  << " (positive-diff-count=" << diffs.size() << ")"
                  << std::endl;
    }
    return true;
}

bool GMTIProcessor::trackFromMT(
    const std::vector<std::vector<double>>& MT_frames,
    const Config& cfg,
    double v_max, double sigma_thresh, int min_len, int max_gap,
    std::vector<Track>& tracks_out, bool bypassTracking)
{
    std::vector<std::vector<std::array<double,2>>> MT_pos;
    std::vector<std::vector<double>> MT_dir;
    std::vector<std::vector<double>> MT_range;
    double T = 0.0;
    if (!buildMTposAndT(MT_frames, cfg, MT_pos, MT_dir, MT_range, T)) return false;
    // T = 6.65625;
    return kalmanTrack(MT_pos, MT_dir, MT_range, T, v_max, sigma_thresh, min_len, max_gap, cfg, tracks_out, bypassTracking);
}

// 2x2 逆
static inline bool inv2(const double S[4], double S_inv[4]) {
    const double det = S[0]*S[3] - S[1]*S[2];
    if (std::fabs(det) < 1e-12) return false;
    const double id = 1.0/det;
    S_inv[0] =  S[3]*id;
    S_inv[1] = -S[1]*id;
    S_inv[2] = -S[2]*id;
    S_inv[3] =  S[0]*id;
    return true;
}

bool GMTIProcessor::kalmanTrack(
    const std::vector<std::vector<std::array<double,2>>>& MT_pos,
    const std::vector<std::vector<double>>& MT_dir,
    const std::vector<std::vector<double>>& MT_range,
    double T, double v_max, double sigma_thresh, int min_len, int max_gap,
    const Config& cfg, std::vector<Track>& tracks_out, bool bypassTracking)
{
    tracks_out.clear();
    const int K = static_cast<int>(MT_pos.size());
    if (K <= 0) return true; // 无数据也算成功

    const bool dbg = (cfg.track_debug_level > 0);
    const int dbgMaxFrames = std::max(1, cfg.track_debug_frames);
    const int dbgMaxItems = std::max(1, cfg.track_debug_points);

    // ==== 在初始化 F, H, Q 之前插入 ====
    if (bypassTracking) {
        std::cout << "Bypass Kalman Tracking: Generating tracks directly from detections without filtering.\n";
        tracks_out.clear();
        int globalId = 1;
        for (size_t k = 0; k < MT_pos.size(); ++k) {
            for (const auto& det : MT_pos[k]) {
                Track tr;
                tr.id = globalId++; 
                tr.pos.push_back(det);        // 存入原始坐标 [x, y]
                tr.kf.push_back(det);         // 此时滤波值等于原始值
                tr.time.push_back(k * T);     // 时间戳
                tr.x = {{ det[0], det[1], 0.0, 0.0 }}; // 速度置 0
                if (k < MT_dir.size() && tr.pos.size() - 1 < MT_dir[k].size()) {
                    tr.direction = MT_dir[k][tr.pos.size() - 1];
                }
                if (k < MT_range.size() && tr.pos.size() - 1 < MT_range[k].size()) {
                    tr.range = MT_range[k][tr.pos.size() - 1];
                }
                
                // 为了兼容后续 write 函数，可能需要初始化一些基础成员
                tr.missed = 0;
                tr.last = static_cast<int>(k + 1);
                
                tracks_out.push_back(tr);
            }
        }
        return true; // 直接返回，跳过后面的卡尔曼循环和筛选
    }

    if (dbg) {
        std::cout << "[TRACK-KF] start K=" << K
                  << " T=" << T
                  << " v_max=" << v_max
                  << " sigma_thresh=" << sigma_thresh
                  << " min_len=" << min_len
                  << " max_gap=" << max_gap
                  << std::endl;
    }

    // F, H, Q, R
    std::array<double,16> F; // [ [1 0 T 0], [0 1 0 T], [0 0 1 0], [0 0 0 1] ]
    for (int i=0;i<16;++i) F[i]=0.0;
    F[0]=1; F[5]=1; F[10]=1; F[15]=1;
    F[2]=T; F[7]=T;

    const double q = 0.1;
    std::array<double,16> Q;
    for (int i=0;i<16;++i) Q[i]=0.0;
    Q[0]=Q[5]=Q[10]=Q[15]=q;

    const double r = 10.0;
    double R[4] = {r,0,0,r};

    // ==== 初始化：第1帧全部作为新轨迹 ====
    const std::vector<std::array<double,2>>& det0 = MT_pos[0];
    tracks_out.reserve(det0.size()+8);
    for (size_t i=0;i<det0.size();++i){
        Track tr;
        tr.x = {{ det0[i][0], det0[i][1], 0.0, 0.0 }};
        // P = 10*I
        for (int k=0;k<16;++k) tr.P[k]=0.0;
        tr.P[0]=tr.P[5]=tr.P[10]=tr.P[15]=10.0;

        tr.pos.push_back(det0[i]);
        tr.kf.push_back(det0[i]);
        tr.x_state.push_back(tr.x);
        tr.time.push_back(0.0);
        tr.last = 1;
        tr.missed = 0;
        tr.id = 0;
        tracks_out.push_back(tr);
    }

    if (dbg) {
        std::cout << "[TRACK-KF] init frame=0 detN=" << det0.size() << " firstXY:";
        print_xy_preview(det0, dbgMaxItems);
        std::cout << std::endl;
    }

    // ==== 帧间循环 ====
    for (int k=2; k<=K; ++k) {
        const std::vector<std::array<double,2>>& det = MT_pos[k-1];
        const int M = static_cast<int>(det.size());
        std::vector<char> assigned(M, 0);
        int frameSkipMaxGap = 0;
        int frameNoDet = 0;
        int frameSingular = 0;
        int frameAccept = 0;
        int frameRejectPos = 0;
        int frameRejectVel = 0;
        int frameRejectBoth = 0;

        if (dbg && k <= dbgMaxFrames) {
            std::cout << "[TRACK-KF] frame=" << (k-1)
                      << " detN=" << M
                      << " activeTracks=" << tracks_out.size()
                      << " firstXY:";
            print_xy_preview(det, dbgMaxItems);
            std::cout << std::endl;
        }

        // ========== 1) 预测步 ==========
        // 对所有活跃轨迹执行Kalman预测（无条件）
        // 关键：即使前一帧miss，该轨迹在本帧仍需要预测
        // - 状态被预测：x_pred = F * x_prev
        // - 协方差膨胀：P_pred = F * P_prev * F^T + Q
        // 这样连续miss时，P矩阵会持续膨胀，直到轨迹重新关联或被删除（max_gap超限）
        // Matlab对齐确认：与MATLAB原版逻辑一致
        for (size_t i=0;i<tracks_out.size();++i){
            // x <- F x
            std::array<double,4> xp;
            mat4_mul_vec4(F, tracks_out[i].x, xp);
            tracks_out[i].x = xp;

            // P <- F P F^T + Q
            std::array<double,16> FP, FPFt;
            mat4_mul(F, tracks_out[i].P, FP);
            // 乘以 F^T
            // F^T 仅 6 个非零： (0,0)=1 (1,1)=1 (2,0)=T (2,2)=1 (3,1)=T (3,3)=1
            for (int r=0;r<4;++r){
                // c=0
                FPFt[r*4+0] = FP[r*4+0]*1.0 + FP[r*4+2]*T;
                // c=1
                FPFt[r*4+1] = FP[r*4+1]*1.0 + FP[r*4+3]*T;
                // c=2
                FPFt[r*4+2] = FP[r*4+2]*1.0;
                // c=3
                FPFt[r*4+3] = FP[r*4+3]*1.0;
            }
            tracks_out[i].P = FPFt;
            mat4_add_inplace(tracks_out[i].P, Q);
        }

        // ========== 2) 数据关联与更新 ==========
        // 对活跃轨迹执行最近邻关联和Kalman更新（有条件）
        // 关键：只有关联成功才更新状态；失败时x和P保持预测值/膨胀值
        for (size_t i=0;i<tracks_out.size();++i){
            if (tracks_out[i].missed > max_gap) {
                frameSkipMaxGap += 1;
                if (dbg && k <= dbgMaxFrames && static_cast<int>(i) < dbgMaxItems) {
                    std::cout << "[TRACK-KF] skip frame=" << (k-1)
                              << " trIdx=" << i
                              << " reason=missed_exceed_max_gap"
                              << " missed=" << tracks_out[i].missed
                              << " max_gap=" << max_gap
                              << std::endl;
                }
                continue;
            }

            if (M <= 0) {
                tracks_out[i].missed += 1;
                tracks_out[i].last = k;
                frameNoDet += 1;
                if (dbg && k <= dbgMaxFrames && static_cast<int>(i) < dbgMaxItems) {
                    std::cout << "[TRACK-KF] miss frame=" << (k-1)
                              << " trIdx=" << i
                              << " reason=no_detection_in_frame"
                              << std::endl;
                }
                continue;
            }

            // 观测预测 z_pred = H x = [x y]
            double z_pred[2] = { tracks_out[i].x[0], tracks_out[i].x[1] };

            // MATLAB对齐：从所有检测中找最近点（不剔除 assigned）
            double best_d = std::numeric_limits<double>::infinity();
            int best_j = -1;
            for (int j=0;j<M;++j){
                const double dx = det[j][0] - z_pred[0];
                const double dy = det[j][1] - z_pred[1];
                const double d = std::sqrt(dx*dx + dy*dy);
                if (d < best_d) { best_d = d; best_j = j; }
            }

            if (dbg && k <= dbgMaxFrames && static_cast<int>(i) < dbgMaxItems) {
                std::cout << "[TRACK-KF] cand frame=" << (k-1)
                          << " trIdx=" << i
                          << " missed=" << tracks_out[i].missed
                          << " pred=(" << z_pred[0] << "," << z_pred[1] << ")"
                          << " unassignedDet=" << (M - std::count(assigned.begin(), assigned.end(), 1))
                          << " list:";
                std::vector<std::pair<double,int>> cand;
                cand.reserve(M);
                for (int j=0; j<M; ++j) {
                    const double dx = det[j][0] - z_pred[0];
                    const double dy = det[j][1] - z_pred[1];
                    const double d = std::sqrt(dx*dx + dy*dy);
                    cand.push_back(std::make_pair(d, j));
                }
                std::sort(cand.begin(), cand.end(),
                          [](const std::pair<double,int>& a, const std::pair<double,int>& b){
                              return a.first < b.first;
                          });
                const int shown = std::min<int>(dbgMaxItems, static_cast<int>(cand.size()));
                for (int c=0; c<shown; ++c) {
                    const int j = cand[c].second;
                    std::cout << " [j=" << j << ",d=" << cand[c].first
                              << (assigned[j] ? "*" : "") << "]";
                }
                if (best_j < 0) std::cout << " none";
                std::cout << std::endl;
            }

            if (best_j >= 0) {
                // 计算隐含速度（新检测相对轨迹最后点）
                const double impl_dx = det[best_j][0] - tracks_out[i].pos.back()[0];
                const double impl_dy = det[best_j][1] - tracks_out[i].pos.back()[1];
                const double implied_v = std::sqrt(impl_dx*impl_dx + impl_dy*impl_dy) / T;

                if (dbg && k <= dbgMaxFrames && static_cast<int>(i) < dbgMaxItems) {
                    std::cout << "[TRACK-KF] assoc frame=" << (k-1)
                              << " trIdx=" << i
                              << " pred=(" << z_pred[0] << "," << z_pred[1] << ")"
                              << " bestDet=" << best_j
                              << " bestD=" << best_d
                              << " impliedV=" << implied_v
                              << std::endl;
                }
                
                // 双判据：位置波门 AND 速度不超标
                if (best_d < sigma_thresh && implied_v < v_max) {
                    // 量测
                    double z[2] = { det[best_j][0], det[best_j][1] };

                    // S = H P H^T + R = P(0:1,0:1) + R
                    const std::array<double,16>& P = tracks_out[i].P;
                    double S[4] = {
                        P[0] + R[0], P[1] + R[1],
                        P[4] + R[2], P[5] + R[3]
                    };
                    double S_inv[4];
                    if (!inv2(S, S_inv)) {
                        // 奇异则跳过更新，当作 missed
                        tracks_out[i].missed += 1;
                        frameSingular += 1;
                        if (dbg && k <= dbgMaxFrames && static_cast<int>(i) < dbgMaxItems) {
                            std::cout << "[TRACK-KF] singular frame=" << (k-1)
                                      << " trIdx=" << i
                                      << " S=[" << S[0] << "," << S[1] << ";"
                                      << S[2] << "," << S[3] << "]"
                                      << std::endl;
                        }
                        continue;
                    }

                    // K = P H^T S^{-1}; H^T 取 P 的前两列 => PHt = [P(:,0) P(:,1)] (4x2)
                    double PHt[8] = {
                        P[0], P[1],
                        P[4], P[5],
                        P[8], P[9],
                        P[12],P[13]
                    };
                    // K = PHt * S_inv (4x2)
                    double K_gain[8];
                    for (int r=0;r<4;++r){
                        K_gain[r*2+0] = PHt[r*2+0]*S_inv[0] + PHt[r*2+1]*S_inv[2];
                        K_gain[r*2+1] = PHt[r*2+0]*S_inv[1] + PHt[r*2+1]*S_inv[3];
                    }

                    // 更新 x：x += K*(z - Hx)
                    double y0 = z[0] - tracks_out[i].x[0];
                    double y1 = z[1] - tracks_out[i].x[1];
                    std::array<double,4> dx;
                    dx[0] = K_gain[0]*y0 + K_gain[1]*y1;
                    dx[1] = K_gain[2]*y0 + K_gain[3]*y1;
                    dx[2] = K_gain[4]*y0 + K_gain[5]*y1;
                    dx[3] = K_gain[6]*y0 + K_gain[7]*y1;
                    vec4_add_inplace(tracks_out[i].x, dx);

                    // 更新 P：P = (I - K H) P  = P - K*(H P)
                    // 先 HP (2x4)：取 P 的前两行
                    double HP[8] = {
                        P[0], P[1], P[2], P[3],
                        P[4], P[5], P[6], P[7]
                    };
                    // KHP (4x4) = K(4x2) * HP(2x4)
                    std::array<double,16> KHP;
                    for (int r=0;r<4;++r){
                        for (int c=0;c<4;++c){
                            KHP[r*4+c] = K_gain[r*2+0]*HP[0*4+c] + K_gain[r*2+1]*HP[1*4+c];
                        }
                    }
                    for (int t=0;t<16;++t) tracks_out[i].P[t] = P[t] - KHP[t];

                    // 记录
                    tracks_out[i].pos.push_back({{ z[0], z[1] }});
                    tracks_out[i].kf.push_back({{ tracks_out[i].x[0], tracks_out[i].x[1] }});
                    tracks_out[i].x_state.push_back(tracks_out[i].x);
                    tracks_out[i].time.push_back((k-1)*T);
                    if (k-1 < MT_dir.size() && best_j < static_cast<int>(MT_dir[k-1].size())) {
                        tracks_out[i].direction = MT_dir[k-1][best_j];
                    }
                    if (k-1 < MT_range.size() && best_j < static_cast<int>(MT_range[k-1].size())) {
                        tracks_out[i].range = MT_range[k-1][best_j];
                    }
                    tracks_out[i].last = k;
                    tracks_out[i].missed = 0;
                    assigned[best_j] = 1;
                    frameAccept += 1;
                    if (dbg && k <= dbgMaxFrames && static_cast<int>(i) < dbgMaxItems) {
                        std::cout << "[TRACK-KF] accept frame=" << (k-1)
                                  << " trIdx=" << i
                                  << " det=" << best_j
                                  << " z=(" << z[0] << "," << z[1] << ")"
                                  << " innov=(" << y0 << "," << y1 << ")"
                                  << " x_post=(" << tracks_out[i].x[0] << "," << tracks_out[i].x[1] << ","
                                  << tracks_out[i].x[2] << "," << tracks_out[i].x[3] << ")"
                                  << std::endl;
                    }
                } else {
                    // miss：不记录pos，只更新missed和last
                    // P矩阵的持续膨胀在预测步已正确处理：P = FPF' + Q
                    // 不在pos中记录预测值是为了确保隐含速度计算基于真实测量值
                    tracks_out[i].missed += 1;
                    tracks_out[i].last = k;
                    const bool failPos = !(best_d < sigma_thresh);
                    const bool failVel = !(implied_v < v_max);
                    if (failPos && failVel) frameRejectBoth += 1;
                    else if (failPos) frameRejectPos += 1;
                    else frameRejectVel += 1;
                    if (dbg && k <= dbgMaxFrames && static_cast<int>(i) < dbgMaxItems) {
                        std::cout << "[TRACK-KF] reject frame=" << (k-1)
                                  << " trIdx=" << i
                                  << " reason="
                                  << (failPos && failVel ? "pos_and_vel" : (failPos ? "pos" : "vel"))
                                  << " bestD=" << best_d
                                  << " gatePos=" << sigma_thresh
                                  << " impliedV=" << implied_v
                                  << " gateV=" << v_max
                                  << std::endl;
                    }
                }
            } else {
                // 理论上 M>0 时不会到这里，保留兜底
                tracks_out[i].missed += 1;
                tracks_out[i].last = k;
                frameNoDet += 1;
            }
        }

        // 3) 未分配检测 → 新建轨迹
        int newCount = 0;
        for (int j=0;j<M;++j){
            if (assigned[j]) continue;
            Track tr;
            tr.x = {{ det[j][0], det[j][1], 0.0, 0.0 }};
            for (int t=0;t<16;++t) tr.P[t]=0.0;
            tr.P[0]=tr.P[5]=tr.P[10]=tr.P[15]=10.0;
            tr.pos.push_back(det[j]);
            tr.kf.push_back(det[j]);
            tr.x_state.push_back(tr.x);
            tr.time.push_back((k-1)*T);
            if (k-1 < MT_dir.size() && j < static_cast<int>(MT_dir[k-1].size())) {
                tr.direction = MT_dir[k-1][j];
            }
            if (k-1 < MT_range.size() && j < static_cast<int>(MT_range[k-1].size())) {
                tr.range = MT_range[k-1][j];
            }
            tr.last = k;
            tr.missed = 0;
            tr.id = 0;
            tracks_out.push_back(tr);
            newCount += 1;

            if (dbg && k <= dbgMaxFrames && newCount <= dbgMaxItems) {
                std::cout << "[TRACK-KF] newTrack frame=" << (k-1)
                          << " det=" << j
                          << " xy=(" << det[j][0] << "," << det[j][1] << ")"
                          << std::endl;
            }
        }

        if (dbg && k <= dbgMaxFrames) {
            std::cout << "[TRACK-KF] frame=" << (k-1)
                      << " newTracks=" << newCount
                      << " accept=" << frameAccept
                      << " rejectPos=" << frameRejectPos
                      << " rejectVel=" << frameRejectVel
                      << " rejectBoth=" << frameRejectBoth
                      << " noDet=" << frameNoDet
                      << " singular=" << frameSingular
                      << " skipMaxGap=" << frameSkipMaxGap
                      << " totalTracksNow=" << tracks_out.size()
                      << std::endl;
        }
    }

    // ==== 末端筛选（长度 / 最大时间间隔 / 速度幅值）====
    std::vector<Track> kept;
    kept.reserve(tracks_out.size());
    for (size_t i=0;i<tracks_out.size();++i){
        const Track& tr = tracks_out[i];
        const size_t L = tr.pos.size();
        if (L < static_cast<size_t>(min_len)) continue;

        // 最大时间间隔
        double max_dt = 0.0;
        for (size_t t=1;t<tr.time.size();++t)
            max_dt = std::max(max_dt, tr.time[t]-tr.time[t-1]);
        if (max_dt > max_gap * T) continue;

        // 速度（按固定 T 近似，与 MATLAB 原逻辑一致）
        if (L >= 2){
            std::vector<double> vmag; vmag.reserve(L-1);
            for (size_t t=1;t<L;++t){
                const double dx = tr.pos[t][0] - tr.pos[t-1][0];
                const double dy = tr.pos[t][1] - tr.pos[t-1][1];
                const double v  = std::sqrt(dx*dx + dy*dy) / T;
                vmag.push_back(v);
            }
            bool too_fast = false;
            double mean=0.0;
            for (size_t t=0;t<vmag.size();++t){ mean += vmag[t]; if (vmag[t] > v_max) { too_fast=true; break; } }
            if (too_fast) continue;
            mean /= (vmag.empty()?1.0:double(vmag.size()));
            double var=0.0;
            for (size_t t=0;t<vmag.size();++t){ double d=vmag[t]-mean; var += d*d; }
            const double sd = (vmag.size()>1) ? std::sqrt(var/double(vmag.size()-1)) : 0.0;
            if (sd > sigma_thresh) continue;
        }

        kept.push_back(tr);
    }

    // 重新编号
    for (size_t i=0;i<kept.size();++i) kept[i].id = int(i+1);
    tracks_out.swap(kept);

    // --- 调试打印：输出所有通过筛选的航迹 ---
    if(dbg) {
        std::cout << "\n>>>>>> [GMTI Tracking Report] Total Tracks: " << tracks_out.size() << " <<<<<<" << std::endl;

        for (const auto& tr : tracks_out) {
            std::cout << "------------------------------------------------" << std::endl;
            std::cout << "Track ID: " << tr.id
                    << " | Points: " << tr.pos.size()
                    << " | Missed Frames: " << tr.missed << std::endl;

            // 打印该航迹经过的所有地理坐标（经纬度）
            std::cout << "Path (Lon, Lat): " << std::endl;
            for (size_t i = 0; i < tr.pos.size(); ++i) {
                double lat = 0.0, lng = 0.0;
                // 调用你代码中的高斯投影反算，将局部 x,y 转回经纬度
                Gaussp3RV(tr.pos[i][0], tr.pos[i][1], cfg.L0, lat, lng);

                printf("  [%02zu] (%.7f, %.7f) | Time: %.2f s\n", i, lng, lat, tr.time[i]);
            }

            // 计算当前速度 (基于最后两个点)
            if (tr.pos.size() >= 2) {
                double dx = tr.pos.back()[0] - tr.pos[tr.pos.size()-2][0];
                double dy = tr.pos.back()[1] - tr.pos[tr.pos.size()-2][1];
                double dt = tr.time.back() - tr.time[tr.time.size()-2];
                double speed = std::sqrt(dx*dx + dy*dy) / (dt > 0 ? dt : T);
                printf("Current Est. Speed: %.2f m/s (%.1f km/h)\n", speed, speed * 3.6);
            }
        }
        std::cout << ">>>>>> [End of Report] <<<<<<" << std::endl;
    }
    
    return true;
}
