#include "truth_writer.h"

#include <iomanip>
#include <sstream>

namespace gmti {
namespace target_injection {

bool TruthWriter::open(const std::string &truth_dir, std::string &err)
{
    if (!ensureDir(truth_dir)) {
        err = "failed to create truth dir: " + truth_dir;
        return false;
    }
    pulse_.open(joinPath(truth_dir, "truth_pulse.csv").c_str());
    summary_.open(joinPath(truth_dir, "truth_beam_summary.csv").c_str());
    if (!pulse_ || !summary_) {
        err = "failed to open truth output";
        return false;
    }
    pulse_ << "period_id,beam_id_0based,beam_id_1based,pulse_id,time_sec,target_id,target_name,"
           << "x_m,y_m,z_m,vx_mps,vy_mps,vz_mps,platform_x_m,platform_y_m,platform_z_m,"
           << "range_m,tau_abs_us,tau_rel_us,range_sample_float,range_sample_int,"
           << "target_azimuth_deg,theta_cmd_deg,theta_true_deg,angle_error_deg,beam_gain,"
           << "visible_by_beam,in_range_window,injection_enabled,local_background_rms,"
           << "target_amplitude,delta_phi_ch_rad,radial_velocity_mps\n";
    summary_ << "period_id,beam_id,target_id,visible_pulse_count,mean_range_m,"
             << "mean_range_sample,mean_target_azimuth_deg,mean_beam_gain,"
             << "mean_radial_velocity_mps,mean_target_amplitude,injected_sample_count\n";
    return true;
}

void TruthWriter::writePulse(const PulseTruth &t)
{
    const GeometrySample &g = t.geom;
    pulse_ << std::setprecision(12)
           << t.period_id << ","
           << t.beam_id << ","
           << (t.beam_id + 1) << ","
           << t.pulse_id << ","
           << g.time_sec << ","
           << t.target_id << ","
           << t.target_name << ","
           << g.target.position.x << ","
           << g.target.position.y << ","
           << g.target.position.z << ","
           << g.target.velocity.x << ","
           << g.target.velocity.y << ","
           << g.target.velocity.z << ","
           << g.platform.position.x << ","
           << g.platform.position.y << ","
           << g.platform.position.z << ","
           << g.range_m << ","
           << g.tau_abs_sec * 1.0e6 << ","
           << g.tau_rel_sec * 1.0e6 << ","
           << g.range_sample_float << ","
           << g.range_sample_int << ","
           << g.target_azimuth_deg << ","
           << g.theta_cmd_deg << ","
           << g.theta_true_deg << ","
           << g.angle_error_deg << ","
           << t.beam_gain << ","
           << (t.visible_by_beam ? 1 : 0) << ","
           << (g.in_range_window ? 1 : 0) << ","
           << (t.injection_enabled ? 1 : 0) << ","
           << t.local_background_rms << ","
           << t.target_amplitude << ","
           << t.delta_phi_ch_rad << ","
           << g.radial_velocity_mps << "\n";

    const std::pair<int, int> key(t.period_id, t.beam_id);
    BeamSummaryAccumulator &a = acc_[key];
    a.period_id = t.period_id;
    a.beam_id = t.beam_id;
    a.target_id = t.target_id;
    ++a.rows;
    if (t.visible_by_beam) ++a.visible_pulse_count;
    a.sum_range_m += g.range_m;
    a.sum_range_sample += g.range_sample_float;
    a.sum_target_azimuth_deg += g.target_azimuth_deg;
    a.sum_beam_gain += t.beam_gain;
    a.sum_radial_velocity_mps += g.radial_velocity_mps;
    a.sum_target_amplitude += t.target_amplitude;
    a.injected_sample_count += t.injected_sample_count;
}

void TruthWriter::writeSummary()
{
    for (std::map<std::pair<int, int>, BeamSummaryAccumulator>::const_iterator it = acc_.begin();
         it != acc_.end(); ++it) {
        const BeamSummaryAccumulator &a = it->second;
        const double denom = a.rows > 0 ? static_cast<double>(a.rows) : 1.0;
        summary_ << std::setprecision(12)
                 << a.period_id << ","
                 << a.beam_id << ","
                 << a.target_id << ","
                 << a.visible_pulse_count << ","
                 << a.sum_range_m / denom << ","
                 << a.sum_range_sample / denom << ","
                 << a.sum_target_azimuth_deg / denom << ","
                 << a.sum_beam_gain / denom << ","
                 << a.sum_radial_velocity_mps / denom << ","
                 << a.sum_target_amplitude / denom << ","
                 << a.injected_sample_count << "\n";
    }
}

void TruthWriter::close()
{
    if (pulse_) pulse_.close();
    if (summary_) summary_.close();
}

bool writeTargetInjectionReport(const std::string &path,
                                const InjectionConfig &cfg,
                                const InjectionStats &stats,
                                const std::string &notes,
                                std::string &err)
{
    std::ofstream out(path.c_str());
    if (!out) {
        err = "failed to write target injection report: " + path;
        return false;
    }
    out << "# 合作目标 DDC 域注入报告\n\n";
    out << "## 输入输出\n\n";
    out << "- input_config: `" << cfg.run.input_config << "`\n";
    out << "- input_data_dir: `" << cfg.run.input_data_dir << "`\n";
    out << "- output_dir: `" << cfg.run.output_dir << "`\n";
    out << "- output_data_file: `" << cfg.radar.output_data_file << "`\n\n";
    out << "## 新系统参数\n\n";
    out << "- fc: " << cfg.radar.fc_hz / 1.0e9 << " GHz\n";
    out << "- Br: " << cfg.radar.br_hz / 1.0e6 << " MHz\n";
    out << "- fs: " << cfg.radar.fs_hz / 1.0e6 << " MHz\n";
    out << "- Tr: " << cfg.radar.tr_sec * 1.0e6 << " us\n";
    out << "- PRF: " << cfg.radar.prf_hz << " Hz\n";
    out << "- DDC_len: " << cfg.radar.pulse_len << "\n";
    out << "- beam_count: " << cfg.radar.beam_count << "\n";
    out << "- pulse_num: " << cfg.radar.pulse_num << "\n";
    out << "- scan_min_deg: " << cfg.radar.scan_min_deg << "\n";
    out << "- scan_step_deg: " << cfg.radar.scan_step_deg << "\n";
    out << "- d_chan_m: " << cfg.radar.d_chan_m << "\n\n";
    out << "## 目标与注入模式\n\n";
    out << "- target_id: " << cfg.target.id << "\n";
    out << "- target_name: " << cfg.target.name << "\n";
    out << "- init_mode: " << cfg.target.init_mode << "\n";
    out << "- p0_m: (" << cfg.target.p0.x << ", " << cfg.target.p0.y << ", " << cfg.target.p0.z << ")\n";
    out << "- v_mps: (" << cfg.target.v.x << ", " << cfg.target.v.y << ", " << cfg.target.v.z << ")\n";
    out << "- visibility_mode: " << cfg.global.visibility_mode << "\n";
    out << "- amplitude_mode: " << cfg.global.amplitude_mode << "\n";
    out << "- target_snr_db: " << cfg.global.target_snr_db << "\n";
    out << "- direct_amplitude: " << cfg.global.direct_amplitude << "\n";
    out << "- chirp_phase_sign: " << cfg.global.chirp_phase_sign << "\n";
    out << "- carrier_phase_sign: " << cfg.global.carrier_phase_sign << "\n";
    out << "- channel_phase_mode: " << cfg.global.channel_phase_mode << "\n\n";
    out << "## 生成统计\n\n";
    out << "- packets_read: " << stats.packets_read << "\n";
    out << "- packets_written: " << stats.packets_written << "\n";
    out << "- pulses_injected: " << stats.pulses_injected << "\n";
    out << "- samples_injected: " << stats.samples_injected << "\n";
    out << "- max_amplitude: " << stats.max_amplitude << "\n";
    out << "- has_nan: " << (stats.has_nan ? "true" : "false") << "\n";
    out << "- has_inf: " << (stats.has_inf ? "true" : "false") << "\n\n";
    out << "## 当前说明\n\n";
    out << "- 合作目标注入位置是 DDC 原始快时间域 `range_sample_float`，不是脉压后 4096 点 bin。\n";
    out << "- 第一版 LFM 延迟约定采用 `0 <= t_fast - tau_rel < Tr`，后续需用单点脉压测试标定固定偏移和符号。\n";
    out << "- 第一版平台模型默认为 `ideal_platform`，真实 POS/导航接入留给下一步。\n";
    if (!notes.empty()) out << "- " << notes << "\n";
    return true;
}

} // namespace target_injection
} // namespace gmti

