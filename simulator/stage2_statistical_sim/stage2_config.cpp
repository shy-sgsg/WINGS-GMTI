#include "stage2_config.h"

#include "../target_injection/radar_geometry.h"
#include "tinyxml.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <set>
#include <sstream>

namespace gmti {
namespace stage2 {

namespace {

std::string readTextFile(const std::string &path)
{
    std::ifstream in(path.c_str());
    std::ostringstream ss;
    ss << in.rdbuf();
    return ss.str();
}

std::string sectionObject(const std::string &txt, const std::string &key)
{
    const std::string needle = "\"" + key + "\"";
    const size_t k = txt.find(needle);
    if (k == std::string::npos) return "";
    const size_t open = txt.find('{', k);
    if (open == std::string::npos) return "";
    int depth = 0;
    for (size_t i = open; i < txt.size(); ++i) {
        if (txt[i] == '{') ++depth;
        else if (txt[i] == '}') {
            --depth;
            if (depth == 0) return txt.substr(open, i - open + 1);
        }
    }
    return "";
}

std::vector<std::string> arrayObjects(const std::string &txt, const std::string &key)
{
    std::vector<std::string> out;
    const std::string needle = "\"" + key + "\"";
    const size_t k = txt.find(needle);
    if (k == std::string::npos) return out;
    const size_t open_arr = txt.find('[', k);
    if (open_arr == std::string::npos) return out;
    int arr_depth = 0;
    int obj_depth = 0;
    size_t obj_start = std::string::npos;
    for (size_t i = open_arr; i < txt.size(); ++i) {
        const char ch = txt[i];
        if (ch == '[') {
            ++arr_depth;
        } else if (ch == ']') {
            --arr_depth;
            if (arr_depth == 0) break;
        } else if (arr_depth > 0 && ch == '{') {
            if (obj_depth == 0) obj_start = i;
            ++obj_depth;
        } else if (arr_depth > 0 && ch == '}') {
            --obj_depth;
            if (obj_depth == 0 && obj_start != std::string::npos) {
                out.push_back(txt.substr(obj_start, i - obj_start + 1));
                obj_start = std::string::npos;
            }
        }
    }
    return out;
}

double jsonDouble(const std::string &obj, const std::string &key, double fallback)
{
    const std::string needle = "\"" + key + "\"";
    const size_t k = obj.find(needle);
    if (k == std::string::npos) return fallback;
    const size_t colon = obj.find(':', k);
    if (colon == std::string::npos) return fallback;
    return std::strtod(obj.c_str() + colon + 1, nullptr);
}

int jsonInt(const std::string &obj, const std::string &key, int fallback)
{
    return static_cast<int>(jsonDouble(obj, key, static_cast<double>(fallback)));
}

bool jsonBool(const std::string &obj, const std::string &key, bool fallback)
{
    const std::string needle = "\"" + key + "\"";
    const size_t k = obj.find(needle);
    if (k == std::string::npos) return fallback;
    const size_t colon = obj.find(':', k);
    if (colon == std::string::npos) return fallback;
    const std::string tail = obj.substr(colon + 1, 8);
    if (tail.find("true") != std::string::npos) return true;
    if (tail.find("false") != std::string::npos) return false;
    return fallback;
}

std::string jsonString(const std::string &obj, const std::string &key, const std::string &fallback)
{
    const std::string needle = "\"" + key + "\"";
    const size_t k = obj.find(needle);
    if (k == std::string::npos) return fallback;
    const size_t colon = obj.find(':', k);
    if (colon == std::string::npos) return fallback;
    const size_t q1 = obj.find('"', colon + 1);
    if (q1 == std::string::npos) return fallback;
    const size_t q2 = obj.find('"', q1 + 1);
    if (q2 == std::string::npos) return fallback;
    return obj.substr(q1 + 1, q2 - q1 - 1);
}

std::string valueArg(int &i, int argc, char **argv)
{
    if (i + 1 >= argc) return "";
    ++i;
    return argv[i];
}

int parseIntArg(const std::string &s, int fallback)
{
    if (s.empty()) return fallback;
    char *end = nullptr;
    const long v = std::strtol(s.c_str(), &end, 10);
    return end == s.c_str() ? fallback : static_cast<int>(v);
}

double parseDoubleArg(const std::string &s, double fallback)
{
    if (s.empty()) return fallback;
    char *end = nullptr;
    const double v = std::strtod(s.c_str(), &end);
    return end == s.c_str() ? fallback : v;
}

void setText(TiXmlElement *parent, const char *name, const std::string &value)
{
    TiXmlElement *e = new TiXmlElement(name);
    e->LinkEndChild(new TiXmlText(value.c_str()));
    parent->LinkEndChild(e);
}

void setDouble(TiXmlElement *parent, const char *name, double value)
{
    std::ostringstream ss;
    ss << value;
    setText(parent, name, ss.str());
}

void setInt(TiXmlElement *parent, const char *name, int value)
{
    std::ostringstream ss;
    ss << value;
    setText(parent, name, ss.str());
}

void setOrReplace(TiXmlElement *parent, const char *name, const std::string &value)
{
    TiXmlElement *e = parent->FirstChildElement(name);
    if (!e) {
        setText(parent, name, value);
        return;
    }
    e->Clear();
    e->LinkEndChild(new TiXmlText(value.c_str()));
}

void setOrReplaceDouble(TiXmlElement *parent, const char *name, double value)
{
    std::ostringstream ss;
    ss << value;
    setOrReplace(parent, name, ss.str());
}

void setOrReplaceInt(TiXmlElement *parent, const char *name, int value)
{
    std::ostringstream ss;
    ss << value;
    setOrReplace(parent, name, ss.str());
}

void writeGeometryXml(TiXmlElement *parent, const gmti::sim_geometry::Stage2GeometryConfig &g)
{
    TiXmlElement *geom = parent->FirstChildElement("simulation_geometry");
    if (!geom) {
        geom = new TiXmlElement("simulation_geometry");
        parent->LinkEndChild(geom);
    }
    setOrReplace(geom, "geometry_config_name", g.geometry_config_name);
    setOrReplace(geom, "local_x_axis", g.local_x_axis);
    setOrReplace(geom, "local_y_axis", g.local_y_axis);
    setOrReplace(geom, "platform_heading_source", g.platform_heading_source);
    setOrReplaceDouble(geom, "platform_heading_deg", g.platform_heading_deg);
    setOrReplace(geom, "beam_angle_reference", g.beam_angle_reference);
    setOrReplace(geom, "beam_zero_direction", g.beam_zero_direction);
    setOrReplace(geom, "beam_positive_direction", g.beam_positive_direction);
    setOrReplaceDouble(geom, "beam_theta_offset_deg", g.beam_theta_offset_deg);
    setOrReplace(geom, "range_geometry", g.range_geometry);
    setOrReplaceInt(geom, "use_ground_range_for_position", g.use_ground_range_for_position ? 1 : 0);
    setOrReplaceDouble(geom, "platform_origin_lat_deg", g.platform_origin_lat_deg);
    setOrReplaceDouble(geom, "platform_origin_lon_deg", g.platform_origin_lon_deg);
    setOrReplaceDouble(geom, "platform_origin_alt_m", g.platform_origin_alt_m);
    setOrReplaceDouble(geom, "projection_ref_lon_deg", g.projection_ref_lon_deg);
    setOrReplaceInt(geom, "squint_side", g.squint_side);
}

bool saveFromTemplate(const Stage2Config &cfg,
                      const std::string &out_xml,
                      const std::string &data_file)
{
    TiXmlDocument doc("outputs/stage1/config/temp_config_stage1_newsystem.xml");
    if (!doc.LoadFile()) return false;
    TiXmlElement *root = doc.RootElement();
    TiXmlElement *p = root ? root->FirstChildElement("GMTI_parameter") : nullptr;
    if (!p) return false;
    std::string result_dir = "outputs/stage2/algorithm_result";
    const std::string marker = "/data/";
    const size_t pos = data_file.find(marker);
    if (pos != std::string::npos) {
        result_dir = data_file.substr(0, pos) + "/algorithm_result";
    }
    setOrReplace(p, "GMTI_data_new", data_file);
    setOrReplace(p, "result_add", result_dir);
    setOrReplaceInt(p, "pulse_len", cfg.radar.pulse_len);
    setOrReplaceInt(p, "rg_len", cfg.radar.range_crop_len);
    setOrReplaceInt(p, "pulse_num", cfg.radar.pulse_num);
    setOrReplaceInt(p, "read_pulse_num", cfg.radar.pulse_num);
    setOrReplaceInt(p, "process_pulse_num", cfg.radar.pulse_num);
    setOrReplaceInt(p, "range_fft_len", cfg.radar.range_fft_len);
    setOrReplaceInt(p, "range_crop_start", cfg.radar.range_crop_start);
    setOrReplaceInt(p, "range_compress_len", cfg.radar.range_crop_len);
    setOrReplaceDouble(p, "fc", cfg.radar.fc_hz / 1.0e9);
    setOrReplaceDouble(p, "Br", cfg.radar.br_hz / 1.0e6);
    setOrReplaceDouble(p, "fs", cfg.radar.fs_hz / 1.0e6);
    setOrReplaceDouble(p, "Tr", cfg.radar.tr_sec * 1.0e6);
    setOrReplaceDouble(p, "PRF", cfg.radar.prf_hz);
    setOrReplaceDouble(p, "d_chan", cfg.radar.d_chan_m);
    setOrReplaceDouble(p, "ref_lon", cfg.projection_ref_lon_deg);
    setOrReplaceDouble(p, "scan_min_deg", cfg.radar.scan_min_deg);
    setOrReplaceDouble(p, "scan_step_deg", cfg.radar.scan_step_deg);
    setOrReplaceDouble(p, "scan_max_deg", cfg.radar.scan_min_deg + cfg.radar.scan_step_deg * (cfg.radar.beam_count - 1));
    setOrReplaceInt(p, "az_count", cfg.radar.beam_count);
    setOrReplaceInt(p, "wavepos_st", 1);
    setOrReplaceInt(p, "wavepos_ed", cfg.radar.beam_count);
    setOrReplaceInt(p, "wavepos_skip", 1);
    setOrReplaceDouble(p, "boshu", cfg.radar.beam_width_deg);
    setOrReplaceDouble(p, "beam_width_deg", cfg.radar.beam_width_deg);
    setOrReplaceDouble(p, "sample_delay_us", cfg.radar.sample_delay_sec * 1.0e6);
    setOrReplaceDouble(p, "sample_window_us", static_cast<double>(cfg.radar.pulse_len) / cfg.radar.fs_hz * 1.0e6);
    setOrReplace(p, "iq_data_type", cfg.radar.iq_data_type);
    setOrReplaceInt(p, "new_protocol_channel_count", cfg.radar.new_protocol_channel_count);
    setOrReplaceInt(p, "new_protocol_read_channel_1", cfg.radar.new_protocol_read_channel_1);
    setOrReplaceInt(p, "new_protocol_read_channel_2", cfg.radar.new_protocol_read_channel_2);
    setOrReplaceInt(p, "estimate_error_angle", 0);
    setOrReplaceDouble(p, "squint_angle", cfg.geometry.beam_theta_offset_deg);
    writeGeometryXml(p, cfg.geometry);
    const size_t data_pos = data_file.find(marker);
    if (data_pos != std::string::npos) {
        const std::string root_dir = data_file.substr(0, data_pos);
        setOrReplaceInt(p, "debug_pc_peak", 1);
        setOrReplace(p, "pc_peak_scene_truth", root_dir + "/truth/scene_truth.csv");
    }
    setOrReplaceDouble(p, "raw_fenbianlv", 25.0);
    setOrReplaceInt(p, "motion_comp_enable", 1);
    setOrReplaceInt(p, "motion_comp_analytic_enable", 1);
    setOrReplaceInt(p, "motion_comp_use_row_doppler", 1);
    setOrReplace(p, "motion_comp_solver", "debug");
    setOrReplaceInt(p, "motion_comp_iter", 8);
    setOrReplaceDouble(p, "motion_comp_iter_tol_mps", 1.0e-4);
    setOrReplaceInt(p, "ati_velocity_sign", 1);
    setOrReplaceInt(p, "ati_phase_to_velocity_sign", 1);
    setOrReplaceInt(p, "motion_doppler_axis_sign", cfg.radar.motion_doppler_axis_sign);
    setOrReplaceDouble(p, "ati_phase_bias_rad", 0.0);
    setOrReplaceDouble(p, "ati_vmax_mps", 1.6);
    setOrReplaceDouble(p, "motion_comp_denom_min", 1.0e-6);
    setOrReplaceDouble(p, "motion_comp_root_grid_step_mps", 0.02);
    setOrReplaceDouble(p, "motion_comp_root_cost_max", 0.25);
    setOrReplaceInt(p, "motion_comp_debug", 1);
    setOrReplaceInt(p, "p38_refit_enable", 1);
    setOrReplaceInt(p, "p38_refit_row_guard_bins", 2);
    setOrReplaceInt(p, "p38_refit_range_guard_bins", 2);
    setOrReplaceDouble(p, "p38_refit_top_power_frac", 0.01);
    setOrReplaceInt(p, "p38_refit_min_sample_count", 8);
    setOrReplaceDouble(p, "p38_refit_min_inlier_ratio", 0.60);
    setOrReplaceDouble(p, "p38_refit_max_rmse_rad", 0.60);
    setOrReplaceDouble(p, "p38_refit_max_delta_k", 0.01);
    setOrReplaceDouble(p, "p38_refit_max_delta_b_rad", 1.50);
    return doc.SaveFile(out_xml.c_str());
}

} // namespace

bool parseStage2CommandLine(int argc, char **argv, Stage2RunOptions &opt)
{
    for (int i = 1; i < argc; ++i) {
        const std::string k = argv[i];
        const std::string v = (i + 1 < argc) ? argv[i + 1] : "";
        if (k == "--config") opt.run_config = valueArg(i, argc, argv);
        else if (!k.empty() && k[0] != '-' && opt.run_config.empty()) opt.run_config = k;
        else if (k == "--stage2-config") { opt.stage2_config = valueArg(i, argc, argv); opt.legacy_override_used = true; }
        else if (k == "--target-config") { opt.target_config = valueArg(i, argc, argv); opt.legacy_override_used = true; }
        else if (k == "--output-dir") { opt.output_dir = valueArg(i, argc, argv); opt.legacy_override_used = true; }
        else if (k == "--scene-mode") { opt.scene_mode = valueArg(i, argc, argv); opt.legacy_override_used = true; }
        else if (k == "--target-enabled") { const std::string bv = valueArg(i, argc, argv); opt.target_enabled = (bv != "false" && bv != "0"); opt.legacy_override_used = true; }
        else if (k == "--period-start") { opt.period_start = parseIntArg(valueArg(i, argc, argv), opt.period_start); opt.legacy_override_used = true; }
        else if (k == "--period-count") { opt.period_count = parseIntArg(valueArg(i, argc, argv), opt.period_count); opt.legacy_override_used = true; }
        else if (k == "--single-scatterer-range") { opt.single_scatterer_range_m = parseDoubleArg(valueArg(i, argc, argv), opt.single_scatterer_range_m); opt.legacy_override_used = true; }
        else if (k == "--single-point-beam-id") { opt.single_point_beam_id_1based = parseIntArg(valueArg(i, argc, argv), opt.single_point_beam_id_1based); opt.legacy_override_used = true; }
        else if (k == "--single-point-expected-bin") { opt.single_point_expected_bin = parseIntArg(valueArg(i, argc, argv), opt.single_point_expected_bin); opt.legacy_override_used = true; }
        else if (k == "--single-scatterer-azimuth") { opt.single_scatterer_azimuth_deg = parseDoubleArg(valueArg(i, argc, argv), opt.single_scatterer_azimuth_deg); opt.legacy_override_used = true; }
        else if (k == "--single-scatterer-amplitude") { opt.single_scatterer_amplitude = parseDoubleArg(valueArg(i, argc, argv), opt.single_scatterer_amplitude); opt.legacy_override_used = true; }
        else if (k == "--moving-target-beam-id") { opt.moving_target_beam_id_1based = parseIntArg(valueArg(i, argc, argv), opt.moving_target_beam_id_1based); opt.legacy_override_used = true; }
        else if (k == "--moving-target-expected-bin") { opt.moving_target_expected_bin = parseIntArg(valueArg(i, argc, argv), opt.moving_target_expected_bin); opt.legacy_override_used = true; }
        else if (k == "--moving-target-speed-mps") { opt.moving_target_speed_mps = parseDoubleArg(valueArg(i, argc, argv), opt.moving_target_speed_mps); opt.legacy_override_used = true; }
        else if (k == "--moving-target-velocity-mode") { opt.moving_target_velocity_mode = valueArg(i, argc, argv); opt.legacy_override_used = true; }
        else if (k == "--moving-target-ve-mps") { opt.moving_target_ve_mps = parseDoubleArg(valueArg(i, argc, argv), opt.moving_target_ve_mps); opt.legacy_override_used = true; }
        else if (k == "--moving-target-vn-mps") { opt.moving_target_vn_mps = parseDoubleArg(valueArg(i, argc, argv), opt.moving_target_vn_mps); opt.legacy_override_used = true; }
        else if (k == "--moving-target-radial-speed-mps") { opt.moving_target_radial_speed_mps = parseDoubleArg(valueArg(i, argc, argv), opt.moving_target_radial_speed_mps); opt.legacy_override_used = true; }
        else if (k == "--moving-target-tangential-speed-mps") { opt.moving_target_tangential_speed_mps = parseDoubleArg(valueArg(i, argc, argv), opt.moving_target_tangential_speed_mps); opt.legacy_override_used = true; }
        else if (k == "--moving-target-rcs-db") { opt.moving_target_rcs_db = parseDoubleArg(valueArg(i, argc, argv), opt.moving_target_rcs_db); opt.legacy_override_used = true; }
        else if (k == "--clutter-amplitude-scale") { opt.clutter_amplitude_scale = parseDoubleArg(valueArg(i, argc, argv), opt.clutter_amplitude_scale); opt.legacy_override_used = true; }
        else if (k == "--area-clutter-scatterer-count") { opt.area_clutter_scatterer_count = parseIntArg(valueArg(i, argc, argv), opt.area_clutter_scatterer_count); opt.legacy_override_used = true; }
        else if (k == "--area-clutter-mean-power") { opt.area_clutter_mean_power = parseDoubleArg(valueArg(i, argc, argv), opt.area_clutter_mean_power); opt.legacy_override_used = true; }
        else if (k == "--area-clutter-texture-sigma") { opt.area_clutter_texture_sigma = parseDoubleArg(valueArg(i, argc, argv), opt.area_clutter_texture_sigma); opt.legacy_override_used = true; }
        else if (k == "--strong-scatterer-count") { opt.strong_scatterer_count = parseIntArg(valueArg(i, argc, argv), opt.strong_scatterer_count); opt.legacy_override_used = true; }
        else if (k == "--strong-rcs-db-min") { opt.strong_rcs_db_min = parseDoubleArg(valueArg(i, argc, argv), opt.strong_rcs_db_min); opt.legacy_override_used = true; }
        else if (k == "--strong-rcs-db-max") { opt.strong_rcs_db_max = parseDoubleArg(valueArg(i, argc, argv), opt.strong_rcs_db_max); opt.legacy_override_used = true; }
        else if (k == "--line-scatterer-count") { opt.line_scatterer_count = parseIntArg(valueArg(i, argc, argv), opt.line_scatterer_count); opt.legacy_override_used = true; }
        else if (k == "--line-points-per-line") { opt.line_points_per_line = parseIntArg(valueArg(i, argc, argv), opt.line_points_per_line); opt.legacy_override_used = true; }
        else if (k == "--line-rcs-db") { opt.line_rcs_db = parseDoubleArg(valueArg(i, argc, argv), opt.line_rcs_db); opt.legacy_override_used = true; }
        else if (k == "--noise-power") { opt.noise_power = parseDoubleArg(valueArg(i, argc, argv), opt.noise_power); opt.legacy_override_used = true; }
        else if (k == "--validate") { const std::string bv = valueArg(i, argc, argv); opt.validate = (bv != "false" && bv != "0"); opt.legacy_override_used = true; }
    }
    if (!opt.run_config.empty() && opt.legacy_override_used) {
        opt.parse_error = "--config/run JSON cannot be mixed with legacy stage2 CLI overrides";
        return false;
    }
    return true;
}

bool loadStage2Config(const std::string &path, Stage2Config &cfg, std::string &err)
{
    std::ifstream test(path.c_str());
    if (!test) {
        if (!writeDefaultStage2Config(path, err)) return false;
    }
    const std::string txt = readTextFile(path);
    const std::string sys = sectionObject(txt, "system");
    const std::string sim = sectionObject(txt, "simulation");
    const std::string scene = sectionObject(txt, "scene");
    const std::string targets = sectionObject(txt, "targets");
    if (sys.empty() || sim.empty() || scene.empty()) {
        err = "stage2_config.json must contain system, simulation, and scene";
        return false;
    }
    cfg.radar.fc_hz = jsonDouble(sys, "fc_ghz", 16.0) * 1.0e9;
    cfg.radar.br_hz = jsonDouble(sys, "bandwidth_mhz", 50.0) * 1.0e6;
    cfg.radar.fs_hz = jsonDouble(sys, "fs_mhz", 60.0) * 1.0e6;
    cfg.radar.tr_sec = jsonDouble(sys, "pulse_width_us", 130.0) * 1.0e-6;
    cfg.radar.prf_hz = jsonDouble(sys, "prf_hz", 1300.0);
    cfg.radar.sample_delay_sec = jsonDouble(sys, "sample_delay_us", 488.0) * 1.0e-6;
    cfg.radar.pulse_len = jsonInt(sys, "ddc_len", 11820);
    cfg.radar.range_fft_len = jsonInt(sys, "fft_len", 12288);
    cfg.radar.range_crop_start = jsonInt(sys, "pc_crop_start", 3864);
    cfg.radar.range_crop_len = jsonInt(sys, "pc_crop_len", 4096);
    cfg.radar.scan_min_deg = jsonDouble(sys, "scan_min_deg", -60.0);
    cfg.radar.scan_step_deg = jsonDouble(sys, "scan_step_deg", 2.0);
    cfg.radar.beam_count = jsonInt(sys, "beam_count", 61);
    cfg.radar.iq_data_type = jsonString(sys, "iq_data_type", cfg.radar.iq_data_type);
    cfg.radar.new_protocol_channel_count = jsonInt(sys, "new_protocol_channel_count", cfg.radar.new_protocol_channel_count);
    cfg.radar.new_protocol_read_channel_1 = jsonInt(sys, "new_protocol_read_channel_1", cfg.radar.new_protocol_read_channel_1);
    cfg.radar.new_protocol_read_channel_2 = jsonInt(sys, "new_protocol_read_channel_2", cfg.radar.new_protocol_read_channel_2);
    cfg.radar.motion_doppler_axis_sign = jsonInt(sys, "motion_doppler_axis_sign", cfg.radar.motion_doppler_axis_sign);
    cfg.radar.beam_width_deg = jsonDouble(sys, "beam_width_deg", 2.28);
    cfg.radar.pulse_num = jsonInt(sys, "pulse_num", 130);
    cfg.radar.d_chan_m = jsonDouble(sys, "d_chan_m", 0.17);
    cfg.platform_height_m = jsonDouble(sys, "platform_height_m", cfg.platform_height_m);
    cfg.platform_speed_mps = jsonDouble(sys, "platform_speed_mps", cfg.platform_speed_mps);
    cfg.platform_origin_lat_deg = jsonDouble(sys, "platform_origin_lat_deg", cfg.platform_origin_lat_deg);
    cfg.platform_origin_lon_deg = jsonDouble(sys, "platform_origin_lon_deg", cfg.platform_origin_lon_deg);
    cfg.platform_origin_alt_m = jsonDouble(sys, "platform_origin_alt_m", cfg.platform_origin_alt_m);
    cfg.projection_ref_lon_deg = jsonDouble(sys, "projection_ref_lon_deg", cfg.projection_ref_lon_deg);
    cfg.geometry.platform_origin_lat_deg = cfg.platform_origin_lat_deg;
    cfg.geometry.platform_origin_lon_deg = cfg.platform_origin_lon_deg;
    cfg.geometry.platform_origin_alt_m = cfg.platform_origin_alt_m;
    cfg.geometry.projection_ref_lon_deg = cfg.projection_ref_lon_deg;

    const std::string geom = sectionObject(txt, "simulation_geometry");
    if (!geom.empty()) {
        cfg.geometry.geometry_config_name = jsonString(geom, "geometry_config_name", cfg.geometry.geometry_config_name);
        cfg.geometry.local_x_axis = jsonString(geom, "local_x_axis", cfg.geometry.local_x_axis);
        cfg.geometry.local_y_axis = jsonString(geom, "local_y_axis", cfg.geometry.local_y_axis);
        cfg.geometry.platform_heading_source = jsonString(geom, "platform_heading_source", cfg.geometry.platform_heading_source);
        cfg.geometry.platform_heading_deg = jsonDouble(geom, "platform_heading_deg", cfg.geometry.platform_heading_deg);
        cfg.geometry.beam_angle_reference = jsonString(geom, "beam_angle_reference", cfg.geometry.beam_angle_reference);
        cfg.geometry.beam_zero_direction = jsonString(geom, "beam_zero_direction", cfg.geometry.beam_zero_direction);
        cfg.geometry.beam_positive_direction = jsonString(geom, "beam_positive_direction", cfg.geometry.beam_positive_direction);
        cfg.geometry.beam_theta_offset_deg = jsonDouble(geom, "beam_theta_offset_deg", cfg.geometry.beam_theta_offset_deg);
        cfg.geometry.range_geometry = jsonString(geom, "range_geometry", cfg.geometry.range_geometry);
        cfg.geometry.use_ground_range_for_position = jsonBool(geom, "use_ground_range_for_position", cfg.geometry.use_ground_range_for_position);
        cfg.geometry.platform_origin_lat_deg = jsonDouble(geom, "platform_origin_lat_deg", cfg.geometry.platform_origin_lat_deg);
        cfg.geometry.platform_origin_lon_deg = jsonDouble(geom, "platform_origin_lon_deg", cfg.geometry.platform_origin_lon_deg);
        cfg.geometry.platform_origin_alt_m = jsonDouble(geom, "platform_origin_alt_m", cfg.geometry.platform_origin_alt_m);
        cfg.geometry.projection_ref_lon_deg = jsonDouble(geom, "projection_ref_lon_deg", cfg.geometry.projection_ref_lon_deg);
        cfg.geometry.squint_side = jsonInt(geom, "squint_side", cfg.geometry.squint_side);
        cfg.platform_origin_lat_deg = cfg.geometry.platform_origin_lat_deg;
        cfg.platform_origin_lon_deg = cfg.geometry.platform_origin_lon_deg;
        cfg.platform_origin_alt_m = cfg.geometry.platform_origin_alt_m;
        cfg.projection_ref_lon_deg = cfg.geometry.projection_ref_lon_deg;
    }

    cfg.sim.period_start = jsonInt(sim, "period_start", cfg.sim.period_start);
    cfg.sim.period_count = jsonInt(sim, "period_count", cfg.sim.period_count);
    cfg.sim.random_seed = static_cast<uint32_t>(jsonInt(sim, "random_seed", static_cast<int>(cfg.sim.random_seed)));
    cfg.sim.platform_mode = jsonString(sim, "platform_mode", cfg.sim.platform_mode);
    cfg.sim.beam_pattern_mode = jsonString(sim, "beam_pattern_mode", cfg.sim.beam_pattern_mode);
    cfg.sim.channel_phase_mode = jsonString(sim, "channel_phase_mode", cfg.sim.channel_phase_mode);
    cfg.sim.chirp_phase_sign = jsonInt(sim, "chirp_phase_sign", cfg.sim.chirp_phase_sign);
    cfg.sim.carrier_phase_sign = jsonInt(sim, "carrier_phase_sign", cfg.sim.carrier_phase_sign);

    cfg.scene.range_min_m = jsonDouble(scene, "range_min_m", cfg.scene.range_min_m);
    cfg.scene.range_max_m = jsonDouble(scene, "range_max_m", cfg.scene.range_max_m);
    cfg.scene.azimuth_min_deg = jsonDouble(scene, "azimuth_min_deg", cfg.scene.azimuth_min_deg);
    cfg.scene.azimuth_max_deg = jsonDouble(scene, "azimuth_max_deg", cfg.scene.azimuth_max_deg);
    cfg.scene.ground_z_m = jsonDouble(scene, "ground_z_m", cfg.scene.ground_z_m);
    cfg.scene.clutter_amplitude_scale = jsonDouble(scene, "clutter_amplitude_scale", cfg.scene.clutter_amplitude_scale);
    const std::string area = sectionObject(scene, "area_clutter");
    cfg.scene.area.enabled = jsonBool(area, "enabled", cfg.scene.area.enabled);
    cfg.scene.area.model = jsonString(area, "model", cfg.scene.area.model);
    cfg.scene.area.scatterer_count = jsonInt(area, "scatterer_count", cfg.scene.area.scatterer_count);
    cfg.scene.area.mean_power = jsonDouble(area, "mean_power", cfg.scene.area.mean_power);
    cfg.scene.area.texture_sigma = jsonDouble(area, "texture_sigma", cfg.scene.area.texture_sigma);
    const std::string strong = sectionObject(scene, "strong_scatterers");
    cfg.scene.strong.enabled = jsonBool(strong, "enabled", cfg.scene.strong.enabled);
    cfg.scene.strong.count = jsonInt(strong, "count", cfg.scene.strong.count);
    cfg.scene.strong.rcs_db_min = jsonDouble(strong, "rcs_db_min", cfg.scene.strong.rcs_db_min);
    cfg.scene.strong.rcs_db_max = jsonDouble(strong, "rcs_db_max", cfg.scene.strong.rcs_db_max);
    const std::string line = sectionObject(scene, "line_scatterers");
    cfg.scene.line.enabled = jsonBool(line, "enabled", cfg.scene.line.enabled);
    cfg.scene.line.line_count = jsonInt(line, "line_count", cfg.scene.line.line_count);
    cfg.scene.line.points_per_line = jsonInt(line, "points_per_line", cfg.scene.line.points_per_line);
    cfg.scene.line.rcs_db = jsonDouble(line, "rcs_db", cfg.scene.line.rcs_db);
    const std::string noise = sectionObject(scene, "thermal_noise");
    cfg.scene.noise.enabled = jsonBool(noise, "enabled", cfg.scene.noise.enabled);
    cfg.scene.noise.noise_power = jsonDouble(noise, "noise_power", cfg.scene.noise.noise_power);
    if (!targets.empty()) {
        cfg.target.enabled = jsonBool(targets, "enabled", cfg.target.enabled);
        cfg.target.target_config_path = jsonString(targets, "target_config_path", cfg.target.target_config_path);
        cfg.target.target_snr_db = jsonDouble(targets, "target_snr_db", cfg.target.target_snr_db);
        cfg.target.amplitude_mode = jsonString(targets, "amplitude_mode", cfg.target.amplitude_mode);
    }
    return true;
}

bool writeDefaultStage2Config(const std::string &path, std::string &err)
{
    std::ofstream out(path.c_str());
    if (!out) {
        err = "failed to create stage2 config: " + path;
        return false;
    }
    out <<
        "{\n"
        "  \"system\": {\n"
        "    \"fc_ghz\": 16.0, \"bandwidth_mhz\": 50.0, \"fs_mhz\": 60.0,\n"
        "    \"pulse_width_us\": 130.0, \"prf_hz\": 1300.0,\n"
        "    \"sample_delay_us\": 488.0, \"sample_window_us\": 197.0,\n"
        "    \"ddc_len\": 11820, \"fft_len\": 12288, \"pc_crop_start\": 3864, \"pc_crop_len\": 4096,\n"
        "    \"scan_min_deg\": -60.0, \"scan_step_deg\": 2.0, \"beam_count\": 61,\n"
        "    \"iq_data_type\": \"float32\", \"new_protocol_channel_count\": 2,\n"
        "    \"new_protocol_read_channel_1\": 1, \"new_protocol_read_channel_2\": 2,\n"
        "    \"beam_width_deg\": 2.28, \"pulse_num\": 130,\n"
        "    \"platform_height_m\": 6000.0, \"platform_speed_mps\": 60.0, \"d_chan_m\": 0.17,\n"
        "    \"platform_origin_lat_deg\": 40.45121057, \"platform_origin_lon_deg\": 116.98377429,\n"
        "    \"platform_origin_alt_m\": 6000.0, \"projection_ref_lon_deg\": 117.0\n"
        "  },\n"
        "  \"simulation_geometry\": {\n"
        "    \"geometry_config_name\": \"algorithm_axis_x_north\",\n"
        "    \"local_x_axis\": \"north\", \"local_y_axis\": \"east\",\n"
        "    \"platform_heading_source\": \"velocity\", \"platform_heading_deg\": 0.0,\n"
        "    \"beam_angle_reference\": \"algorithm\",\n"
        "    \"beam_zero_direction\": \"algorithm\",\n"
        "    \"beam_positive_direction\": \"algorithm\",\n"
        "    \"beam_theta_offset_deg\": 0.0,\n"
        "    \"range_geometry\": \"algorithm\",\n"
        "    \"use_ground_range_for_position\": true,\n"
        "    \"platform_origin_lat_deg\": 40.4512107203,\n"
        "    \"platform_origin_lon_deg\": 116.985931582,\n"
        "    \"platform_origin_alt_m\": 6000.0,\n"
        "    \"projection_ref_lon_deg\": 117.0,\n"
        "    \"squint_side\": 1\n"
        "  },\n"
        "  \"simulation\": {\n"
        "    \"period_start\": 0, \"period_count\": 1, \"random_seed\": 202606,\n"
        "    \"platform_mode\": \"ideal_straight\", \"beam_pattern_mode\": \"gaussian\",\n"
        "    \"channel_phase_mode\": \"baseline_approx\", \"chirp_phase_sign\": 1,\n"
        "    \"carrier_phase_sign\": -1, \"output_format\": \"compatible_with_algorithm\"\n"
        "  },\n"
        "  \"scene\": {\n"
        "    \"range_min_m\": 82800.0, \"range_max_m\": 93000.0,\n"
        "    \"azimuth_min_deg\": -60.0, \"azimuth_max_deg\": 60.0, \"ground_z_m\": 0.0,\n"
        "    \"clutter_amplitude_scale\": 1.0,\n"
        "    \"area_clutter\": {\"enabled\": true, \"model\": \"rayleigh_lognormal_texture\", \"scatterer_count\": 1000, \"mean_power\": 1.0, \"texture_sigma\": 0.4, \"spatial_cell_m\": 30.0},\n"
        "    \"strong_scatterers\": {\"enabled\": true, \"count\": 20, \"rcs_db_min\": 10.0, \"rcs_db_max\": 30.0},\n"
        "    \"line_scatterers\": {\"enabled\": true, \"line_count\": 1, \"points_per_line\": 50, \"rcs_db\": 12.0},\n"
        "    \"thermal_noise\": {\"enabled\": true, \"noise_power\": 0.01}\n"
        "  },\n"
        "  \"targets\": {\"enabled\": true, \"target_config_path\": \"targets.json\", \"target_snr_db\": 25.0, \"amplitude_mode\": \"snr_db\"}\n"
        "}\n";
    return true;
}

bool writeStage2OutputConfig(const Stage2Config &cfg,
                             const std::string &out_xml,
                             const std::string &data_file,
                             std::string &err)
{
    if (saveFromTemplate(cfg, out_xml, data_file)) {
        return true;
    }
    TiXmlDocument doc;
    TiXmlDeclaration *decl = new TiXmlDeclaration("1.0", "UTF-8", "");
    doc.LinkEndChild(decl);
    TiXmlElement *root = new TiXmlElement("GMTI");
    TiXmlElement *p = new TiXmlElement("GMTI_parameter");
    root->LinkEndChild(p);
    doc.LinkEndChild(root);
    setText(p, "GMTI_data_new", data_file);
    setInt(p, "pulse_len", cfg.radar.pulse_len);
    setInt(p, "pulse_num", cfg.radar.pulse_num);
    setInt(p, "read_pulse_num", cfg.radar.pulse_num);
    setInt(p, "process_pulse_num", cfg.radar.pulse_num);
    setDouble(p, "fc", cfg.radar.fc_hz / 1.0e9);
    setDouble(p, "Br", cfg.radar.br_hz / 1.0e6);
    setDouble(p, "fs", cfg.radar.fs_hz / 1.0e6);
    setDouble(p, "Tr", cfg.radar.tr_sec * 1.0e6);
    setDouble(p, "PRF", cfg.radar.prf_hz);
    setDouble(p, "d_chan", cfg.radar.d_chan_m);
    setDouble(p, "ref_lon", cfg.projection_ref_lon_deg);
    setDouble(p, "scan_min_deg", cfg.radar.scan_min_deg);
    setDouble(p, "scan_step_deg", cfg.radar.scan_step_deg);
    setDouble(p, "scan_max_deg", cfg.radar.scan_min_deg + cfg.radar.scan_step_deg * (cfg.radar.beam_count - 1));
    setInt(p, "az_count", cfg.radar.beam_count);
    setDouble(p, "beam_width_deg", cfg.radar.beam_width_deg);
    setInt(p, "range_fft_len", cfg.radar.range_fft_len);
    setInt(p, "range_crop_start", cfg.radar.range_crop_start);
    setInt(p, "range_compress_len", cfg.radar.range_crop_len);
    setDouble(p, "sample_delay_us", cfg.radar.sample_delay_sec * 1.0e6);
    setDouble(p, "sample_window_us", static_cast<double>(cfg.radar.pulse_len) / cfg.radar.fs_hz * 1.0e6);
    setText(p, "iq_data_type", cfg.radar.iq_data_type);
    setInt(p, "new_protocol_channel_count", cfg.radar.new_protocol_channel_count);
    setInt(p, "new_protocol_read_channel_1", cfg.radar.new_protocol_read_channel_1);
    setInt(p, "new_protocol_read_channel_2", cfg.radar.new_protocol_read_channel_2);
    setInt(p, "estimate_error_angle", 0);
    setDouble(p, "squint_angle", cfg.geometry.beam_theta_offset_deg);
    setDouble(p, "raw_fenbianlv", 25.0);
    setInt(p, "motion_comp_enable", 1);
    setInt(p, "motion_comp_analytic_enable", 1);
    setInt(p, "motion_comp_use_row_doppler", 1);
    setOrReplace(p, "motion_comp_solver", "debug");
    setInt(p, "motion_comp_iter", 8);
    setDouble(p, "motion_comp_iter_tol_mps", 1.0e-4);
    setInt(p, "ati_velocity_sign", 1);
    setInt(p, "ati_phase_to_velocity_sign", 1);
    setInt(p, "motion_doppler_axis_sign", cfg.radar.motion_doppler_axis_sign);
    setDouble(p, "ati_phase_bias_rad", 0.0);
    setDouble(p, "ati_vmax_mps", 1.6);
    setDouble(p, "motion_comp_denom_min", 1.0e-6);
    setDouble(p, "motion_comp_root_grid_step_mps", 0.02);
    setDouble(p, "motion_comp_root_cost_max", 0.25);
    setInt(p, "motion_comp_debug", 1);
    setInt(p, "p38_refit_enable", 1);
    setInt(p, "p38_refit_row_guard_bins", 2);
    setInt(p, "p38_refit_range_guard_bins", 2);
    setDouble(p, "p38_refit_top_power_frac", 0.01);
    setInt(p, "p38_refit_min_sample_count", 8);
    setDouble(p, "p38_refit_min_inlier_ratio", 0.60);
    setDouble(p, "p38_refit_max_rmse_rad", 0.60);
    setDouble(p, "p38_refit_max_delta_k", 0.01);
    setDouble(p, "p38_refit_max_delta_b_rad", 1.50);
    writeGeometryXml(p, cfg.geometry);
    if (!doc.SaveFile(out_xml.c_str())) {
        err = "failed to save stage2 output xml";
        return false;
    }
    return true;
}

gmti::target_injection::TargetGlobalConfig makeTargetGlobal(const Stage2Config &cfg)
{
    gmti::target_injection::TargetGlobalConfig g;
    g.visibility_mode = cfg.sim.beam_pattern_mode == "hard_gate" ? "hard_gate" : "gaussian";
    g.amplitude_mode = cfg.target.amplitude_mode;
    g.channel_phase_mode = cfg.sim.channel_phase_mode;
    g.target_snr_db = cfg.target.target_snr_db;
    g.chirp_phase_sign = cfg.sim.chirp_phase_sign;
    g.carrier_phase_sign = cfg.sim.carrier_phase_sign;
    g.platform_speed_mps = cfg.platform_speed_mps;
    g.platform_height_m = cfg.platform_height_m;
    g.platform_origin_lat_deg = cfg.platform_origin_lat_deg;
    g.platform_origin_lon_deg = cfg.platform_origin_lon_deg;
    g.platform_origin_alt_m = cfg.platform_origin_alt_m;
    g.projection_ref_lon_deg = cfg.projection_ref_lon_deg;
    g.geometry = cfg.geometry;
    return g;
}

namespace {

std::string jsonType(const std::string &obj, const std::string &fallback)
{
    return jsonString(obj, "type", fallback);
}

bool resolveRunTarget(Stage2RunConfig &run, Stage2RunTarget &rt, std::string &err)
{
    const Stage2Config &cfg = run.cfg;
    const int beam0 = rt.beam_id - run.beam_index_base;
    if (beam0 < 0 || beam0 >= cfg.radar.beam_count) {
        err = "target " + rt.target_id + " has invalid beam_id";
        return false;
    }
    rt.theta_cmd_deg = cfg.radar.scan_min_deg + static_cast<double>(beam0) * cfg.radar.scan_step_deg;
    if (rt.init_type == "beam_bin_azimuth_offset") {
        rt.azimuth_deg = rt.theta_cmd_deg + rt.azimuth_offset_deg;
    }
    const double range_sample_float = static_cast<double>(cfg.radar.range_crop_start + rt.expected_bin);
    const int range_sample_int = static_cast<int>(std::floor(range_sample_float + 0.5));
    const double range_m = 0.5 * gmti::target_injection::kC *
        (cfg.radar.sample_delay_sec + range_sample_float / cfg.radar.fs_hz);
    const int ref_pulse_idx = std::max(0, cfg.radar.pulse_num / 2);
    const double ref_time_s =
        gmti::target_injection::pulseTimeSec(cfg.radar, cfg.sim.period_start, beam0, ref_pulse_idx);
    const gmti::target_injection::PlatformState ref_platform =
        gmti::target_injection::evaluatePlatformState(run.global, ref_time_s);
    const gmti::sim_geometry::LocalPoint target_ref_local =
        gmti::sim_geometry::makePointFromRangeAzimuth(
            gmti::sim_geometry::LocalPoint(ref_platform.position.x,
                                           ref_platform.position.y,
                                           ref_platform.position.z),
            gmti::sim_geometry::LocalVelocity(ref_platform.velocity.x,
                                              ref_platform.velocity.y,
                                              ref_platform.velocity.z),
            range_m,
            rt.azimuth_deg,
            ref_platform.position.z,
            0.0,
            cfg.geometry);
    const gmti::target_injection::Vec3 target_ref(target_ref_local.x, target_ref_local.y, 0.0);
    const gmti::sim_geometry::ENUPoint ref_enu =
        gmti::sim_geometry::localToEnu(
            gmti::sim_geometry::LocalPoint(ref_platform.position.x,
                                           ref_platform.position.y,
                                           ref_platform.position.z),
            cfg.geometry);
    const gmti::sim_geometry::ENUPoint target_enu =
        gmti::sim_geometry::localToEnu(target_ref_local, cfg.geometry);
    const double de = target_enu.e - ref_enu.e;
    const double dn = target_enu.n - ref_enu.n;
    const double ground_range = std::sqrt(de * de + dn * dn);
    const gmti::sim_geometry::ENUVelocity ref_vel_en =
        gmti::sim_geometry::localVelocityToEnu(
            gmti::sim_geometry::LocalVelocity(ref_platform.velocity.x,
                                              ref_platform.velocity.y,
                                              ref_platform.velocity.z),
            cfg.geometry);
    const gmti::sim_geometry::LookVectorEN look =
        gmti::sim_geometry::makeAlgorithmLookVectorEN(ref_vel_en.ve,
                                                      ref_vel_en.vn,
                                                      rt.azimuth_deg,
                                                      cfg.geometry);
    const double ur_e = ground_range > 1.0e-9 ? de / ground_range : look.east;
    const double ur_n = ground_range > 1.0e-9 ? dn / ground_range : look.north;
    const double ut_e = -ur_n;
    const double ut_n = ur_e;
    const double target_vr_self = rt.ve_mps * ur_e + rt.vn_mps * ur_n;
    const double target_vt_self = rt.ve_mps * ut_e + rt.vn_mps * ut_n;
    const double lambda = gmti::target_injection::kC / cfg.radar.fc_hz;
    const double af_motion_truth = lambda > 0.0
        ? static_cast<double>(cfg.radar.motion_doppler_axis_sign) * 2.0 * target_vr_self / lambda
        : 0.0;
    const gmti::sim_geometry::LocalVelocity target_vel_local =
        gmti::sim_geometry::enuVelocityToLocal(
            gmti::sim_geometry::ENUVelocity(rt.ve_mps, rt.vn_mps, 0.0),
            cfg.geometry);

    gmti::target_injection::TargetConfig target;
    target.id = static_cast<int>(run.targets.size()) + 1;
    target.name = rt.target_id;
    target.enabled = rt.enabled;
    target.start_period = cfg.sim.period_start;
    target.end_period = cfg.sim.period_start + std::max(0, cfg.sim.period_count - 1);
    target.height_m = 0.0;
    target.motion_model = "constant_velocity";
    target.init_mode = "project_local";
    target.v = gmti::target_injection::Vec3(target_vel_local.vx, target_vel_local.vy, target_vel_local.vz);
    target.p0 = target_ref - target.v * ref_time_s;
    target.slant_range_m = range_m;
    target.azimuth_deg = rt.azimuth_deg;
    target.has_ref_geometry = true;
    target.ref_beam_id = beam0;
    target.ref_pulse_idx = ref_pulse_idx;
    target.ref_time_s = ref_time_s;
    target.ref_platform = ref_platform.position;
    target.ref_target = target_ref;
    target.ref_range_m = range_m;
    target.ref_range_sample_float = range_sample_float;
    target.ref_range_sample_int = range_sample_int;
    target.echo_delay_sample_center_used =
        (2.0 * gmti::target_injection::norm(target_ref - ref_platform.position) /
             gmti::target_injection::kC -
         cfg.radar.sample_delay_sec) *
        cfg.radar.fs_hz;
    target.target_snr_db = rt.snr_db;
    target.velocity_mode = rt.motion_type;
    target.target_ve_mps = rt.ve_mps;
    target.target_vn_mps = rt.vn_mps;
    target.target_vr_self_mps = target_vr_self;
    target.target_vt_self_mps = target_vt_self;
    target.af_motion_truth_hz = af_motion_truth;
    target.override_speed_mps = std::sqrt(rt.ve_mps * rt.ve_mps + rt.vn_mps * rt.vn_mps);
    target.visibility_mode = rt.visibility_type;
    rt.target = target;
    (void)range_sample_int;
    return true;
}

} // namespace

bool loadStage2RunConfig(const std::string &path, Stage2RunConfig &run, std::string &err)
{
    const std::string txt = readTextFile(path);
    if (txt.empty()) {
        err = "failed to read run config: " + path;
        return false;
    }
    Stage2Config cfg;

    run.case_id = jsonString(txt, "case_id", "");
    run.output_dir = jsonString(txt, "output_dir", "");
    run.scene_mode = jsonString(txt, "scene_mode", run.scene_mode);
    run.truth_output = jsonBool(txt, "truth_output", true);
    Stage2RunOptions scene_opt;

    const std::string waveform = sectionObject(txt, "waveform");
    if (!waveform.empty()) {
        cfg.radar.fc_hz = jsonDouble(waveform, "fc_ghz", cfg.radar.fc_hz / 1.0e9) * 1.0e9;
        cfg.radar.br_hz = jsonDouble(waveform, "bandwidth_mhz", cfg.radar.br_hz / 1.0e6) * 1.0e6;
        cfg.radar.fs_hz = jsonDouble(waveform, "fs_mhz", cfg.radar.fs_hz / 1.0e6) * 1.0e6;
        cfg.radar.tr_sec = jsonDouble(waveform, "tr_us", cfg.radar.tr_sec * 1.0e6) * 1.0e-6;
        cfg.radar.prf_hz = jsonDouble(waveform, "prf_hz", cfg.radar.prf_hz);
        cfg.radar.pulse_len = jsonInt(waveform, "pulse_len", cfg.radar.pulse_len);
        cfg.radar.pulse_num = jsonInt(waveform, "pulse_num", cfg.radar.pulse_num);
        cfg.radar.d_chan_m = jsonDouble(waveform, "d_chan_m", cfg.radar.d_chan_m);
        cfg.radar.iq_data_type = jsonString(waveform, "iq_data_type", cfg.radar.iq_data_type);
        cfg.radar.new_protocol_channel_count = jsonInt(waveform, "new_protocol_channel_count", cfg.radar.new_protocol_channel_count);
        cfg.radar.new_protocol_read_channel_1 = jsonInt(waveform, "new_protocol_read_channel_1", cfg.radar.new_protocol_read_channel_1);
        cfg.radar.new_protocol_read_channel_2 = jsonInt(waveform, "new_protocol_read_channel_2", cfg.radar.new_protocol_read_channel_2);
    }
    const std::string range_processing = sectionObject(txt, "range_processing");
    if (!range_processing.empty()) {
        cfg.radar.range_fft_len = jsonInt(range_processing, "range_fft_len", cfg.radar.range_fft_len);
        cfg.radar.range_crop_start = jsonInt(range_processing, "range_crop_start", cfg.radar.range_crop_start);
        cfg.radar.range_crop_len = jsonInt(range_processing, "range_crop_len", cfg.radar.range_crop_len);
        cfg.radar.sample_delay_sec =
            jsonDouble(range_processing, "sample_delay_us", cfg.radar.sample_delay_sec * 1.0e6) * 1.0e-6;
    }
    const std::string scan = sectionObject(txt, "scan");
    if (!scan.empty()) {
        cfg.radar.scan_min_deg = jsonDouble(scan, "scan_min_deg", cfg.radar.scan_min_deg);
        cfg.radar.scan_step_deg = jsonDouble(scan, "scan_step_deg", cfg.radar.scan_step_deg);
        cfg.radar.beam_count = jsonInt(scan, "beam_count", cfg.radar.beam_count);
        cfg.radar.beam_width_deg = jsonDouble(scan, "beam_width_deg", cfg.radar.beam_width_deg);
        run.beam_index_base = jsonInt(scan, "beam_index_base", 1);
    }
    const std::string platform = sectionObject(txt, "platform");
    if (!platform.empty()) {
        cfg.platform_speed_mps = jsonDouble(platform, "speed_mps", cfg.platform_speed_mps);
        cfg.platform_height_m = jsonDouble(platform, "height_m", cfg.platform_height_m);
        cfg.platform_origin_lat_deg = jsonDouble(platform, "origin_lat_deg", cfg.platform_origin_lat_deg);
        cfg.platform_origin_lon_deg = jsonDouble(platform, "origin_lon_deg", cfg.platform_origin_lon_deg);
        cfg.platform_origin_alt_m = jsonDouble(platform, "origin_alt_m", cfg.platform_origin_alt_m);
        cfg.projection_ref_lon_deg = jsonDouble(platform, "projection_ref_lon_deg", cfg.projection_ref_lon_deg);
        cfg.geometry.squint_side = jsonInt(platform, "squint_side", cfg.geometry.squint_side);
    }
    const std::string random = sectionObject(txt, "random");
    if (!random.empty()) {
        cfg.sim.period_start = jsonInt(random, "period_start", cfg.sim.period_start);
        cfg.sim.period_count = jsonInt(random, "period_count", cfg.sim.period_count);
        cfg.sim.random_seed = static_cast<uint32_t>(jsonInt(random, "random_seed", static_cast<int>(cfg.sim.random_seed)));
    }
    const std::string scene = sectionObject(txt, "scene");
    if (!scene.empty()) {
        run.scene_mode = jsonString(scene, "mode", run.scene_mode);
        cfg.scene.range_min_m = jsonDouble(scene, "range_min_m", cfg.scene.range_min_m);
        cfg.scene.range_max_m = jsonDouble(scene, "range_max_m", cfg.scene.range_max_m);
        cfg.scene.azimuth_min_deg = jsonDouble(scene, "azimuth_min_deg", cfg.scene.azimuth_min_deg);
        cfg.scene.azimuth_max_deg = jsonDouble(scene, "azimuth_max_deg", cfg.scene.azimuth_max_deg);
        cfg.scene.ground_z_m = jsonDouble(scene, "ground_z_m", cfg.scene.ground_z_m);
        cfg.scene.clutter_amplitude_scale =
            jsonDouble(scene, "clutter_amplitude_scale", cfg.scene.clutter_amplitude_scale);

        const std::string single = sectionObject(scene, "single_point");
        if (!single.empty()) {
            scene_opt.single_scatterer_range_m =
                jsonDouble(single, "range_m", scene_opt.single_scatterer_range_m);
            scene_opt.single_point_beam_id_1based =
                jsonInt(single, "beam_id", scene_opt.single_point_beam_id_1based);
            scene_opt.single_point_expected_bin =
                jsonInt(single, "expected_bin", scene_opt.single_point_expected_bin);
            scene_opt.single_scatterer_azimuth_deg =
                jsonDouble(single, "azimuth_deg", scene_opt.single_scatterer_azimuth_deg);
            scene_opt.single_scatterer_amplitude =
                jsonDouble(single, "amplitude", scene_opt.single_scatterer_amplitude);
        }

        const std::string area = sectionObject(scene, "area_clutter");
        if (!area.empty()) {
            cfg.scene.area.enabled = jsonBool(area, "enabled", cfg.scene.area.enabled);
            cfg.scene.area.model = jsonString(area, "model", cfg.scene.area.model);
            cfg.scene.area.scatterer_count =
                jsonInt(area, "scatterer_count", cfg.scene.area.scatterer_count);
            cfg.scene.area.mean_power = jsonDouble(area, "mean_power", cfg.scene.area.mean_power);
            cfg.scene.area.texture_sigma =
                jsonDouble(area, "texture_sigma", cfg.scene.area.texture_sigma);
            cfg.scene.area.spatial_cell_m =
                jsonDouble(area, "spatial_cell_m", cfg.scene.area.spatial_cell_m);
        }

        const std::string strong = sectionObject(scene, "strong_scatterers");
        if (!strong.empty()) {
            cfg.scene.strong.enabled = jsonBool(strong, "enabled", cfg.scene.strong.enabled);
            cfg.scene.strong.count = jsonInt(strong, "count", cfg.scene.strong.count);
            cfg.scene.strong.rcs_db_min =
                jsonDouble(strong, "rcs_db_min", cfg.scene.strong.rcs_db_min);
            cfg.scene.strong.rcs_db_max =
                jsonDouble(strong, "rcs_db_max", cfg.scene.strong.rcs_db_max);
        }

        const std::string line = sectionObject(scene, "line_scatterers");
        if (!line.empty()) {
            cfg.scene.line.enabled = jsonBool(line, "enabled", cfg.scene.line.enabled);
            cfg.scene.line.line_count = jsonInt(line, "line_count", cfg.scene.line.line_count);
            cfg.scene.line.points_per_line =
                jsonInt(line, "points_per_line", cfg.scene.line.points_per_line);
            cfg.scene.line.rcs_db = jsonDouble(line, "rcs_db", cfg.scene.line.rcs_db);
        }

        const std::string noise = sectionObject(scene, "thermal_noise");
        if (!noise.empty()) {
            cfg.scene.noise.enabled = jsonBool(noise, "enabled", cfg.scene.noise.enabled);
            cfg.scene.noise.noise_power =
                jsonDouble(noise, "noise_power", cfg.scene.noise.noise_power);
        }
    }
    cfg.geometry.platform_origin_lat_deg = cfg.platform_origin_lat_deg;
    cfg.geometry.platform_origin_lon_deg = cfg.platform_origin_lon_deg;
    cfg.geometry.platform_origin_alt_m = cfg.platform_origin_alt_m;
    cfg.geometry.projection_ref_lon_deg = cfg.projection_ref_lon_deg;
    cfg.target.enabled = true;
    run.cfg = cfg;
    run.global = makeTargetGlobal(cfg);
    scene_opt.output_dir = run.output_dir;
    scene_opt.scene_mode = run.scene_mode;
    scene_opt.target_enabled = true;
    run.legacy_scene_options = scene_opt;
    run.use_legacy_scene_options = true;

    const std::vector<std::string> target_objs = arrayObjects(txt, "targets");
    for (size_t i = 0; i < target_objs.size(); ++i) {
        const std::string &obj = target_objs[i];
        Stage2RunTarget rt;
        rt.target_id = jsonString(obj, "target_id", "target_" + std::to_string(i + 1));
        rt.enabled = jsonBool(obj, "enabled", true);
        const std::string init = sectionObject(obj, "init");
        const std::string motion = sectionObject(obj, "motion");
        const std::string amplitude = sectionObject(obj, "amplitude");
        const std::string visibility = sectionObject(obj, "visibility");
        rt.init_type = jsonType(init, "");
        rt.beam_id = jsonInt(init, "beam_id", -1);
        rt.expected_bin = jsonInt(init, "expected_bin", -1);
        rt.azimuth_deg = jsonDouble(init, "azimuth_deg", std::numeric_limits<double>::quiet_NaN());
        rt.azimuth_offset_deg = jsonDouble(init, "azimuth_offset_deg", std::numeric_limits<double>::quiet_NaN());
        rt.motion_type = jsonType(motion, "static");
        if (rt.motion_type == "enu_velocity") {
            rt.ve_mps = jsonDouble(motion, "ve_mps", 0.0);
            rt.vn_mps = jsonDouble(motion, "vn_mps", 0.0);
        }
        rt.amplitude_type = jsonType(amplitude, "snr_db");
        rt.snr_db = jsonDouble(amplitude, "snr_db", cfg.target.target_snr_db);
        rt.visibility_type = jsonType(visibility, run.global.visibility_mode);
        if (!resolveRunTarget(run, rt, err)) {
            return false;
        }
        run.targets.push_back(rt);
    }
    return validateStage2RunConfig(run, err);
}

bool validateStage2RunConfig(const Stage2RunConfig &run, std::string &err)
{
    if (run.case_id.empty()) { err = "case_id is required"; return false; }
    if (run.output_dir.empty()) { err = "output_dir is required"; return false; }
    if (run.targets.empty()) { err = "targets must not be empty"; return false; }
    if (run.cfg.radar.scan_step_deg == 0.0) { err = "scan_step_deg must be non-zero"; return false; }
    if (run.beam_index_base != 0 && run.beam_index_base != 1) { err = "beam_index_base must be 0 or 1"; return false; }
    std::set<std::string> ids;
    for (size_t i = 0; i < run.targets.size(); ++i) {
        const Stage2RunTarget &t = run.targets[i];
        if (!ids.insert(t.target_id).second) { err = "duplicate target_id: " + t.target_id; return false; }
        const int beam0 = t.beam_id - run.beam_index_base;
        if (beam0 < 0 || beam0 >= run.cfg.radar.beam_count) { err = "invalid beam_id for " + t.target_id; return false; }
        if (t.expected_bin < 0 || t.expected_bin >= run.cfg.radar.range_crop_len) { err = "invalid expected_bin for " + t.target_id; return false; }
        if (t.init_type != "beam_bin_azimuth" && t.init_type != "beam_bin_azimuth_offset") { err = "invalid init.type for " + t.target_id; return false; }
        if (t.motion_type != "static" && t.motion_type != "enu_velocity") { err = "invalid motion.type for " + t.target_id; return false; }
        if (t.amplitude_type != "snr_db") { err = "invalid amplitude.type for " + t.target_id; return false; }
    }
    return true;
}

bool makeLegacyStage2RunConfig(const Stage2RunOptions &opt, Stage2RunConfig &run, std::string &err)
{
    Stage2Config cfg;
    if (!loadStage2Config(opt.stage2_config, cfg, err)) {
        return false;
    }
    if (opt.period_start >= 0) cfg.sim.period_start = opt.period_start;
    if (opt.period_count > 0) cfg.sim.period_count = opt.period_count;
    run.case_id = "legacy_stage2";
    run.output_dir = opt.output_dir;
    run.scene_mode = opt.scene_mode;
    run.beam_index_base = 1;
    run.cfg = cfg;
    run.global = makeTargetGlobal(cfg);
    run.legacy_scene_options = opt;
    run.use_legacy_scene_options = true;
    if (opt.target_enabled && (opt.moving_target_beam_id_1based >= 1 || opt.moving_target_expected_bin >= 0)) {
        Stage2RunTarget rt;
        rt.target_id = "legacy_moving_target";
        rt.init_type = "beam_bin_azimuth";
        rt.beam_id = opt.moving_target_beam_id_1based >= 1 ? opt.moving_target_beam_id_1based : 1;
        rt.expected_bin = opt.moving_target_expected_bin >= 0 ? opt.moving_target_expected_bin : 0;
        const int beam0 = rt.beam_id - 1;
        rt.azimuth_deg = cfg.radar.scan_min_deg + static_cast<double>(beam0) * cfg.radar.scan_step_deg;
        rt.motion_type = "enu_velocity";
        rt.ve_mps = std::isfinite(opt.moving_target_ve_mps) ? opt.moving_target_ve_mps : -opt.moving_target_speed_mps;
        rt.vn_mps = std::isfinite(opt.moving_target_vn_mps) ? opt.moving_target_vn_mps : 0.0;
        rt.amplitude_type = "snr_db";
        rt.snr_db = cfg.target.target_snr_db;
        rt.visibility_type = run.global.visibility_mode;
        if (!resolveRunTarget(run, rt, err)) return false;
        run.targets.push_back(rt);
    }
    return true;
}

} // namespace stage2
} // namespace gmti
