#include "target_config.h"

#include "tinyxml.h"

#include <cstdlib>
#include <fstream>
#include <sstream>

namespace gmti {
namespace target_injection {

namespace {

std::string readTextFile(const std::string &path)
{
    std::ifstream in(path.c_str());
    std::ostringstream ss;
    ss << in.rdbuf();
    return ss.str();
}

std::string textOf(TiXmlElement *root, const char *name)
{
    if (!root) return "";
    TiXmlElement *e = root->FirstChildElement(name);
    if (!e || !e->GetText()) return "";
    return e->GetText();
}

double toDouble(const std::string &s, double fallback)
{
    if (s.empty()) return fallback;
    char *end = nullptr;
    const double v = std::strtod(s.c_str(), &end);
    return end == s.c_str() ? fallback : v;
}

int toInt(const std::string &s, int fallback)
{
    if (s.empty()) return fallback;
    char *end = nullptr;
    const long v = std::strtol(s.c_str(), &end, 10);
    return end == s.c_str() ? fallback : static_cast<int>(v);
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

std::string firstArrayObject(const std::string &txt, const std::string &key)
{
    const std::string needle = "\"" + key + "\"";
    const size_t k = txt.find(needle);
    if (k == std::string::npos) return "";
    const size_t open_arr = txt.find('[', k);
    if (open_arr == std::string::npos) return "";
    const size_t open = txt.find('{', open_arr);
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

std::string valueArg(int &i, int argc, char **argv)
{
    if (i + 1 >= argc) return "";
    ++i;
    return argv[i];
}

void initTargetFromRangeAzimuth(const RadarConfig &radar,
                                const TargetGlobalConfig &global,
                                TargetConfig &target)
{
    const gmti::sim_geometry::LocalPoint ref_platform(0.0, 0.0, global.platform_height_m);
    const gmti::sim_geometry::LocalVelocity ref_velocity(global.platform_speed_mps, 0.0, 0.0);
    const gmti::sim_geometry::LocalPoint target_local =
        gmti::sim_geometry::makePointFromRangeAzimuth(ref_platform,
                                                      ref_velocity,
                                                      target.slant_range_m,
                                                      target.azimuth_deg,
                                                      global.platform_height_m,
                                                      target.height_m,
                                                      global.geometry);
    target.p0.x = target_local.x;
    target.p0.y = target_local.y;
    target.p0.z = target.height_m;
    const gmti::sim_geometry::ENUVelocity ref_enu_velocity =
        gmti::sim_geometry::localVelocityToEnu(ref_velocity, global.geometry);
    const gmti::sim_geometry::LookVectorEN look =
        gmti::sim_geometry::makeAlgorithmLookVectorEN(ref_enu_velocity.ve,
                                                      ref_enu_velocity.vn,
                                                      target.azimuth_deg,
                                                      global.geometry);
    const gmti::sim_geometry::ENUVelocity target_enu_velocity(
        target.radial_velocity_mps * look.east,
        target.radial_velocity_mps * look.north,
        0.0);
    const gmti::sim_geometry::LocalVelocity target_local_velocity =
        gmti::sim_geometry::enuVelocityToLocal(target_enu_velocity, global.geometry);
    target.v.x = target_local_velocity.vx;
    target.v.y = target_local_velocity.vy;
    target.v.z = target_local_velocity.vz;
    (void)radar;
}

} // namespace

bool parseCommandLine(int argc, char **argv, RunOptions &opt)
{
    for (int i = 1; i < argc; ++i) {
        const std::string k = argv[i];
        const std::string v = (i + 1 < argc) ? argv[i + 1] : "";
        if (k == "--input-config") opt.input_config = valueArg(i, argc, argv);
        else if (k == "--input-data-dir") opt.input_data_dir = valueArg(i, argc, argv);
        else if (k == "--target-config") opt.target_config = valueArg(i, argc, argv);
        else if (k == "--output-dir") opt.output_dir = valueArg(i, argc, argv);
        else if (k == "--background-mode") opt.background_mode = valueArg(i, argc, argv);
        else if (k == "--visibility-mode") opt.visibility_mode = valueArg(i, argc, argv);
        else if (k == "--amplitude-mode") opt.amplitude_mode = valueArg(i, argc, argv);
        else if (k == "--period-start") opt.period_start = toInt(valueArg(i, argc, argv), opt.period_start);
        else if (k == "--period-count") opt.period_count = toInt(valueArg(i, argc, argv), opt.period_count);
        else if (k == "--beam-start") opt.beam_start = toInt(valueArg(i, argc, argv), opt.beam_start);
        else if (k == "--beam-count") opt.beam_count = toInt(valueArg(i, argc, argv), opt.beam_count);
        else if (k == "--pulse-start") opt.pulse_start = toInt(valueArg(i, argc, argv), opt.pulse_start);
        else if (k == "--pulse-count") opt.pulse_count = toInt(valueArg(i, argc, argv), opt.pulse_count);
        else if (k == "--target-snr-db") opt.target_snr_db = toDouble(valueArg(i, argc, argv), opt.target_snr_db);
        else if (k == "--direct-amplitude") opt.direct_amplitude = toDouble(valueArg(i, argc, argv), opt.direct_amplitude);
        else if (k == "--write-truth") opt.write_truth = (v != "false" && v != "0");
        else if (k == "--validate") opt.validate = (v != "false" && v != "0");
        else if (k == "--write-config") opt.write_config = (v != "false" && v != "0");
        else if (k == "--help") return false;
    }
    return true;
}

bool loadRadarConfig(const std::string &xml_path, RadarConfig &cfg, std::string &err)
{
    TiXmlDocument doc(xml_path.c_str());
    if (!doc.LoadFile()) {
        err = "failed to load xml: " + xml_path;
        return false;
    }
    TiXmlElement *root = doc.RootElement();
    TiXmlElement *param = root ? root->FirstChildElement("GMTI_parameter") : nullptr;
    if (!param && root && std::string(root->Value()) == "GMTI_parameter") param = root;
    if (!param) {
        err = "missing GMTI_parameter in xml: " + xml_path;
        return false;
    }
    cfg.input_config = xml_path;
    cfg.input_data_file = textOf(param, "GMTI_data_new");
    cfg.pulse_len = toInt(textOf(param, "pulse_len"), cfg.pulse_len);
    cfg.pulse_num = toInt(textOf(param, "pulse_num"), cfg.pulse_num);
    cfg.beam_count = toInt(textOf(param, "az_count"), cfg.beam_count);
    cfg.range_fft_len = toInt(textOf(param, "range_fft_len"), cfg.range_fft_len);
    cfg.range_crop_start = toInt(textOf(param, "range_crop_start"), cfg.range_crop_start);
    cfg.range_crop_len = toInt(textOf(param, "range_compress_len"), cfg.range_crop_len);
    cfg.scan_min_deg = toDouble(textOf(param, "scan_min_deg"), cfg.scan_min_deg);
    cfg.scan_step_deg = toDouble(textOf(param, "scan_step_deg"), cfg.scan_step_deg);
    cfg.beam_width_deg = toDouble(textOf(param, "beam_width_deg"), cfg.beam_width_deg);
    cfg.fc_hz = toDouble(textOf(param, "fc"), cfg.fc_hz / 1.0e9) * 1.0e9;
    cfg.br_hz = toDouble(textOf(param, "Br"), cfg.br_hz / 1.0e6) * 1.0e6;
    cfg.fs_hz = toDouble(textOf(param, "fs"), cfg.fs_hz / 1.0e6) * 1.0e6;
    cfg.tr_sec = toDouble(textOf(param, "Tr"), cfg.tr_sec * 1.0e6) * 1.0e-6;
    cfg.prf_hz = toDouble(textOf(param, "PRF"), cfg.prf_hz);
    cfg.sample_delay_sec = toDouble(textOf(param, "sample_delay_us"), cfg.sample_delay_sec * 1.0e6) * 1.0e-6;
    cfg.d_chan_m = toDouble(textOf(param, "d_chan"), cfg.d_chan_m);
    return true;
}

bool loadTargetConfig(const std::string &json_path,
                      const RadarConfig &radar,
                      TargetGlobalConfig &global,
                      TargetConfig &target,
                      std::string &err)
{
    std::ifstream in(json_path.c_str());
    if (!in) {
        if (!writeDefaultTargetsJson(json_path, err)) return false;
        in.open(json_path.c_str());
    }
    const std::string txt = readTextFile(json_path);
    const std::string g = sectionObject(txt, "global");
    const std::string t = firstArrayObject(txt, "targets");
    if (g.empty() || t.empty()) {
        err = "targets.json must contain global and one target";
        return false;
    }
    global.coordinate_mode = jsonString(g, "coordinate_mode", global.coordinate_mode);
    global.visibility_mode = jsonString(g, "visibility_mode", global.visibility_mode);
    global.amplitude_mode = jsonString(g, "amplitude_mode", global.amplitude_mode);
    global.channel_phase_mode = jsonString(g, "channel_phase_mode", global.channel_phase_mode);
    global.platform_mode = jsonString(g, "platform_mode", global.platform_mode);
    global.target_snr_db = jsonDouble(g, "target_snr_db", global.target_snr_db);
    global.direct_amplitude = jsonDouble(g, "direct_amplitude", global.direct_amplitude);
    global.beam_gain_threshold = jsonDouble(g, "beam_gain_threshold", global.beam_gain_threshold);
    global.rms_floor = jsonDouble(g, "rms_floor", global.rms_floor);
    global.chirp_phase_sign = jsonInt(g, "chirp_phase_sign", global.chirp_phase_sign);
    global.carrier_phase_sign = jsonInt(g, "carrier_phase_sign", global.carrier_phase_sign);
    global.platform_speed_mps = jsonDouble(g, "platform_speed_mps", global.platform_speed_mps);
    global.platform_height_m = jsonDouble(g, "platform_height_m", global.platform_height_m);
    global.platform_origin_lat_deg = jsonDouble(g, "platform_origin_lat_deg", global.platform_origin_lat_deg);
    global.platform_origin_lon_deg = jsonDouble(g, "platform_origin_lon_deg", global.platform_origin_lon_deg);
    global.platform_origin_alt_m = jsonDouble(g, "platform_origin_alt_m", global.platform_origin_alt_m);
    global.projection_ref_lon_deg = jsonDouble(g, "projection_ref_lon_deg", global.projection_ref_lon_deg);
    global.geometry.platform_origin_lat_deg = global.platform_origin_lat_deg;
    global.geometry.platform_origin_lon_deg = global.platform_origin_lon_deg;
    global.geometry.platform_origin_alt_m = global.platform_origin_alt_m;
    global.geometry.projection_ref_lon_deg = global.projection_ref_lon_deg;
    const std::string geom = sectionObject(txt, "simulation_geometry");
    if (!geom.empty()) {
        global.geometry.geometry_config_name = jsonString(geom, "geometry_config_name", global.geometry.geometry_config_name);
        global.geometry.local_x_axis = jsonString(geom, "local_x_axis", global.geometry.local_x_axis);
        global.geometry.local_y_axis = jsonString(geom, "local_y_axis", global.geometry.local_y_axis);
        global.geometry.platform_heading_source = jsonString(geom, "platform_heading_source", global.geometry.platform_heading_source);
        global.geometry.platform_heading_deg = jsonDouble(geom, "platform_heading_deg", global.geometry.platform_heading_deg);
        global.geometry.beam_angle_reference = jsonString(geom, "beam_angle_reference", global.geometry.beam_angle_reference);
        global.geometry.beam_zero_direction = jsonString(geom, "beam_zero_direction", global.geometry.beam_zero_direction);
        global.geometry.beam_positive_direction = jsonString(geom, "beam_positive_direction", global.geometry.beam_positive_direction);
        global.geometry.beam_theta_offset_deg = jsonDouble(geom, "beam_theta_offset_deg", global.geometry.beam_theta_offset_deg);
        global.geometry.range_geometry = jsonString(geom, "range_geometry", global.geometry.range_geometry);
        global.geometry.use_ground_range_for_position = jsonBool(geom, "use_ground_range_for_position", global.geometry.use_ground_range_for_position);
        global.geometry.platform_origin_lat_deg = jsonDouble(geom, "platform_origin_lat_deg", global.geometry.platform_origin_lat_deg);
        global.geometry.platform_origin_lon_deg = jsonDouble(geom, "platform_origin_lon_deg", global.geometry.platform_origin_lon_deg);
        global.geometry.platform_origin_alt_m = jsonDouble(geom, "platform_origin_alt_m", global.geometry.platform_origin_alt_m);
        global.geometry.projection_ref_lon_deg = jsonDouble(geom, "projection_ref_lon_deg", global.geometry.projection_ref_lon_deg);
        global.geometry.squint_side = jsonInt(geom, "squint_side", global.geometry.squint_side);
        global.platform_origin_lat_deg = global.geometry.platform_origin_lat_deg;
        global.platform_origin_lon_deg = global.geometry.platform_origin_lon_deg;
        global.platform_origin_alt_m = global.geometry.platform_origin_alt_m;
        global.projection_ref_lon_deg = global.geometry.projection_ref_lon_deg;
    }

    target.id = jsonInt(t, "id", target.id);
    target.name = jsonString(t, "name", target.name);
    target.motion_model = jsonString(t, "motion_model", target.motion_model);
    target.init_mode = jsonString(t, "init_mode", target.init_mode);
    target.enabled = jsonBool(t, "enabled", target.enabled);
    target.start_period = jsonInt(t, "start_period", target.start_period);
    target.end_period = jsonInt(t, "end_period", target.end_period);
    target.p0.x = jsonDouble(t, "x0_m", target.p0.x);
    target.p0.y = jsonDouble(t, "y0_m", target.p0.y);
    target.p0.z = jsonDouble(t, "z0_m", target.p0.z);
    target.v.x = jsonDouble(t, "vx_mps", target.v.x);
    target.v.y = jsonDouble(t, "vy_mps", target.v.y);
    target.v.z = jsonDouble(t, "vz_mps", target.v.z);
    target.slant_range_m = jsonDouble(t, "slant_range_m", target.slant_range_m);
    target.azimuth_deg = jsonDouble(t, "azimuth_deg", target.azimuth_deg);
    target.height_m = jsonDouble(t, "height_m", target.height_m);
    target.radial_velocity_mps = jsonDouble(t, "radial_velocity_mps", target.radial_velocity_mps);
    target.cross_velocity_mps = jsonDouble(t, "cross_velocity_mps", target.cross_velocity_mps);
    if (target.init_mode == "range_azimuth") {
        initTargetFromRangeAzimuth(radar, global, target);
    }
    return true;
}

void applyRunOverrides(const RunOptions &run, TargetGlobalConfig &global)
{
    if (!run.visibility_mode.empty()) global.visibility_mode = run.visibility_mode;
    if (!run.amplitude_mode.empty()) global.amplitude_mode = run.amplitude_mode;
    if (run.target_snr_db > -900.0) global.target_snr_db = run.target_snr_db;
    if (run.direct_amplitude >= 0.0) global.direct_amplitude = run.direct_amplitude;
}

bool writeOutputConfig(const std::string &input_xml,
                       const std::string &out_xml,
                       const std::string &out_data_file,
                       std::string &err)
{
    TiXmlDocument doc(input_xml.c_str());
    if (!doc.LoadFile()) {
        err = "failed to load input xml for output config";
        return false;
    }
    TiXmlElement *root = doc.RootElement();
    TiXmlElement *param = root ? root->FirstChildElement("GMTI_parameter") : nullptr;
    if (!param && root && std::string(root->Value()) == "GMTI_parameter") param = root;
    if (!param) {
        err = "missing GMTI_parameter for output config";
        return false;
    }
    TiXmlElement *node = param->FirstChildElement("GMTI_data_new");
    if (!node) {
        node = new TiXmlElement("GMTI_data_new");
        param->LinkEndChild(node);
    }
    node->Clear();
    node->LinkEndChild(new TiXmlText(out_data_file.c_str()));

    auto setTextNode = [&](const char *name, const char *value) {
        TiXmlElement *e = param->FirstChildElement(name);
        if (!e) {
            e = new TiXmlElement(name);
            param->LinkEndChild(e);
        }
        e->Clear();
        e->LinkEndChild(new TiXmlText(value));
    };
    setTextNode("raw_fenbianlv", "25");
    setTextNode("motion_comp_enable", "1");
    setTextNode("motion_comp_analytic_enable", "1");
    setTextNode("motion_comp_use_row_doppler", "1");
    setTextNode("motion_comp_iter", "3");
    setTextNode("ati_velocity_sign", "1");
    setTextNode("ati_phase_to_velocity_sign", "1");
    setTextNode("motion_doppler_axis_sign", "-1");
    setTextNode("ati_phase_bias_rad", "0.0");
    setTextNode("ati_vmax_mps", "60.0");
    setTextNode("motion_comp_denom_min", "1e-6");
    setTextNode("motion_comp_debug", "1");
    if (!doc.SaveFile(out_xml.c_str())) {
        err = "failed to save output config";
        return false;
    }
    return true;
}

bool writeDefaultTargetsJson(const std::string &path, std::string &err)
{
    std::ofstream out(path.c_str());
    if (!out) {
        err = "failed to create default target config: " + path;
        return false;
    }
    out <<
        "{\n"
        "  \"global\": {\n"
        "    \"coordinate_mode\": \"project_local\",\n"
        "    \"visibility_mode\": \"hard_gate\",\n"
        "    \"amplitude_mode\": \"snr_db\",\n"
        "    \"target_snr_db\": 30.0,\n"
        "    \"direct_amplitude\": 1.0,\n"
        "    \"chirp_phase_sign\": 1,\n"
        "    \"carrier_phase_sign\": -1,\n"
        "    \"channel_phase_mode\": \"baseline_approx\",\n"
        "    \"platform_mode\": \"ideal_platform\",\n"
        "    \"platform_speed_mps\": 60.0,\n"
        "    \"platform_height_m\": 6000.0,\n"
        "    \"platform_origin_lat_deg\": 40.45121057,\n"
        "    \"platform_origin_lon_deg\": 116.98377429,\n"
        "    \"platform_origin_alt_m\": 6000.0,\n"
        "    \"projection_ref_lon_deg\": 117.0,\n"
        "    \"beam_gain_threshold\": 0.05,\n"
        "    \"rms_floor\": 1e-6\n"
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
        "  \"targets\": [\n"
        "    {\n"
        "      \"id\": 1,\n"
        "      \"name\": \"strong_cooperative_target\",\n"
        "      \"motion_model\": \"constant_velocity\",\n"
        "      \"init_mode\": \"range_azimuth\",\n"
        "      \"slant_range_m\": 85000.0,\n"
        "      \"azimuth_deg\": 0.0,\n"
        "      \"height_m\": 0.0,\n"
        "      \"radial_velocity_mps\": 10.0,\n"
        "      \"cross_velocity_mps\": 0.0,\n"
        "      \"start_period\": 0,\n"
        "      \"end_period\": 10,\n"
        "      \"enabled\": true\n"
        "    }\n"
        "  ]\n"
        "}\n";
    return true;
}

} // namespace target_injection
} // namespace gmti
