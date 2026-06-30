#include "stage2_config.h"

#include "tinyxml.h"

#include <cstdlib>
#include <fstream>
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
    return doc.SaveFile(out_xml.c_str());
}

} // namespace

bool parseStage2CommandLine(int argc, char **argv, Stage2RunOptions &opt)
{
    for (int i = 1; i < argc; ++i) {
        const std::string k = argv[i];
        const std::string v = (i + 1 < argc) ? argv[i + 1] : "";
        if (k == "--stage2-config") opt.stage2_config = valueArg(i, argc, argv);
        else if (k == "--target-config") opt.target_config = valueArg(i, argc, argv);
        else if (k == "--output-dir") opt.output_dir = valueArg(i, argc, argv);
        else if (k == "--scene-mode") opt.scene_mode = valueArg(i, argc, argv);
        else if (k == "--target-enabled") opt.target_enabled = (v != "false" && v != "0");
        else if (k == "--period-start") opt.period_start = parseIntArg(valueArg(i, argc, argv), opt.period_start);
        else if (k == "--period-count") opt.period_count = parseIntArg(valueArg(i, argc, argv), opt.period_count);
        else if (k == "--single-scatterer-range") opt.single_scatterer_range_m = parseDoubleArg(valueArg(i, argc, argv), opt.single_scatterer_range_m);
        else if (k == "--single-scatterer-azimuth") opt.single_scatterer_azimuth_deg = parseDoubleArg(valueArg(i, argc, argv), opt.single_scatterer_azimuth_deg);
        else if (k == "--single-scatterer-amplitude") opt.single_scatterer_amplitude = parseDoubleArg(valueArg(i, argc, argv), opt.single_scatterer_amplitude);
        else if (k == "--validate") opt.validate = (v != "false" && v != "0");
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
    cfg.radar.beam_width_deg = jsonDouble(sys, "beam_width_deg", 2.28);
    cfg.radar.pulse_num = jsonInt(sys, "pulse_num", 130);
    cfg.radar.d_chan_m = jsonDouble(sys, "d_chan_m", 0.17);
    cfg.platform_height_m = jsonDouble(sys, "platform_height_m", cfg.platform_height_m);
    cfg.platform_speed_mps = jsonDouble(sys, "platform_speed_mps", cfg.platform_speed_mps);

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
        "    \"beam_width_deg\": 2.28, \"pulse_num\": 130,\n"
        "    \"platform_height_m\": 6000.0, \"platform_speed_mps\": 60.0, \"d_chan_m\": 0.17\n"
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
    return g;
}

} // namespace stage2
} // namespace gmti
