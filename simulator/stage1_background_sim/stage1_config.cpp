#include "stage1_config.h"

#include "tinyxml.h"

#include <cerrno>
#include <cstdlib>
#include <cstring>
#include <sstream>
#include <sys/stat.h>
#include <sys/types.h>
#include <vector>

namespace gmti {
namespace sim_stage1 {

namespace {

std::string textOf(TiXmlElement *parent, const char *name)
{
    if (!parent) return "";
    TiXmlElement *e = parent->FirstChildElement(name);
    if (!e || !e->GetText()) return "";
    return e->GetText();
}

int toInt(const std::string &s, int dflt)
{
    if (s.empty()) return dflt;
    return std::atoi(s.c_str());
}

double toDouble(const std::string &s, double dflt)
{
    if (s.empty()) return dflt;
    return std::atof(s.c_str());
}

bool mkdirOne(const std::string &p, std::string &err)
{
    if (p.empty()) return true;
    if (::mkdir(p.c_str(), 0775) == 0 || errno == EEXIST) return true;
    err = "mkdir failed for " + p + ": " + std::strerror(errno);
    return false;
}

bool mkdirRecursive(const std::string &p, std::string &err)
{
    std::string cur;
    for (size_t i = 0; i < p.size(); ++i) {
        cur.push_back(p[i]);
        if (p[i] == '/' || i + 1 == p.size()) {
            if (cur == "/" || cur.empty()) continue;
            if (cur[cur.size() - 1] == '/' && cur.size() > 1) {
                cur.resize(cur.size() - 1);
            }
            if (!cur.empty() && !mkdirOne(cur, err)) return false;
            if (i + 1 < p.size() && (cur.empty() || cur[cur.size() - 1] != '/')) {
                cur.push_back('/');
            }
        }
    }
    return true;
}

} // namespace

double Stage1NewSystemConfig::periodTimeSec() const
{
    return static_cast<double>(beam_count) * static_cast<double>(pulse_num_new) / prf_new_hz;
}

double Stage1NewSystemConfig::realtime80Sec() const
{
    return periodTimeSec() * 0.8;
}

bool parseOldConfigXml(const std::string &xml_path, Stage1OldSystemConfig &cfg, std::string &err)
{
    TiXmlDocument doc;
    if (!doc.LoadFile(xml_path.c_str())) {
        err = "load XML failed: " + xml_path;
        return false;
    }
    TiXmlElement *root = doc.FirstChildElement("GMTI");
    TiXmlElement *param = root ? root->FirstChildElement("GMTI_parameter") : nullptr;
    if (!param) {
        err = "XML missing <GMTI><GMTI_parameter>";
        return false;
    }

    cfg.xml_path = xml_path;
    cfg.data_ch1 = textOf(param, "GMTI_data");
    cfg.data_ch2 = textOf(param, "GMTI_data2");
    cfg.pos_path = textOf(param, "Plane_POS");
    cfg.result_add = textOf(param, "result_add");
    cfg.info_len = toInt(textOf(param, "info_len"), cfg.info_len);
    cfg.pulse_len = toInt(textOf(param, "pulse_len"), cfg.pulse_len);
    cfg.rg_len = toInt(textOf(param, "rg_len"), cfg.rg_len);
    cfg.pulse_num = toInt(textOf(param, "pulse_num"), cfg.pulse_num);
    cfg.read_pulse_num = toInt(textOf(param, "read_pulse_num"), cfg.read_pulse_num);
    cfg.read_pulse_offset = toInt(textOf(param, "read_pulse_offset"), cfg.read_pulse_offset);
    cfg.pulse_dec = toInt(textOf(param, "pulse_dec"), cfg.pulse_dec);
    cfg.fc_ghz = toDouble(textOf(param, "fc"), cfg.fc_ghz);
    cfg.br_mhz = toDouble(textOf(param, "Br"), cfg.br_mhz);
    cfg.fs_mhz = toDouble(textOf(param, "fs"), cfg.fs_mhz);
    cfg.tr_us = toDouble(textOf(param, "Tr"), cfg.tr_us);
    cfg.prf_hz = toDouble(textOf(param, "PRF"), cfg.prf_hz);
    cfg.scan_min_deg = toDouble(textOf(param, "scan_min_deg"), cfg.scan_min_deg);
    cfg.scan_max_deg = toDouble(textOf(param, "scan_max_deg"), cfg.scan_max_deg);
    cfg.az_count = toInt(textOf(param, "az_count"), cfg.az_count);
    cfg.wavepos_st = toInt(textOf(param, "wavepos_st"), cfg.wavepos_st);
    cfg.wavepos_ed = toInt(textOf(param, "wavepos_ed"), cfg.wavepos_ed);
    cfg.wavepos_skip = toInt(textOf(param, "wavepos_skip"), cfg.wavepos_skip);
    cfg.skip_pulses = toInt(textOf(param, "skip_pulses"), cfg.skip_pulses);
    cfg.week_offset = toInt(textOf(param, "week_offset"), cfg.week_offset);
    cfg.sec_bias = toInt(textOf(param, "secBias"), cfg.sec_bias);
    cfg.d_chan = toDouble(textOf(param, "d_chan"), cfg.d_chan);
    cfg.rmin_m = toDouble(textOf(param, "Rmin"), cfg.rmin_m);

    if (cfg.data_ch1.empty() || cfg.data_ch2.empty()) {
        err = "XML missing GMTI_data or GMTI_data2";
        return false;
    }
    return true;
}

bool parseCommandLine(int argc, char **argv, Stage1RunOptions &opt, Stage1NewSystemConfig &new_cfg)
{
    for (int i = 1; i < argc; ++i) {
        const std::string k = argv[i];
        const char *v = (i + 1 < argc) ? argv[i + 1] : nullptr;
        if (k == "--old-config" && v) opt.old_config_path = argv[++i];
        else if (k == "--output-dir" && v) opt.output_dir = argv[++i];
        else if (k == "--scan-min" && v) new_cfg.scan_min_deg = std::atof(argv[++i]);
        else if (k == "--scan-max" && v) new_cfg.scan_max_deg = std::atof(argv[++i]);
        else if (k == "--scan-step" && v) new_cfg.scan_step_deg = std::atof(argv[++i]);
        else if (k == "--beam-count" && v) new_cfg.beam_count = std::atoi(argv[++i]);
        else if (k == "--pulse-resample-mode" && v) opt.pulse_resample_mode = argv[++i];
        else if (k == "--range-resize-mode" && v) opt.range_resize_mode = argv[++i];
        else if (k == "--period-start" && v) opt.period_start = std::atoi(argv[++i]);
        else if (k == "--period-count" && v) opt.period_count = std::atoi(argv[++i]);
        else if (k == "--beam-start" && v) opt.beam_start = std::atoi(argv[++i]);
        else if (k == "--generate-beam-count" && v) opt.beam_count = std::atoi(argv[++i]);
        else if (k == "--channel" && v) opt.channel_mode = argv[++i];
        else if (k == "--write-config" && v) opt.write_config = (std::string(argv[++i]) != "false");
        else if (k == "--validate" && v) opt.validate = (std::string(argv[++i]) != "false");
        else if (k == "--generate-data" && v) opt.generate_data = (std::string(argv[++i]) != "false");
        else {
            std::ostringstream oss;
            oss << "unknown or incomplete option: " << k;
            throw std::runtime_error(oss.str());
        }
    }
    return true;
}

bool ensureStage1Dirs(const std::string &output_dir, std::string &err)
{
    std::vector<std::string> dirs;
    dirs.push_back(output_dir);
    dirs.push_back(pathJoin(output_dir, "config"));
    dirs.push_back(pathJoin(output_dir, "data"));
    dirs.push_back(pathJoin(output_dir, "debug"));
    dirs.push_back(pathJoin(output_dir, "logs"));
    dirs.push_back(pathJoin(output_dir, "reports"));
    for (size_t i = 0; i < dirs.size(); ++i) {
        if (!mkdirRecursive(dirs[i], err)) return false;
    }
    return true;
}

std::string pathJoin(const std::string &a, const std::string &b)
{
    if (a.empty()) return b;
    if (a[a.size() - 1] == '/') return a + b;
    return a + "/" + b;
}

std::string boolText(bool v)
{
    return v ? "true" : "false";
}

} // namespace sim_stage1
} // namespace gmti
