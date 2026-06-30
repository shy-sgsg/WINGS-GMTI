#include "stage1_pulse_resampler.h"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>

namespace gmti {
namespace sim_stage1 {

bool buildPulseMap(const Stage1OldSystemConfig &old_cfg,
                   const Stage1NewSystemConfig &new_cfg,
                   const std::string &mode,
                   std::vector<PulseMapEntry> &map,
                   std::string &err)
{
    map.clear();
    map.reserve(static_cast<size_t>(new_cfg.pulse_num_new));
    for (int m = 0; m < new_cfg.pulse_num_new; ++m) {
        double pos = 0.0;
        if (mode == "physical_time") {
            pos = static_cast<double>(m) * old_cfg.prf_hz / new_cfg.prf_new_hz;
        } else if (mode == "normalized") {
            pos = static_cast<double>(m) * static_cast<double>(old_cfg.pulse_num - 1) /
                  static_cast<double>(new_cfg.pulse_num_new - 1);
        } else {
            err = "unknown pulse_resample_mode: " + mode;
            return false;
        }
        if (pos < -1e-6 || pos > static_cast<double>(old_cfg.pulse_num - 1) + 1e-6) {
            std::ostringstream oss;
            oss << "old pulse position out of range: m=" << m << " pos=" << pos;
            err = oss.str();
            return false;
        }
        const int left = static_cast<int>(std::floor(pos));
        int right = left + 1;
        if (right >= old_cfg.pulse_num) right = old_cfg.pulse_num - 1;
        PulseMapEntry e;
        e.new_pulse_index = m;
        e.old_pulse_pos = pos;
        e.old_left_index = left;
        e.old_right_index = right;
        e.weight = pos - static_cast<double>(left);
        if (std::fabs(e.weight) < 1e-6) {
            e.old_right_index = e.old_left_index;
            e.weight = 0.0;
        }
        map.push_back(e);
    }
    return true;
}

bool writePulseMap(const std::string &output_dir,
                   const std::string &mode,
                   const std::vector<PulseMapEntry> &map,
                   std::string &err)
{
    const std::string path = pathJoin(pathJoin(output_dir, "debug"),
                                      "pulse_resample_map_" + mode + ".csv");
    std::ofstream out(path.c_str());
    if (!out) {
        err = "failed to open " + path;
        return false;
    }
    out << "new_pulse_index,old_pulse_pos,old_left_index,old_right_index,weight\n";
    out << std::fixed << std::setprecision(8);
    for (size_t i = 0; i < map.size(); ++i) {
        const PulseMapEntry &e = map[i];
        out << e.new_pulse_index << "," << e.old_pulse_pos << ","
            << e.old_left_index << "," << e.old_right_index << "," << e.weight << "\n";
    }
    return true;
}

} // namespace sim_stage1
} // namespace gmti

