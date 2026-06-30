#include "stage1_beam_mapper.h"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>

namespace gmti {
namespace sim_stage1 {

namespace {

double wrapToOldAngleRange(double theta)
{
    double ref = theta;
    while (ref > 25.0) ref -= 52.0;
    while (ref < -25.0) ref += 52.0;
    return ref;
}

} // namespace

std::vector<double> buildOldAngles()
{
    std::vector<double> a(26);
    for (int i = 0; i < 26; ++i) a[static_cast<size_t>(i)] = -25.0 + 2.0 * i;
    return a;
}

std::vector<double> buildNewAngles(const Stage1NewSystemConfig &cfg)
{
    std::vector<double> a(static_cast<size_t>(cfg.beam_count));
    for (int i = 0; i < cfg.beam_count; ++i) {
        a[static_cast<size_t>(i)] = cfg.scan_min_deg + cfg.scan_step_deg * i;
    }
    return a;
}

bool buildBeamMap(const Stage1NewSystemConfig &cfg, std::vector<BeamMapEntry> &map, std::string &err)
{
    const std::vector<double> new_angles = buildNewAngles(cfg);
    map.clear();
    map.reserve(new_angles.size());
    for (size_t k = 0; k < new_angles.size(); ++k) {
        const double theta = new_angles[k];
        const double ref = wrapToOldAngleRange(theta);
        const bool wrapped = std::fabs(ref - theta) > 1e-6;
        if (ref < -25.000001 || ref > 25.000001) {
            std::ostringstream oss;
            oss << "background_ref_angle out of old range for theta=" << theta << " ref=" << ref;
            err = oss.str();
            return false;
        }

        const double pos = (ref + 25.0) / 2.0;
        const int nearest = static_cast<int>(std::floor(pos + 0.5));
        BeamMapEntry e;
        e.new_beam_index = static_cast<int>(k);
        e.theta_new_deg = theta;
        e.background_ref_angle_deg = ref;
        e.is_wrapped = wrapped;
        if (std::fabs(pos - nearest) < 1e-6) {
            if (nearest < 0 || nearest > 25) {
                err = "direct beam index out of range";
                return false;
            }
            e.source_mode = "direct";
            e.source_left_beam_index = nearest;
            e.source_right_beam_index = nearest;
            e.source_left_angle_deg = -25.0 + 2.0 * nearest;
            e.source_right_angle_deg = e.source_left_angle_deg;
            e.interp_weight = 0.0;
        } else {
            int left = static_cast<int>(std::floor(pos));
            int right = left + 1;
            if (left < 0) {
                left = 0;
                right = 0;
            }
            if (right > 25) {
                left = 25;
                right = 25;
            }
            e.source_mode = "interp";
            e.source_left_beam_index = left;
            e.source_right_beam_index = right;
            e.source_left_angle_deg = -25.0 + 2.0 * left;
            e.source_right_angle_deg = -25.0 + 2.0 * right;
            e.interp_weight = (right == left) ? 0.0 : (ref - e.source_left_angle_deg) / 2.0;
            if (e.interp_weight < -1e-6 || e.interp_weight > 1.000001) {
                err = "interp weight out of range";
                return false;
            }
        }
        map.push_back(e);
    }
    return true;
}

bool writeBeamTables(const std::string &output_dir, const std::vector<BeamMapEntry> &map, std::string &err)
{
    const std::string debug_dir = pathJoin(output_dir, "debug");
    std::ofstream angle(pathJoin(debug_dir, "beam_angle_table.csv").c_str());
    std::ofstream beam(pathJoin(debug_dir, "old_to_new_beam_map.csv").c_str());
    if (!angle || !beam) {
        err = "failed to open beam map csv outputs";
        return false;
    }
    angle << "new_beam_index,theta_new_deg\n";
    beam << "new_beam_index,theta_new_deg,background_ref_angle_deg,source_mode,"
         << "source_left_beam_index,source_left_angle_deg,source_right_beam_index,"
         << "source_right_angle_deg,interp_weight,is_wrapped\n";
    angle << std::fixed << std::setprecision(6);
    beam << std::fixed << std::setprecision(6);
    for (size_t i = 0; i < map.size(); ++i) {
        const BeamMapEntry &e = map[i];
        angle << e.new_beam_index << "," << e.theta_new_deg << "\n";
        beam << e.new_beam_index << "," << e.theta_new_deg << ","
             << e.background_ref_angle_deg << "," << e.source_mode << ","
             << e.source_left_beam_index << "," << e.source_left_angle_deg << ","
             << e.source_right_beam_index << "," << e.source_right_angle_deg << ","
             << e.interp_weight << "," << boolText(e.is_wrapped) << "\n";
    }
    return true;
}

} // namespace sim_stage1
} // namespace gmti
