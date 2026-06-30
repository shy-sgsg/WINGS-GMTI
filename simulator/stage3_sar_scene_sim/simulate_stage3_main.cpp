#include "stage3_config.h"
#include "stage3_validator.h"

#include <fstream>
#include <iostream>

using namespace gmti::stage3;

int main(int argc, char **argv)
{
    Stage3RunOptions opt;
    std::string err;
    if (!parseStage3CommandLine(argc, argv, opt, err)) {
        std::cerr << "[stage3][ERR] " << err << "\n";
        return 2;
    }
    if (!ensureStage3Dirs(opt.output_dir, err)) {
        std::cerr << "[stage3][ERR] " << err << "\n";
        return 2;
    }

    const std::string report_path = pathJoin(pathJoin(opt.output_dir, "reports"), "stage3_sar_scene_report.md");
    if (!writeStage3PlanningReport(report_path, opt, err)) {
        std::cerr << "[stage3][ERR] " << err << "\n";
        return 3;
    }

    const std::string cfg_path = pathJoin(pathJoin(opt.output_dir, "config"), "stage3_config.example.json");
    if (!writeDefaultStage3Config(cfg_path, err)) {
        std::cerr << "[stage3][ERR] " << err << "\n";
        return 4;
    }

    std::cout << "[stage3] SAR-scene simulator scaffold ready. report=" << report_path << "\n";
    std::cout << "[stage3] Full SAR decoding, scatterer extraction, and DDC forward generation are tracked in simulator/docs/stage3_sar_scene_sim_implementation.md\n";
    return 0;
}
