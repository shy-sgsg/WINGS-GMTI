#include "stage3_validator.h"

#include <fstream>

namespace gmti {
namespace stage3 {

bool writeStage3PlanningReport(const std::string &path,
                               const Stage3RunOptions &opt,
                               std::string &err)
{
    std::ofstream md(path.c_str());
    if (!md) {
        err = "failed to open " + path;
        return false;
    }
    md << "# Stage 3 SAR scene simulation report\n\n";
    md << "This scaffold records the adjusted third-stage plan: high-resolution SAR image scatterer inversion plus new-system 61-beam dual-channel DDC forward generation.\n\n";
    md << "## Requested run\n\n";
    md << "- stage3_config: " << opt.stage3_config << "\n";
    md << "- sar_image: " << opt.sar_image << "\n";
    md << "- georef_path: " << opt.georef_path << "\n";
    md << "- scatterer_csv: " << opt.scatterer_csv << "\n";
    md << "- output_dir: " << opt.output_dir << "\n";
    md << "- scene_mode: " << opt.scene_mode << "\n";
    md << "- extract_only: " << (opt.extract_only ? "true" : "false") << "\n";
    md << "- forward_only: " << (opt.forward_only ? "true" : "false") << "\n";
    md << "- target_enabled: " << (opt.target_enabled ? "true" : "false") << "\n\n";
    md << "## Required limitations to retain in final validation\n\n";
    md << "1. Intensity-only SAR images cannot recover true coherent phase; random phase is a statistical approximation.\n";
    md << "2. SAR imaging geometry differs from the new 61-beam scan geometry, so extracted scatterers are a reflectivity-scene estimate, not raw echo inversion.\n";
    md << "3. Georeferencing accuracy directly controls spatial correctness.\n";
    md << "4. More scatterers improve background realism but increase compute cost quickly.\n";
    md << "5. The first platform model is ideal_straight; navigation disturbance and attitude error are future extensions.\n";
    md << "6. The first beam model is Gaussian; measured antenna patterns can replace it later.\n";
    md << "7. The first dual-channel phase model is baseline_approx; a baseline-vector model is the next upgrade.\n";
    return true;
}

} // namespace stage3
} // namespace gmti
