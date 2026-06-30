#include "sar_scene_forward.h"

namespace gmti {
namespace stage3 {

bool forwardSarSceneToDdc(const Stage3Config &,
                          const std::vector<SarScatterer> &,
                          const std::string &,
                          std::string &err)
{
    err = "SAR scatterer DDC forward generation is not implemented yet; reuse stage2 LFM kernels with beam/range prefilters.";
    return false;
}

} // namespace stage3
} // namespace gmti
