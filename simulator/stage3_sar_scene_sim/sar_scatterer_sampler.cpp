#include "sar_scatterer_sampler.h"

namespace gmti {
namespace stage3 {

bool limitScattererCount(int max_scatterers,
                         std::vector<SarScatterer> &scatterers,
                         std::string &)
{
    if (max_scatterers >= 0 && static_cast<int>(scatterers.size()) > max_scatterers) {
        scatterers.resize(static_cast<size_t>(max_scatterers));
    }
    return true;
}

} // namespace stage3
} // namespace gmti
