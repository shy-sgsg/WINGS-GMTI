#include "sar_radiometric_normalizer.h"

#include <algorithm>
#include <cmath>

namespace gmti {
namespace stage3 {

bool normalizeSarRadiometry(const SarImage &input,
                            const RadiometricNormalizeConfig &,
                            SarImage &normalized,
                            std::string &err)
{
    if (input.value.empty() || input.width <= 0 || input.height <= 0) {
        err = "empty SAR image";
        return false;
    }
    normalized = input;
    float maxv = 0.0f;
    for (size_t i = 0; i < input.value.size(); ++i) maxv = std::max(maxv, input.value[i]);
    if (maxv <= 0.0f || !std::isfinite(maxv)) {
        err = "invalid SAR intensity range";
        return false;
    }
    for (size_t i = 0; i < normalized.value.size(); ++i) normalized.value[i] = std::max(0.0f, normalized.value[i] / maxv);
    return true;
}

} // namespace stage3
} // namespace gmti
