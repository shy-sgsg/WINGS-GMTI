#pragma once

#include "sar_image_reader.h"

#include <string>

namespace gmti {
namespace stage3 {

struct RadiometricNormalizeConfig {
    double p_low = 1.0;
    double p_high = 99.9;
    std::string image_value_type = "intensity";
};

bool normalizeSarRadiometry(const SarImage &input,
                            const RadiometricNormalizeConfig &cfg,
                            SarImage &normalized,
                            std::string &err);

} // namespace stage3
} // namespace gmti
