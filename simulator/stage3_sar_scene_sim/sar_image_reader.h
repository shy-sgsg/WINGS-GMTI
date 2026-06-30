#pragma once

#include <string>
#include <vector>

namespace gmti {
namespace stage3 {

struct SarImage {
    int width = 0;
    int height = 0;
    std::vector<float> value;
};

bool readSarImage(const std::string &path, const std::string &input_type, SarImage &image, std::string &err);

} // namespace stage3
} // namespace gmti
