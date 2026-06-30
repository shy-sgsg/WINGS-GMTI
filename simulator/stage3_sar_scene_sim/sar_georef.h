#pragma once

#include <string>

namespace gmti {
namespace stage3 {

struct GroundPoint {
    double x_m = 0.0;
    double y_m = 0.0;
    double z_m = 0.0;
};

struct ImageCornerGeoref {
    GroundPoint top_left;
    GroundPoint top_right;
    GroundPoint bottom_left;
    GroundPoint bottom_right;
    double ground_z_m = 0.0;
};

GroundPoint bilinearImageToGround(const ImageCornerGeoref &georef,
                                  int width,
                                  int height,
                                  int row,
                                  int col);

bool loadImageCornerGeoref(const std::string &path, ImageCornerGeoref &georef, std::string &err);

} // namespace stage3
} // namespace gmti
