#include "sar_georef.h"

namespace gmti {
namespace stage3 {

GroundPoint bilinearImageToGround(const ImageCornerGeoref &g,
                                  int width,
                                  int height,
                                  int row,
                                  int col)
{
    const double u = (width <= 1) ? 0.0 : static_cast<double>(col) / static_cast<double>(width - 1);
    const double v = (height <= 1) ? 0.0 : static_cast<double>(row) / static_cast<double>(height - 1);
    GroundPoint p;
    p.x_m = (1.0 - u) * (1.0 - v) * g.top_left.x_m + u * (1.0 - v) * g.top_right.x_m +
            (1.0 - u) * v * g.bottom_left.x_m + u * v * g.bottom_right.x_m;
    p.y_m = (1.0 - u) * (1.0 - v) * g.top_left.y_m + u * (1.0 - v) * g.top_right.y_m +
            (1.0 - u) * v * g.bottom_left.y_m + u * v * g.bottom_right.y_m;
    p.z_m = g.ground_z_m;
    return p;
}

bool loadImageCornerGeoref(const std::string &, ImageCornerGeoref &, std::string &err)
{
    err = "external georef JSON parsing is not implemented yet.";
    return false;
}

} // namespace stage3
} // namespace gmti
