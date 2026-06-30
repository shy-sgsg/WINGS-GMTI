#include "sar_image_reader.h"

namespace gmti {
namespace stage3 {

bool readSarImage(const std::string &, const std::string &, SarImage &, std::string &err)
{
    err = "SAR image decoding is not implemented yet; planned backends are GeoTIFF and PNG/TIF with external georef.";
    return false;
}

} // namespace stage3
} // namespace gmti
