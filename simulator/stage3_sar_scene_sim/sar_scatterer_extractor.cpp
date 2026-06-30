#include "sar_scatterer_extractor.h"

#include <cmath>

namespace gmti {
namespace stage3 {

bool extractSarScatterers(const SarImage &normalized,
                          const ImageCornerGeoref &georef,
                          const Stage3Config &cfg,
                          std::vector<SarScatterer> &scatterers,
                          std::string &err)
{
    if (normalized.value.empty()) {
        err = "empty normalized SAR image";
        return false;
    }
    scatterers.clear();
    const int stride = std::max(1, static_cast<int>(std::sqrt(static_cast<double>(normalized.value.size()) /
                                                              std::max(1, cfg.extraction.max_scatterers))));
    for (int r = 0; r < normalized.height; r += stride) {
        for (int c = 0; c < normalized.width; c += stride) {
            const size_t idx = static_cast<size_t>(r) * static_cast<size_t>(normalized.width) + static_cast<size_t>(c);
            const double inorm = normalized.value[idx];
            if (inorm <= 0.0) continue;
            SarScatterer s;
            s.scatterer_id = static_cast<int>(scatterers.size());
            s.source_row = r;
            s.source_col = c;
            s.pos = bilinearImageToGround(georef, normalized.width, normalized.height, r, c);
            s.sar_intensity_raw = inorm;
            s.sar_intensity_norm = inorm;
            s.amplitude = cfg.extraction.amplitude_scale * std::pow(inorm, cfg.extraction.amplitude_gamma);
            scatterers.push_back(s);
            if (static_cast<int>(scatterers.size()) >= cfg.extraction.max_scatterers) return true;
        }
    }
    return true;
}

} // namespace stage3
} // namespace gmti
