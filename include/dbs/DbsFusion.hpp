#ifndef DBS_FUSION_HPP
#define DBS_FUSION_HPP

#include "DbsFusionTypes.hpp"

bool estimateBeamPointingBiasByCenterBeams(const std::vector<int> &periodList,
                                           const std::vector<FusionBeamMeta> &beamMeta,
                                           double &biasDeg);

bool applyBeamPointingBiasToFusionContext(double biasDeg,
                                          FusionGroupContext &ctx);

bool relocateFusionDetections(const FusionGroupContext &ctx,
                              const Config &cfg,
                              std::vector<GMTIOutput> &results);
bool runDbsFusionImaging(const FusionGroupContext &ctx,
                         const Config &cfg,
                         bool useGpu = true);

#endif // DBS_FUSION_HPP
