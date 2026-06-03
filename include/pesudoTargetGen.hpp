#include <vector>
#include <cmath>
#include <ctime>
#include <chrono>
#include <geo/geoProj.hpp>
#include <config_structs.hpp> // 包含 Config / GMTIOutput::Result 声明

void pesudoTarget_Gen(GMTIOutput& resRef,
                      double lat_st, double lon_st,
                      double lat_ed, double lon_ed,
                      double v_mps, double L0_deg,
                    const Config& cfg);