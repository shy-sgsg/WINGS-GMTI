#include <iostream>

extern "C" bool gmtiUploadTrigLutToDevice(const float*,
                                          const float*,
                                          const float*,
                                          const float*,
                                          int)
{
    std::cerr << "[TRIG-LUT][WARN] device upload requested in a host-only target" << std::endl;
    return false;
}

extern "C" bool gmtiSetDeviceTrigMode(int)
{
    return true;
}
