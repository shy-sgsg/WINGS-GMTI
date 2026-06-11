#pragma once

// 这三个常量编译进程序。
// 实际交付时，通常流程是：
// 1. 在目标机器运行 hw_bind_client --collect
// 2. 用 hw_bind_tool 根据采集结果生成 64 位十六进制摘要
// 3. 把下面的 HW_BIND_EXPECTED_SM3_HEX 替换成目标值后重新编译

inline constexpr const char* HW_BIND_SOFT_TYPE = "MyApp";
inline constexpr const char* HW_BIND_SOFT_VERSION = "1.0.0";
inline constexpr const char* HW_BIND_EXPECTED_SM3_HEX = "53743eded24c191368208e5763da146627fd4b2855675e713d403cb7a361f9f7";
