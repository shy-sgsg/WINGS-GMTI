#pragma once
#include <string>

namespace client {

struct HardwareInfo {
    std::string cpu_serial;
    std::string disk_serial;
    std::string board_serial;
    std::string cpu_source;
    std::string disk_source;
    std::string board_source;
};

// 采集当前机器硬件信息（需要 root 或适当权限）
HardwareInfo collect_hardware_info();

// 序列化为单行字符串，用于传输和密钥派生
std::string serialize(const HardwareInfo& hw);

// 从序列化字符串解析
HardwareInfo deserialize(const std::string& s);

} // namespace client
