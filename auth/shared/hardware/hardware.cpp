#include "hardware.h"

#include <algorithm>
#include <cstdio>
#include <fstream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace client {

namespace {

struct DetectedValue {
    std::string value;
    std::string source;
};

std::string trim(const std::string& s) {
    const auto b = s.find_first_not_of(" \t\r\n");
    if (b == std::string::npos) {
        return "";
    }
    const auto e = s.find_last_not_of(" \t\r\n");
    return s.substr(b, e - b + 1);
}

std::string read_first_line(const std::string& path) {
    std::ifstream f(path);
    if (!f) {
        return "";
    }
    std::string line;
    std::getline(f, line);
    return trim(line);
}

std::string exec_cmd(const std::string& cmd) {
    FILE* pipe = popen(cmd.c_str(), "r");
    if (!pipe) {
        return "";
    }

    char buf[256] = {};
    std::string result;
    while (fgets(buf, sizeof(buf), pipe)) {
        result += buf;
    }
    pclose(pipe);
    return result;
}

std::vector<std::string> split_lines(const std::string& text) {
    std::vector<std::string> lines;
    std::istringstream iss(text);
    std::string line;
    while (std::getline(iss, line)) {
        line = trim(line);
        if (!line.empty()) {
            lines.push_back(line);
        }
    }
    return lines;
}

std::string join(const std::vector<std::string>& values, const std::string& sep) {
    std::ostringstream oss;
    for (size_t i = 0; i < values.size(); ++i) {
        if (i != 0) {
            oss << sep;
        }
        oss << values[i];
    }
    return oss.str();
}

DetectedValue get_cpu_serial() {
    std::ifstream f("/proc/cpuinfo");
    if (!f) {
        return {"", "/proc/cpuinfo"};
    }

    std::string line;
    while (std::getline(f, line)) {
        if (line.find("Serial") == std::string::npos) {
            continue;
        }
        const auto p = line.find(':');
        if (p == std::string::npos) {
            continue;
        }
        const std::string value = trim(line.substr(p + 1));
        if (!value.empty()) {
            return {value, "/proc/cpuinfo"};
        }
    }
    return {"", "/proc/cpuinfo"};
}

DetectedValue get_disk_serial() {
    const std::string output = exec_cmd("lsblk -d -o serial 2>/dev/null");
    const std::vector<std::string> lines = split_lines(output);

    std::set<std::string> serials;
    for (const auto& line : lines) {
        if (line == "SERIAL") {
            continue;
        }
        serials.insert(line);
    }

    std::vector<std::string> sorted(serials.begin(), serials.end());
    return {join(sorted, ","), "lsblk"};
}

DetectedValue get_machine_id() {
    return {read_first_line("/etc/machine-id"), "/etc/machine-id"};
}

void require_detected(const DetectedValue& value, const std::string& name) {
    if (value.value.empty()) {
        throw std::runtime_error("missing hardware identity field: " + name +
                                 " (" + value.source + ")");
    }
}

} // namespace

HardwareInfo collect_hardware_info() {
    const DetectedValue cpu = get_cpu_serial();
    const DetectedValue disk = get_disk_serial();
    const DetectedValue machine = get_machine_id();

    require_detected(cpu, "cpu.serial");
    require_detected(disk, "disk.serial");
    require_detected(machine, "system.machine_id");

    HardwareInfo hw;
    hw.cpu_serial = cpu.value;
    hw.cpu_source = cpu.source;
    hw.disk_serial = disk.value;
    hw.disk_source = disk.source;
    hw.board_serial = machine.value;
    hw.board_source = machine.source;
    return hw;
}

std::string serialize(const HardwareInfo& hw) {
    // Keep the existing public format. In this implementation board_serial stores /etc/machine-id.
    return "cpu=" + hw.cpu_serial +
           "|disk=" + hw.disk_serial +
           "|board=" + hw.board_serial;
}

HardwareInfo deserialize(const std::string& s) {
    HardwareInfo hw;
    auto get_field = [&](const std::string& key) -> std::string {
        auto pos = s.find(key + "=");
        if (pos == std::string::npos) {
            return "";
        }
        pos += key.size() + 1;
        const auto end = s.find('|', pos);
        return s.substr(pos, (end == std::string::npos) ? std::string::npos : end - pos);
    };

    hw.cpu_serial = get_field("cpu");
    hw.disk_serial = get_field("disk");
    hw.board_serial = get_field("board");
    return hw;
}

} // namespace client
