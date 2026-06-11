#include "binding.h"
#include "binding_config.h"

#include "../../shared/crypto/sm3.h"

#include <algorithm>
#include <cctype>
#include <iomanip>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace hwbind {

namespace {

// 默认校验器会在多个 checkpoint 之间共享缓存结果，
// 用全局锁保证并发触发 checkpoint 时状态一致。
std::mutex g_verifier_mutex;

std::string to_hex(const std::vector<uint8_t>& bytes) {
    std::ostringstream oss;
    oss << std::hex << std::setfill('0');
    for (uint8_t b : bytes) {
        oss << std::setw(2) << static_cast<int>(b);
    }
    return oss.str();
}

std::string trim(std::string s) {
    const auto begin = s.find_first_not_of(" \t\r\n");
    if (begin == std::string::npos) {
        return "";
    }
    const auto end = s.find_last_not_of(" \t\r\n");
    return s.substr(begin, end - begin + 1);
}

std::string normalize_hex(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) {
        return static_cast<char>(std::tolower(c));
    });
    return s;
}

// 摘要比较避免使用普通字符串早停比较，减少比较过程泄露差异位置。
bool constant_time_equal(const std::string& a, const std::string& b) {
    if (a.size() != b.size()) {
        return false;
    }
    unsigned char diff = 0;
    for (size_t i = 0; i < a.size(); ++i) {
        diff |= static_cast<unsigned char>(a[i] ^ b[i]);
    }
    return diff == 0;
}

// hwbind_v3 是严格三字段策略：字段缺失直接拒绝生成绑定画像。
std::string require_token(const std::string& name, const std::string& value) {
    const std::string token = trim(value);
    if (token.empty()) {
        throw std::runtime_error("missing hardware identity field: " + name);
    }
    return token;
}

// 绑定串格式是摘要输入的一部分，修改 schema 或字段顺序都会改变摘要。
std::string build_canonical_hardware(FingerprintProfile& profile) {
    profile.cpu_component = require_token("cpu.serial", profile.hardware.cpu_serial);
    profile.disk_component = require_token("disk.serial", profile.hardware.disk_serial);
    profile.board_component = require_token("system.machine_id", profile.hardware.board_serial);

    return "schema=hwbind_v3|cpu=" + profile.cpu_component +
           "|disk=" + profile.disk_component +
           "|machine_id=" + profile.board_component;
}

// 真正的一次校验：采集硬件画像、计算当前摘要、和编译期目标摘要比较。
CheckResult compute_check_result(const std::string& soft_type,
                                 const std::string& soft_version,
                                 const std::string& expected_digest_hex,
                                 std::string checkpoint) {
    CheckResult result;
    result.checkpoint = std::move(checkpoint);
    result.expected_digest_hex = normalize_hex(expected_digest_hex);
    result.configured = !result.expected_digest_hex.empty() &&
                        result.expected_digest_hex.size() == 64 &&
                        std::all_of(result.expected_digest_hex.begin(),
                                    result.expected_digest_hex.end(),
                                    [](unsigned char c) { return std::isxdigit(c) != 0; });
    result.fingerprint = collect_fingerprint_profile();
    result.actual_digest_hex = digest_hex_for_profile(result.fingerprint, soft_type, soft_version);
    if (!result.configured) {
        return result;
    }
    result.matched = constant_time_equal(result.actual_digest_hex, result.expected_digest_hex);
    return result;
}

} // namespace

FingerprintProfile build_fingerprint_profile(const client::HardwareInfo& hardware) {
    FingerprintProfile profile;
    profile.hardware = hardware;
    profile.raw_serialized_hardware = client::serialize(hardware);
    profile.canonical_hardware = build_canonical_hardware(profile);
    // 当前版本没有弱绑定 fallback；能构建成功就认为三字段完整。
    profile.strength = FingerprintStrength::STRONG;
    profile.stability_score = 100;
    profile.used_machine_id_fallback = false;
    return profile;
}

FingerprintProfile collect_fingerprint_profile() {
    return build_fingerprint_profile(client::collect_hardware_info());
}

std::string build_binding_material(const std::string& canonical_hardware,
                                   const std::string& soft_type,
                                   const std::string& soft_version) {
    // 软件类型和绑定版本参与摘要，避免同一硬件摘要跨产品/策略复用。
    return canonical_hardware + "|" + soft_type + "|" + soft_version;
}

std::string digest_hex_for_binding_hardware(const std::string& canonical_hardware,
                                            const std::string& soft_type,
                                            const std::string& soft_version) {
    const std::string material = build_binding_material(canonical_hardware, soft_type, soft_version);
    const std::vector<uint8_t> data(material.begin(), material.end());
    return to_hex(crypto::SM3::hash(data));
}

std::string digest_hex_for_profile(const FingerprintProfile& profile,
                                   const std::string& soft_type,
                                   const std::string& soft_version) {
    return digest_hex_for_binding_hardware(profile.canonical_hardware, soft_type, soft_version);
}

bool is_configured() {
    const std::string expected = normalize_hex(HW_BIND_EXPECTED_SM3_HEX);
    if (expected.size() != 64) {
        return false;
    }
    return std::all_of(expected.begin(), expected.end(), [](unsigned char c) {
        return std::isxdigit(c) != 0;
    });
}

RuntimeVerifier::RuntimeVerifier(std::string soft_type,
                                 std::string soft_version,
                                 std::string expected_digest_hex)
    : soft_type_(std::move(soft_type)),
      soft_version_(std::move(soft_version)),
      expected_digest_hex_(normalize_hex(std::move(expected_digest_hex))) {}

CheckResult RuntimeVerifier::verify(bool force_refresh, std::string checkpoint) {
    std::lock_guard<std::mutex> lock(g_verifier_mutex);
    if (!force_refresh && has_cached_result_) {
        // 缓存命中时只更新 checkpoint 标签，硬件摘要结果仍复用第一次校验。
        cached_result_.from_cache = true;
        if (!checkpoint.empty()) {
            cached_result_.checkpoint = std::move(checkpoint);
        }
        return cached_result_;
    }

    // force_refresh=true 时重新采集硬件并重新计算摘要，PGA 等后段检查点会用到。
    cached_result_ = compute_check_result(soft_type_, soft_version_, expected_digest_hex_, std::move(checkpoint));
    cached_result_.from_cache = false;
    has_cached_result_ = true;
    return cached_result_;
}

bool RuntimeVerifier::require_match(const std::string& checkpoint,
                                    bool force_refresh) {
    return verify(force_refresh, checkpoint).matched;
}

RuntimeVerifier& default_verifier() {
    static RuntimeVerifier verifier(
        HW_BIND_SOFT_TYPE,
        HW_BIND_SOFT_VERSION,
        HW_BIND_EXPECTED_SM3_HEX
    );
    return verifier;
}

CheckResult verify_current_machine(bool force_refresh,
                                   std::string checkpoint) {
    return default_verifier().verify(force_refresh, std::move(checkpoint));
}

bool require_current_machine(const std::string& checkpoint,
                             bool force_refresh) {
    return default_verifier().require_match(checkpoint, force_refresh);
}

const char* strength_to_string(FingerprintStrength strength) {
    switch (strength) {
        case FingerprintStrength::STRONG:
            return "strong";
        case FingerprintStrength::MEDIUM:
            return "medium";
        case FingerprintStrength::WEAK:
            return "weak";
    }
    return "unknown";
}

} // namespace hwbind
