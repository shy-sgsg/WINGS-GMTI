#pragma once

#include "../../shared/hardware/hardware.h"

#include <string>
#include <vector>

namespace hwbind {

enum class FingerprintStrength {
    WEAK = 0,
    MEDIUM = 1,
    STRONG = 2,
};

// 一次硬件采集与规范化后的完整画像。
// 当前 hwbind_v3 策略要求 cpu/disk/machine_id 三个组件都可用。
struct FingerprintProfile {
    client::HardwareInfo hardware;
    std::string raw_serialized_hardware;
    std::string canonical_hardware;
    std::string cpu_component;
    std::string disk_component;
    std::string board_component;
    FingerprintStrength strength = FingerprintStrength::WEAK;
    int stability_score = 0;
    bool used_machine_id_fallback = false;
    std::vector<std::string> advisories;
};

// 单次运行时校验结果。gate 层会读取 matched/configured/from_cache
// 来判断当前 checkpoint 是否通过，以及是否复用了缓存。
struct CheckResult {
    bool configured = false;
    bool matched = false;
    bool from_cache = false;
    std::string checkpoint;
    FingerprintProfile fingerprint;
    std::string actual_digest_hex;
    std::string expected_digest_hex;
};

FingerprintProfile build_fingerprint_profile(const client::HardwareInfo& hardware);
FingerprintProfile collect_fingerprint_profile();

std::string build_binding_material(const std::string& canonical_hardware,
                                   const std::string& soft_type,
                                   const std::string& soft_version);

std::string digest_hex_for_binding_hardware(const std::string& canonical_hardware,
                                            const std::string& soft_type,
                                            const std::string& soft_version);

std::string digest_hex_for_profile(const FingerprintProfile& profile,
                                   const std::string& soft_type,
                                   const std::string& soft_version);

bool is_configured();

// 运行期校验器会缓存第一次硬件采集和摘要比较结果。
// 多个 checkpoint 默认复用缓存；需要重新采集时传 force_refresh=true。
class RuntimeVerifier {
public:
    RuntimeVerifier(std::string soft_type,
                    std::string soft_version,
                    std::string expected_digest_hex);

    CheckResult verify(bool force_refresh = false,
                       std::string checkpoint = {});

    bool require_match(const std::string& checkpoint,
                       bool force_refresh = false);

private:
    std::string soft_type_;
    std::string soft_version_;
    std::string expected_digest_hex_;
    bool has_cached_result_ = false;
    CheckResult cached_result_;
};

RuntimeVerifier& default_verifier();

CheckResult verify_current_machine(bool force_refresh = false,
                                   std::string checkpoint = {});

bool require_current_machine(const std::string& checkpoint,
                             bool force_refresh = false);

const char* strength_to_string(FingerprintStrength strength);

} // namespace hwbind
