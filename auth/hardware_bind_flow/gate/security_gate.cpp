#include "security_gate.h"

#include "../core/binding.h"

#include <exception>
#include <mutex>
#include <string>

bool run_checkpoint(security::CheckpointId id, const char* name)
{
    const security::GateCode code = security::checkpoint(id);
    std::cout << "[gate] " << name << ": "
              << security::gate_code_to_string(code) << "\n";
    return code == security::GateCode::Ok;
}

namespace security {

namespace {

std::mutex g_gate_mutex;
GateSnapshot g_state;

// 将一次 checkpoint 结果折叠进全局状态。
// 这里同时检查流程顺序，避免后面的检查点先出现再回跳前面的检查点。
void record_result(GateSnapshot& state,
                   CheckpointId id,
                   bool passed,
                   bool refreshed,
                   GateCode code) {
    const std::uint32_t bit = checkpoint_bit(id);
    const std::uint8_t order = checkpoint_order(id);

    if (state.has_checkpoint && order < state.highest_order) {
        state.flow_error = true;
        if (code == GateCode::Ok) {
            code = GateCode::BadFlow;
            passed = false;
        }
    }

    state.has_checkpoint = true;
    if (order > state.highest_order) {
        state.highest_order = order;
    }

    state.seen_mask |= bit;
    if (refreshed) {
        state.refresh_mask |= bit;
    }
    if (passed) {
        state.pass_mask |= bit;
    } else {
        state.fail_mask |= bit;
    }
    state.last_code = code;
}

} // namespace

void reset_gate() {
    std::lock_guard<std::mutex> lock(g_gate_mutex);
    g_state = GateSnapshot{};
}

GateCode checkpoint(CheckpointId id) {
    const DecodedCheckpoint decoded = decode_checkpoint(id);
    const bool force_refresh = checkpoint_requires_refresh(id);

    GateCode code = GateCode::Ok;
    bool matched = false;
    try {
        // real hwbind 仍然只判断“当前机器是否匹配内置摘要”；
        // gate 层负责把多个 checkpoint 的结果组合成流程级判定。
        const hwbind::CheckResult result =
            hwbind::verify_current_machine(force_refresh, std::string(decoded.c_str()));
        if (!result.configured) {
            code = GateCode::NotConfigured;
        } else if (!result.matched) {
            code = GateCode::HardwareMismatch;
        }
        matched = code == GateCode::Ok;
    } catch (const std::exception&) {
        code = GateCode::HardwareMismatch;
    }

    std::lock_guard<std::mutex> lock(g_gate_mutex);
    record_result(g_state, id, matched, force_refresh, code);
    return g_state.last_code;
}

bool final_decision(std::uint32_t required_mask,
                    std::uint32_t required_refresh_mask) {
    std::lock_guard<std::mutex> lock(g_gate_mutex);
    // deferred 决策点：不关心单次调用位置是否立即失败，
    // 而是统一检查必经点、通过点、失败点和强制刷新点。
    return !g_state.flow_error &&
           (g_state.seen_mask & required_mask) == required_mask &&
           (g_state.pass_mask & required_mask) == required_mask &&
           (g_state.fail_mask & required_mask) == 0 &&
           (g_state.refresh_mask & required_refresh_mask) == required_refresh_mask;
}

GateSnapshot gate_snapshot() {
    std::lock_guard<std::mutex> lock(g_gate_mutex);
    return g_state;
}

const char* gate_code_to_string(GateCode code) {
    switch (code) {
        case GateCode::Ok:
            return "ok";
        case GateCode::NotConfigured:
            return "not_configured";
        case GateCode::HardwareMismatch:
            return "hardware_mismatch";
        case GateCode::BadFlow:
            return "bad_flow";
    }
    return "unknown";
}

} // namespace security
