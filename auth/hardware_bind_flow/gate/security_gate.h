#pragma once

#include "checkpoint_codec.h"

#include <cstdint>
#include <iostream>

namespace security {

enum class GateCode {
    Ok = 0,
    NotConfigured,
    HardwareMismatch,
    BadFlow,
};

// combined/deferred gate 的运行期状态。
// seen/pass/fail/refresh 使用位图记录七个 checkpoint 的出现和结果。
struct GateSnapshot {
    std::uint32_t seen_mask = 0;
    std::uint32_t pass_mask = 0;
    std::uint32_t fail_mask = 0;
    std::uint32_t refresh_mask = 0;
    std::uint8_t highest_order = 0;
    bool has_checkpoint = false;
    bool flow_error = false;
    GateCode last_code = GateCode::Ok;
};

// 新任务或新流程开始前调用，清空上一轮 checkpoint 状态。
void reset_gate();

// 触发一个 checkpoint：内部解码名称、调用 real hwbind、记录状态。
GateCode checkpoint(CheckpointId id);

// 综合判定：默认要求七个点都出现、都通过，并且 PGA 完成强制刷新。
bool final_decision(std::uint32_t required_mask = kAllCheckpointMask,
                    std::uint32_t required_refresh_mask = kPgaCheckpointMask);

// 用于日志或测试观察当前 gate 状态，不改变内部状态。
GateSnapshot gate_snapshot();

const char* gate_code_to_string(GateCode code);

} // namespace security


bool run_checkpoint(security::CheckpointId id, const char* name);

