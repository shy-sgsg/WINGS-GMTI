#pragma once

#include <array>
#include <cstdint>

namespace security {

// 正式成像流程的七个检测点。枚举顺序也用于流程顺序检查。
enum class CheckpointId : std::uint8_t {
    ProgramStartup = 0,
    SarParameterInit = 1,
    ImagingKernel = 2,
    RangeCmp = 3,
    MoCo = 4,
    AzComp = 5,
    PGA = 6,
};

// 七个 checkpoint 的位图掩码；PGA 单独作为强制刷新要求。
inline constexpr std::uint32_t kAllCheckpointMask = 0x7fu;
inline constexpr std::uint32_t kPgaCheckpointMask = 0x40u;

// 解码后的 checkpoint 文本只在 gate 内部临时使用，避免业务层散落明文字符串。
struct DecodedCheckpoint {
    std::array<char, 32> text{};

    const char* c_str() const {
        return text.data();
    }
};

DecodedCheckpoint decode_checkpoint(CheckpointId id);

std::uint32_t checkpoint_bit(CheckpointId id);
std::uint8_t checkpoint_order(CheckpointId id);
bool checkpoint_requires_refresh(CheckpointId id);

} // namespace security
