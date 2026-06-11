#include "checkpoint_codec.h"

#include <cstddef>

namespace security {

namespace {

// 使用 volatile 防止编译器把 XOR 解码完全常量折叠回明文字符串。
volatile std::uint8_t g_decode_key = 0x5au;

// 以下数组是 checkpoint 明文与 g_decode_key 异或后的结果。
// 这样二进制中不会直接出现 program_startup/range_cmp 等明文。
constexpr std::array<std::uint8_t, 16> kProgramStartup = {
    0x2a, 0x28, 0x35, 0x3d, 0x28, 0x3b, 0x37, 0x05,
    0x29, 0x2e, 0x3b, 0x28, 0x2e, 0x2f, 0x2a, 0x5a
};

constexpr std::array<std::uint8_t, 19> kSarParameterInit = {
    0x29, 0x3b, 0x28, 0x05, 0x2a, 0x3b, 0x28, 0x3b,
    0x37, 0x3f, 0x2e, 0x3f, 0x28, 0x05, 0x33, 0x34,
    0x33, 0x2e, 0x5a
};

constexpr std::array<std::uint8_t, 15> kImagingKernel = {
    0x33, 0x37, 0x3b, 0x3d, 0x33, 0x34, 0x3d, 0x05,
    0x31, 0x3f, 0x28, 0x34, 0x3f, 0x36, 0x5a
};

constexpr std::array<std::uint8_t, 10> kRangeCmp = {
    0x28, 0x3b, 0x34, 0x3d, 0x3f, 0x05, 0x39, 0x37,
    0x2a, 0x5a
};

constexpr std::array<std::uint8_t, 5> kMoCo = {
    0x37, 0x35, 0x39, 0x35, 0x5a
};

constexpr std::array<std::uint8_t, 8> kAzComp = {
    0x3b, 0x20, 0x05, 0x39, 0x35, 0x37, 0x2a, 0x5a
};

constexpr std::array<std::uint8_t, 4> kPga = {
    0x2a, 0x3d, 0x3b, 0x5a
};

template <std::size_t N>
DecodedCheckpoint decode_bytes(const std::array<std::uint8_t, N>& encoded) {
    DecodedCheckpoint decoded;
    for (std::size_t i = 0; i < N && i < decoded.text.size(); ++i) {
        decoded.text[i] = static_cast<char>(encoded[i] ^ g_decode_key);
    }
    return decoded;
}

} // namespace

DecodedCheckpoint decode_checkpoint(CheckpointId id) {
    switch (id) {
        case CheckpointId::ProgramStartup:
            return decode_bytes(kProgramStartup);
        case CheckpointId::SarParameterInit:
            return decode_bytes(kSarParameterInit);
        case CheckpointId::ImagingKernel:
            return decode_bytes(kImagingKernel);
        case CheckpointId::RangeCmp:
            return decode_bytes(kRangeCmp);
        case CheckpointId::MoCo:
            return decode_bytes(kMoCo);
        case CheckpointId::AzComp:
            return decode_bytes(kAzComp);
        case CheckpointId::PGA:
            return decode_bytes(kPga);
    }
    return {};
}

std::uint32_t checkpoint_bit(CheckpointId id) {
    return 1u << checkpoint_order(id);
}

std::uint8_t checkpoint_order(CheckpointId id) {
    return static_cast<std::uint8_t>(id);
}

bool checkpoint_requires_refresh(CheckpointId id) {
    // ProgramStartup\PGA 检查点重新采集硬件，避免整条流程只依赖启动时缓存。
    return id == CheckpointId::ProgramStartup || id == CheckpointId::PGA;
}

} // namespace security
