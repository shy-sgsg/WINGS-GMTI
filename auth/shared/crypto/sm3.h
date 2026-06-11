#pragma once
#include <cstdint>
#include <string>
#include <vector>

namespace crypto {

class SM3 {
public:
    static constexpr size_t DIGEST_SIZE = 32;

    SM3();
    void update(const uint8_t* data, size_t len);
    void update(const std::string& s);
    void final(uint8_t digest[DIGEST_SIZE]);

    // 一次性计算
    static void hash(const uint8_t* data, size_t len, uint8_t digest[DIGEST_SIZE]);
    static void hash(const std::string& s, uint8_t digest[DIGEST_SIZE]);
    static std::vector<uint8_t> hash(const std::vector<uint8_t>& data);

private:
    uint32_t state_[8];
    uint8_t  buf_[64];
    uint64_t total_;
    size_t   buf_len_;

    void compress(const uint8_t block[64]);
    void reset();
};

} // namespace crypto
