#include "sm3.h"
#include <cstring>
#include <stdexcept>

namespace crypto {

// 常量
static const uint32_t T[64] = {
    0x79cc4519,0x79cc4519,0x79cc4519,0x79cc4519,0x79cc4519,0x79cc4519,0x79cc4519,0x79cc4519,
    0x79cc4519,0x79cc4519,0x79cc4519,0x79cc4519,0x79cc4519,0x79cc4519,0x79cc4519,0x79cc4519,
    0x7a879d8a,0x7a879d8a,0x7a879d8a,0x7a879d8a,0x7a879d8a,0x7a879d8a,0x7a879d8a,0x7a879d8a,
    0x7a879d8a,0x7a879d8a,0x7a879d8a,0x7a879d8a,0x7a879d8a,0x7a879d8a,0x7a879d8a,0x7a879d8a,
    0x7a879d8a,0x7a879d8a,0x7a879d8a,0x7a879d8a,0x7a879d8a,0x7a879d8a,0x7a879d8a,0x7a879d8a,
    0x7a879d8a,0x7a879d8a,0x7a879d8a,0x7a879d8a,0x7a879d8a,0x7a879d8a,0x7a879d8a,0x7a879d8a,
    0x7a879d8a,0x7a879d8a,0x7a879d8a,0x7a879d8a,0x7a879d8a,0x7a879d8a,0x7a879d8a,0x7a879d8a,
    0x7a879d8a,0x7a879d8a,0x7a879d8a,0x7a879d8a,0x7a879d8a,0x7a879d8a,0x7a879d8a,0x7a879d8a
};

static inline uint32_t rotl(uint32_t x, int n) {
    return (x << n) | (x >> (32 - n));
}
static inline uint32_t P0(uint32_t x) { return x ^ rotl(x,9)  ^ rotl(x,17); }
static inline uint32_t P1(uint32_t x) { return x ^ rotl(x,15) ^ rotl(x,23); }
static inline uint32_t FF(uint32_t x, uint32_t y, uint32_t z, int j) {
    return (j < 16) ? (x ^ y ^ z) : ((x & y) | (x & z) | (y & z));
}
static inline uint32_t GG(uint32_t x, uint32_t y, uint32_t z, int j) {
    return (j < 16) ? (x ^ y ^ z) : ((x & y) | (~x & z));
}
static inline uint32_t be32(const uint8_t* p) {
    return ((uint32_t)p[0]<<24)|((uint32_t)p[1]<<16)|((uint32_t)p[2]<<8)|(uint32_t)p[3];
}
static inline void store_be32(uint8_t* p, uint32_t v) {
    p[0]=(v>>24)&0xff; p[1]=(v>>16)&0xff; p[2]=(v>>8)&0xff; p[3]=v&0xff;
}

// ── 初始值 ────────────────────────────────────────────────────────────────
static const uint32_t IV[8] = {
    0x7380166f,0x4914b2b9,0x172442d7,0xda8a0600,
    0xa96f30bc,0x163138aa,0xe38dee4d,0xb0fb0e4e
};

SM3::SM3() { reset(); }

void SM3::reset() {
    memcpy(state_, IV, sizeof(IV));
    total_   = 0;
    buf_len_ = 0;
}

void SM3::compress(const uint8_t block[64]) {
    uint32_t W[68], W1[64];
    for (int i = 0; i < 16; ++i) W[i] = be32(block + i*4);
    for (int i = 16; i < 68; ++i)
        W[i] = P1(W[i-16] ^ W[i-9] ^ rotl(W[i-3],15)) ^ rotl(W[i-13],7) ^ W[i-6];
    for (int i = 0; i < 64; ++i)
        W1[i] = W[i] ^ W[i+4];

    uint32_t A=state_[0], B=state_[1], C=state_[2], D=state_[3],
             E=state_[4], F=state_[5], G=state_[6], H=state_[7];

    for (int j = 0; j < 64; ++j) {
        uint32_t Tj = (j<16) ? 0x79cc4519u : 0x7a879d8au;
        uint32_t SS1 = rotl(rotl(A,12) + E + rotl(Tj,j%32), 7);
        uint32_t SS2 = SS1 ^ rotl(A,12);
        uint32_t TT1 = FF(A,B,C,j) + D + SS2 + W1[j];
        uint32_t TT2 = GG(E,F,G,j) + H + SS1 + W[j];
        D = C; C = rotl(B,9); B = A; A = TT1;
        H = G; G = rotl(F,19); F = E; E = P0(TT2);
    }
    state_[0]^=A; state_[1]^=B; state_[2]^=C; state_[3]^=D;
    state_[4]^=E; state_[5]^=F; state_[6]^=G; state_[7]^=H;
}

void SM3::update(const uint8_t* data, size_t len) {
    total_ += (uint64_t)len * 8;
    size_t off = 0;
    if (buf_len_ > 0) {
        size_t need = 64 - buf_len_;
        if (len < need) {
            memcpy(buf_ + buf_len_, data, len);
            buf_len_ += len;
            return;
        }
        memcpy(buf_ + buf_len_, data, need);
        compress(buf_);
        buf_len_ = 0;
        off = need;
    }
    while (off + 64 <= len) {
        compress(data + off);
        off += 64;
    }
    if (off < len) {
        buf_len_ = len - off;
        memcpy(buf_, data + off, buf_len_);
    }
}

void SM3::update(const std::string& s) {
    update(reinterpret_cast<const uint8_t*>(s.data()), s.size());
}

void SM3::final(uint8_t digest[DIGEST_SIZE]) {
    buf_[buf_len_++] = 0x80;
    if (buf_len_ > 56) {
        memset(buf_ + buf_len_, 0, 64 - buf_len_);
        compress(buf_);
        buf_len_ = 0;
    }
    memset(buf_ + buf_len_, 0, 56 - buf_len_);
    // append bit-length (big-endian 64-bit)
    for (int i = 7; i >= 0; --i)
        buf_[56 + (7-i)] = (total_ >> (i*8)) & 0xff;
    compress(buf_);
    for (int i = 0; i < 8; ++i)
        store_be32(digest + i*4, state_[i]);
    reset();
}

void SM3::hash(const uint8_t* data, size_t len, uint8_t digest[DIGEST_SIZE]) {
    SM3 h; h.update(data, len); h.final(digest);
}
void SM3::hash(const std::string& s, uint8_t digest[DIGEST_SIZE]) {
    hash(reinterpret_cast<const uint8_t*>(s.data()), s.size(), digest);
}
std::vector<uint8_t> SM3::hash(const std::vector<uint8_t>& data) {
    std::vector<uint8_t> d(DIGEST_SIZE);
    hash(data.data(), data.size(), d.data());
    return d;
}

} // namespace crypto
