#include "dbs/NewProtocolLayout.hpp"

#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include <cstdlib>

namespace {

std::vector<uint8_t> readAt(const std::string &path, std::size_t off, std::size_t n)
{
    std::ifstream in(path.c_str(), std::ios::binary);
    if (!in) {
        throw std::runtime_error("failed to open input file: " + path);
    }
    in.seekg(static_cast<std::streamoff>(off), std::ios::beg);
    std::vector<uint8_t> buf(n);
    in.read(reinterpret_cast<char *>(buf.data()), static_cast<std::streamsize>(buf.size()));
    if (in.gcount() != static_cast<std::streamsize>(buf.size())) {
        throw std::runtime_error("failed to read requested bytes");
    }
    return buf;
}

std::size_t fileSize(const std::string &path)
{
    std::ifstream in(path.c_str(), std::ios::binary | std::ios::ate);
    if (!in) {
        throw std::runtime_error("failed to open input file: " + path);
    }
    const std::streamoff size = in.tellg();
    if (size < 0) {
        throw std::runtime_error("failed to get file size");
    }
    return static_cast<std::size_t>(size);
}

void printHeader(const std::string &path, std::size_t packet_bytes, std::size_t pulse_idx)
{
    const std::vector<uint8_t> hdr =
        readAt(path, pulse_idx * packet_bytes, gmti::new_protocol::kHeaderBytes);
    const gmti::new_protocol::HeaderSample s =
        gmti::new_protocol::readHeaderSample(hdr.data());
    std::cout << pulse_idx << ","
              << std::setprecision(12)
              << s.utc << ","
              << s.lat_deg << ","
              << s.lon_deg << ","
              << s.height_m << ","
              << s.vn_mps << ","
              << s.ve_mps << ","
              << s.vd_mps << ","
              << s.theta_cmd_deg << "\n";
}

} // namespace

int main(int argc, char **argv)
{
    if (argc < 2 || argc > 4) {
        std::cerr << "usage: dump_newprotocol_header <stage2_statistical_newprotocol.bin> [channel_count] [iq_data_type]\n";
        return 2;
    }
    try {
        const std::string path = argv[1];
        const std::size_t channel_count = (argc >= 3)
            ? static_cast<std::size_t>(std::max(1, std::atoi(argv[2])))
            : 2U;
        const std::string iq_data_type = (argc >= 4) ? argv[3] : "float32";
        const std::size_t sample_bytes =
            gmti::new_protocol::sampleBytes(channel_count, iq_data_type);
        const std::vector<uint8_t> first =
            readAt(path, 0, gmti::new_protocol::kHeaderBytes);
        const std::size_t packet_bytes =
            static_cast<std::size_t>(gmti::new_protocol::loadU32LE(
                first.data() + gmti::new_protocol::kOffPrtLen));
        if (packet_bytes < gmti::new_protocol::kHeaderBytes ||
            ((packet_bytes - gmti::new_protocol::kHeaderBytes) %
             sample_bytes) != 0U) {
            throw std::runtime_error("invalid prt_len in first header");
        }
        const std::size_t size = fileSize(path);
        if (size == 0 || size % packet_bytes != 0U) {
            throw std::runtime_error("file size is not an integer multiple of prt_len");
        }
        const std::size_t total = size / packet_bytes;
        std::cout << "# packet_bytes=" << packet_bytes
                  << " channel_count=" << channel_count
                  << " iq_data_type=" << iq_data_type
                  << " samples_per_prt="
                  << ((packet_bytes - gmti::new_protocol::kHeaderBytes) /
                      sample_bytes)
                  << " total_prt=" << total << "\n";
        std::cout << "pulse_idx,utc,lat,lon,h,vn,ve,vd,theta_cmd_deg\n";
        printHeader(path, packet_bytes, 0);
        printHeader(path, packet_bytes, total / 2);
        printHeader(path, packet_bytes, total - 1);
    } catch (const std::exception &e) {
        std::cerr << "[ERR] " << e.what() << "\n";
        return 1;
    }
    return 0;
}
