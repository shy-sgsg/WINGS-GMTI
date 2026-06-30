#include "RingBuffer.h"
#include <iostream>
#include <chrono>
#include <thread>
#include <memory>
#include <cstring>
#include <cstdint>
#include <unistd.h>

int main() {
    const size_t chunk_size = 1ULL * 1024 * 1024;   // 1 MB

    try {
        RingBufferManager* manager = nullptr;
        while (true) {
            try {
                manager = new RingBufferManager("/my_ring_buffer", RingBufferManager::RECEIVER);
                break;
            } catch (const std::exception&) {
                std::cerr << "Waiting for sender to create shared memory...\n";
                std::this_thread::sleep_for(std::chrono::seconds(1));
            }
        }
        std::unique_ptr<RingBufferManager> mgr(manager);
        RingBuffer& ring = mgr->get_ring();

        char* buffer = new char[chunk_size];
        std::cout << "Receiver started, reading " << chunk_size / (1024*1024) << " MB per chunk\n";

        uint64_t total_read = 0;
        auto start_time = std::chrono::steady_clock::now();
        auto last_print_time = start_time;
        uint64_t bytes_since_last_print = 0;

        while (true) {
            size_t io_len = chunk_size;
            bool ok = ring.read(buffer, io_len);
            if (!ok) {
                std::cerr << "Read failed\n";
                break;
            }

            total_read += io_len;
            bytes_since_last_print += io_len;

            // ----- -----
            auto now = std::chrono::steady_clock::now();
            double sec_since_last = std::chrono::duration<double>(now - last_print_time).count();
            if (sec_since_last >= 1.0) {
                double current_speed = bytes_since_last_print / sec_since_last;
                std::cout << "[Receiver] Current speed: " << current_speed / 1e9 << " GB/s"
                          << "  (total: " << total_read / 1e9 << " GB)\n";
                last_print_time = now;
                bytes_since_last_print = 0;
            }
        }

        delete[] buffer;
    } catch (const std::exception& e) {
        std::cerr << "Exception: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
