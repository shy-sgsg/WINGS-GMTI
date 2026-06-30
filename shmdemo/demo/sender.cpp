#include "RingBuffer.h"
#include <iostream>
#include <chrono>
#include <thread>
#include <cstring>
#include <cstdint>
#include <unistd.h>
#include <csignal>
#include <atomic>

std::atomic<bool> running{true};

void signal_handler(int sig) {
    if (sig == SIGINT || sig == SIGTERM) {
        running = false;
    }
}

int main() {
    // 注册信号处理函数
    std::signal(SIGINT, signal_handler);
    std::signal(SIGTERM, signal_handler);

    const size_t buffer_size = 64ULL * 1024 * 1024;   // 64 MB
    const size_t chunk_size  = 1ULL * 1024 * 1024;    // 1 MB
    const double target_speed = 1.8 * 1024 * 1024 * 1024; // 1.8 GB/s

    try {
        RingBufferManager manager("/my_ring_buffer", RingBufferManager::SENDER, buffer_size);
        RingBuffer& ring = manager.get_ring();

        char* data = new char[chunk_size];
        memset(data, 'A', chunk_size);

        std::cout << "Sender started.\n"
                  << "  Chunk size : " << chunk_size / (1024*1024) << " MB\n"
                  << "  Target speed: " << target_speed / 1e9 << " GB/s\n";

        uint64_t total_sent = 0;
        auto start_time = std::chrono::steady_clock::now();
        auto last_print_time = start_time;
        uint64_t bytes_since_last_print = 0;

        while (running) {
            auto before = std::chrono::steady_clock::now();

            bool ok = ring.write(data, chunk_size);
            if (!ok) {
                std::cerr << "Write failed\n";
                break;
            }

            auto after = std::chrono::steady_clock::now();
            total_sent += chunk_size;
            bytes_since_last_print += chunk_size;

            // 
            auto elapsed_us = std::chrono::duration_cast<std::chrono::microseconds>(after - before).count();
            double expected_us = chunk_size * 1e6 / target_speed;
            if (elapsed_us < expected_us) {
                std::this_thread::sleep_for(
                    std::chrono::microseconds(static_cast<long long>(expected_us - elapsed_us))
                );
            }

            // ----- -----
            auto now = std::chrono::steady_clock::now();
            double sec_since_last = std::chrono::duration<double>(now - last_print_time).count();
            if (sec_since_last >= 1.0) {
                double current_speed = bytes_since_last_print / sec_since_last;
                std::cout << "[Sender] Current speed: " << current_speed / 1e9 << " GB/s"
                          << "  (total: " << total_sent / 1e9 << " GB)\n";
                // 
                last_print_time = now;
                bytes_since_last_print = 0;
            }
        }

        delete[] data;
    } catch (const std::exception& e) {
        std::cerr << "Exception: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
