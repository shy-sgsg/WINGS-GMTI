#include "RingBuffer.h"
#include <iostream>
#include <chrono>
#include <cstring>
#include <signal.h>
#include <atomic>

std::atomic<bool> running{true};

void signal_handler(int) { running = false; }

int main() {
    const size_t buffer_size = 64ULL * 1024 * 1024;   // 64 MB
    const size_t chunk_size  = 4ULL * 1024 * 1024;    // 4 MB

    try {
        RingBufferManager manager("/my_ring_buffer", RingBufferManager::SENDER, buffer_size);
        RingBuffer& ring = manager.get_ring();

        char* data = new char[chunk_size];
        memset(data, 'A', chunk_size);

        std::cout << "Sender running in speed test mode (no rate limit). Press Ctrl+C to stop.\n";

        signal(SIGINT, signal_handler);

        uint64_t total_sent = 0;
        auto start_time = std::chrono::steady_clock::now();

        while (running) {
            if (!ring.write(data, chunk_size)) {
                std::cerr << "Write failed (buffer full?)\n";
                break;
            }
            total_sent += chunk_size;

            // 
            auto now = std::chrono::steady_clock::now();
            static auto last_print = start_time;
            if (std::chrono::duration<double>(now - last_print).count() >= 1.0) {
                double elapsed = std::chrono::duration<double>(now - start_time).count();
                double speed = total_sent / elapsed;
                std::cout << "Sent: " << total_sent / 1e9 << " GB, "
                          << "speed: " << speed / 1e9 << " GB/s\n";
                last_print = now;
            }
        }

        auto end_time = std::chrono::steady_clock::now();
        double elapsed = std::chrono::duration<double>(end_time - start_time).count();
        double avg_speed = total_sent / elapsed;
        std::cout << "\nTotal sent: " << total_sent / 1e9 << " GB\n"
                  << "Elapsed: " << elapsed << " s\n"
                  << "Average speed: " << avg_speed / 1e9 << " GB/s\n";

        delete[] data;
    } catch (const std::exception& e) {
        std::cerr << "Exception: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}
