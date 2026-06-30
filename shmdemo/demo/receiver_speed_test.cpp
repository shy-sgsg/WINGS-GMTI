#include "RingBuffer.h"
#include <iostream>
#include <chrono>
#include <memory>
#include <signal.h>
#include <atomic>
#include <thread>

std::atomic<bool> running{true};
void signal_handler(int) { running = false; }

int main() {
    const size_t chunk_size = 4ULL * 1024 * 1024;   // 4 MB

    try {
        RingBufferManager* manager = nullptr;
        while (true) {
            try {
                manager = new RingBufferManager("/my_ring_buffer", RingBufferManager::RECEIVER);
                break;
            } catch (const std::exception&) {
                std::cerr << "Waiting for sender...\n";
                std::this_thread::sleep_for(std::chrono::seconds(1));
            }
        }
        std::unique_ptr<RingBufferManager> mgr(manager);
        RingBuffer& ring = mgr->get_ring();

        char* buffer = new char[chunk_size];
        std::cout << "Receiver running in speed test mode. Press Ctrl+C to stop.\n";
        signal(SIGINT, signal_handler);

        uint64_t total_read = 0;
        auto start_time = std::chrono::steady_clock::now();

        while (running) {
            size_t io_len = chunk_size;
            if (ring.try_read(buffer, io_len)) {
                total_read += io_len;

                auto now = std::chrono::steady_clock::now();
                static auto last_print = start_time;
                if (std::chrono::duration<double>(now - last_print).count() >= 1.0) {
                    double elapsed = std::chrono::duration<double>(now - start_time).count();
                    double speed = total_read / elapsed;
                    std::cout << "Received: " << total_read / 1e9 << " GB, "
                              << "speed: " << speed / 1e9 << " GB/s\n";
                    last_print = now;
                }
            } else {
                std::this_thread::sleep_for(std::chrono::microseconds(100));
            }
        }

        auto end_time = std::chrono::steady_clock::now();
        double elapsed = std::chrono::duration<double>(end_time - start_time).count();
        double avg_speed = total_read / elapsed;
        std::cout << "\nTotal received: " << total_read / 1e9 << " GB\n"
                  << "Elapsed: " << elapsed << " s\n"
                  << "Average speed: " << avg_speed / 1e9 << " GB/s\n";

        delete[] buffer;
    } catch (const std::exception& e) {
        std::cerr << "Exception: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}
