// m2_pcie_load.cpp
// Linux C++17
//
// 功能：
//   1. 独立读线程：从 M.2 文件/块设备持续读取，可限速/不限速
//   2. 独立写线程：向 M.2 文件/块设备持续写入，可限速/不限速
//   3. 每秒输出实时读写速率
//
// 编译：
//   g++ -O2 -std=c++17 -pthread m2_pcie_load.cpp -o m2_pcie_load
//
// 注意：
//   - 默认使用 O_DIRECT，要求 block size 通常为 4096 的整数倍。
//   - 推荐写普通测试文件，不推荐直接写 /dev/nvme0n1。
//   - --read-rate 0 表示读不限速。
//   - --write-rate 0 表示写不限速。

#include <atomic>
#include <chrono>
#include <csignal>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <fcntl.h>
#include <iomanip>
#include <iostream>
#include <linux/fs.h>
#include <string>
#include <sys/ioctl.h>
#include <sys/stat.h>
#include <thread>
#include <unistd.h>

static std::atomic<bool> g_running{true};

static void signal_handler(int)
{
    g_running.store(false);
}

struct Config {
    std::string read_path;
    std::string write_path;

    double read_rate_mib = -1.0;   // -1 表示未启用读线程；0 表示不限速
    double write_rate_mib = -1.0;  // -1 表示未启用写线程；0 表示不限速

    uint64_t block_size = 1024ULL * 1024ULL;       // 默认 1 MiB
    uint64_t write_size = 16ULL * 1024ULL * 1024ULL * 1024ULL; // 默认 16 GiB
    int duration_sec = 60;                         // 0 表示一直运行
    bool use_direct = true;
    bool allow_raw_write = false;
};

struct Counters {
    std::atomic<uint64_t> read_bytes{0};
    std::atomic<uint64_t> write_bytes{0};
    std::atomic<uint64_t> read_ops{0};
    std::atomic<uint64_t> write_ops{0};
    std::atomic<uint64_t> read_errors{0};
    std::atomic<uint64_t> write_errors{0};
};

static void print_help(const char* prog)
{
    std::cout << R"(Usage:
  )" << prog << R"( [options]

Options:
  --read-path PATH          读取源，可以是普通文件，也可以是块设备，例如 /dev/nvme0n1
  --write-path PATH         写入目标，强烈建议使用普通测试文件，例如 /mnt/nvme/write_test.bin

  --read-rate MiBps         读速率，单位 MiB/s；0 表示不限速；不设置则不启动读线程
  --write-rate MiBps        写速率，单位 MiB/s；0 表示不限速；不设置则不启动写线程

  --bs SIZE                 单次 IO 块大小，默认 1M。支持 K/M/G 后缀，例如 4K、1M、8M
  --write-size SIZE         写测试文件大小，默认 16G。写到末尾后从头循环覆盖
  --duration SEC            测试时长，默认 60 秒；0 表示一直运行直到 Ctrl+C

  --direct 0|1              是否使用 O_DIRECT，默认 1
  --allow-raw-write         允许直接写块设备，危险，会破坏磁盘数据
  --help                    显示帮助

Examples:
  # 同时读 800 MiB/s，写 400 MiB/s，运行 300 秒
  )" << prog << R"( \
    --read-path /mnt/nvme/input.bin \
    --write-path /mnt/nvme/write_test.bin \
    --read-rate 800 \
    --write-rate 400 \
    --bs 1M \
    --write-size 64G \
    --duration 300

  # 读不限速，写 500 MiB/s
  )" << prog << R"( \
    --read-path /mnt/nvme/input.bin \
    --write-path /mnt/nvme/write_test.bin \
    --read-rate 0 \
    --write-rate 500

  # 只读，1 GiB/s
  )" << prog << R"( \
    --read-path /mnt/nvme/input.bin \
    --read-rate 1024

  # 只写，不限速
  )" << prog << R"( \
    --write-path /mnt/nvme/write_test.bin \
    --write-rate 0 \
    --write-size 64G

)";
}

static uint64_t parse_size(const std::string& s)
{
    if (s.empty()) {
        throw std::runtime_error("empty size string");
    }

    char suffix = s.back();
    double value = 0.0;

    if (suffix == 'K' || suffix == 'k' ||
        suffix == 'M' || suffix == 'm' ||
        suffix == 'G' || suffix == 'g') {
        value = std::stod(s.substr(0, s.size() - 1));
    } else {
        suffix = 0;
        value = std::stod(s);
    }

    double mul = 1.0;
    if (suffix == 'K' || suffix == 'k') mul = 1024.0;
    if (suffix == 'M' || suffix == 'm') mul = 1024.0 * 1024.0;
    if (suffix == 'G' || suffix == 'g') mul = 1024.0 * 1024.0 * 1024.0;

    if (value <= 0.0) {
        throw std::runtime_error("size must be positive");
    }

    return static_cast<uint64_t>(value * mul);
}

static bool get_fd_size(int fd, uint64_t& size_bytes, bool& is_block_device)
{
    struct stat st {};
    if (fstat(fd, &st) != 0) {
        return false;
    }

    is_block_device = S_ISBLK(st.st_mode);

    if (S_ISREG(st.st_mode)) {
        size_bytes = static_cast<uint64_t>(st.st_size);
        return size_bytes > 0;
    }

    if (S_ISBLK(st.st_mode)) {
        uint64_t bytes = 0;
        if (ioctl(fd, BLKGETSIZE64, &bytes) == 0 && bytes > 0) {
            size_bytes = bytes;
            return true;
        }
    }

    return false;
}

class RateLimiter {
public:
    explicit RateLimiter(double mib_per_sec)
    {
        if (mib_per_sec <= 0.0) {
            rate_bytes_per_sec_ = 0.0; // 0 表示不限速
        } else {
            rate_bytes_per_sec_ = mib_per_sec * 1024.0 * 1024.0;
        }
        start_ = std::chrono::steady_clock::now();
    }

    void add_and_wait(uint64_t bytes)
    {
        if (rate_bytes_per_sec_ <= 0.0) {
            return;
        }

        bytes_done_ += bytes;

        const double expected_sec =
            static_cast<double>(bytes_done_) / rate_bytes_per_sec_;

        auto now = std::chrono::steady_clock::now();
        const double elapsed_sec =
            std::chrono::duration<double>(now - start_).count();

        if (expected_sec > elapsed_sec) {
            const double sleep_sec = expected_sec - elapsed_sec;
            std::this_thread::sleep_for(
                std::chrono::duration<double>(sleep_sec));
        }
    }

private:
    double rate_bytes_per_sec_ = 0.0;
    uint64_t bytes_done_ = 0;
    std::chrono::steady_clock::time_point start_;
};

static void* alloc_aligned_buffer(uint64_t block_size)
{
    void* ptr = nullptr;
    const size_t alignment = 4096;

    int ret = posix_memalign(&ptr, alignment, block_size);
    if (ret != 0 || ptr == nullptr) {
        return nullptr;
    }

    std::memset(ptr, 0x5A, block_size);
    return ptr;
}

static void read_worker(const Config& cfg, Counters& cnt)
{
    int flags = O_RDONLY;
    if (cfg.use_direct) {
        flags |= O_DIRECT;
    }

    int fd = open(cfg.read_path.c_str(), flags);
    if (fd < 0) {
        std::cerr << "[READ] open failed: " << cfg.read_path
                  << ", errno=" << errno
                  << ", " << std::strerror(errno) << std::endl;
        g_running.store(false);
        return;
    }

    uint64_t source_size = 0;
    bool is_block = false;
    if (!get_fd_size(fd, source_size, is_block)) {
        std::cerr << "[READ] failed to get source size: " << cfg.read_path << std::endl;
        close(fd);
        g_running.store(false);
        return;
    }

    void* buf = alloc_aligned_buffer(cfg.block_size);
    if (!buf) {
        std::cerr << "[READ] aligned buffer allocation failed" << std::endl;
        close(fd);
        g_running.store(false);
        return;
    }

    std::cout << "[READ] path=" << cfg.read_path
              << ", size=" << static_cast<double>(source_size) / 1024.0 / 1024.0 / 1024.0
              << " GiB"
              << ", rate=" << cfg.read_rate_mib << " MiB/s"
              << ", direct=" << cfg.use_direct
              << std::endl;

    RateLimiter limiter(cfg.read_rate_mib);

    uint64_t offset = 0;

    while (g_running.load()) {
        if (source_size >= cfg.block_size) {
            if (offset + cfg.block_size > source_size) {
                offset = 0;
            }
        } else {
            offset = 0;
        }

        ssize_t n = pread(fd, buf, cfg.block_size, static_cast<off_t>(offset));
        if (n > 0) {
            cnt.read_bytes.fetch_add(static_cast<uint64_t>(n), std::memory_order_relaxed);
            cnt.read_ops.fetch_add(1, std::memory_order_relaxed);
            offset += static_cast<uint64_t>(n);
            limiter.add_and_wait(static_cast<uint64_t>(n));
        } else if (n == 0) {
            offset = 0;
        } else {
            cnt.read_errors.fetch_add(1, std::memory_order_relaxed);

            if (errno == EINTR || errno == EAGAIN) {
                continue;
            }

            std::cerr << "[READ] pread failed, errno=" << errno
                      << ", " << std::strerror(errno) << std::endl;
            break;
        }
    }

    free(buf);
    close(fd);
}

static void write_worker(const Config& cfg, Counters& cnt)
{
    int flags = O_RDWR | O_CREAT;
    if (cfg.use_direct) {
        flags |= O_DIRECT;
    }

    int fd = open(cfg.write_path.c_str(), flags, 0644);
    if (fd < 0) {
        std::cerr << "[WRITE] open failed: " << cfg.write_path
                  << ", errno=" << errno
                  << ", " << std::strerror(errno) << std::endl;
        g_running.store(false);
        return;
    }

    uint64_t target_size = 0;
    bool is_block = false;

    struct stat st {};
    if (fstat(fd, &st) != 0) {
        std::cerr << "[WRITE] fstat failed" << std::endl;
        close(fd);
        g_running.store(false);
        return;
    }

    is_block = S_ISBLK(st.st_mode);

    if (is_block) {
        if (!cfg.allow_raw_write) {
            std::cerr << "[WRITE] Refuse to write raw block device: "
                      << cfg.write_path << std::endl;
            std::cerr << "[WRITE] If you really want to destroy data on this device, "
                      << "add --allow-raw-write" << std::endl;
            close(fd);
            g_running.store(false);
            return;
        }

        uint64_t bytes = 0;
        if (ioctl(fd, BLKGETSIZE64, &bytes) != 0 || bytes == 0) {
            std::cerr << "[WRITE] failed to get block device size" << std::endl;
            close(fd);
            g_running.store(false);
            return;
        }
        target_size = bytes;
    } else {
        target_size = cfg.write_size;

        if (ftruncate(fd, static_cast<off_t>(target_size)) != 0) {
            std::cerr << "[WRITE] ftruncate failed: errno=" << errno
                      << ", " << std::strerror(errno) << std::endl;
            close(fd);
            g_running.store(false);
            return;
        }
    }

    if (target_size < cfg.block_size) {
        std::cerr << "[WRITE] target size is smaller than block size" << std::endl;
        close(fd);
        g_running.store(false);
        return;
    }

    void* buf = alloc_aligned_buffer(cfg.block_size);
    if (!buf) {
        std::cerr << "[WRITE] aligned buffer allocation failed" << std::endl;
        close(fd);
        g_running.store(false);
        return;
    }

    std::cout << "[WRITE] path=" << cfg.write_path
              << ", loop_size=" << static_cast<double>(target_size) / 1024.0 / 1024.0 / 1024.0
              << " GiB"
              << ", rate=" << cfg.write_rate_mib << " MiB/s"
              << ", direct=" << cfg.use_direct
              << std::endl;

    RateLimiter limiter(cfg.write_rate_mib);

    uint64_t offset = 0;

    while (g_running.load()) {
        if (offset + cfg.block_size > target_size) {
            offset = 0;
        }

        ssize_t n = pwrite(fd, buf, cfg.block_size, static_cast<off_t>(offset));
        if (n > 0) {
            cnt.write_bytes.fetch_add(static_cast<uint64_t>(n), std::memory_order_relaxed);
            cnt.write_ops.fetch_add(1, std::memory_order_relaxed);
            offset += static_cast<uint64_t>(n);
            limiter.add_and_wait(static_cast<uint64_t>(n));
        } else {
            cnt.write_errors.fetch_add(1, std::memory_order_relaxed);

            if (errno == EINTR || errno == EAGAIN) {
                continue;
            }

            std::cerr << "[WRITE] pwrite failed, errno=" << errno
                      << ", " << std::strerror(errno) << std::endl;
            break;
        }
    }

    // 测试压力时一般不建议每次都 fsync，否则会引入额外语义。
    // 这里退出时同步一下，避免写入长期停留在设备/系统队列中。
    fsync(fd);

    free(buf);
    close(fd);
}

static bool parse_args(int argc, char** argv, Config& cfg)
{
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];

        auto need_value = [&](const std::string& opt) -> std::string {
            if (i + 1 >= argc) {
                throw std::runtime_error("missing value for " + opt);
            }
            return argv[++i];
        };

        if (a == "--help" || a == "-h") {
            print_help(argv[0]);
            return false;
        } else if (a == "--read-path") {
            cfg.read_path = need_value(a);
        } else if (a == "--write-path") {
            cfg.write_path = need_value(a);
        } else if (a == "--read-rate") {
            cfg.read_rate_mib = std::stod(need_value(a));
        } else if (a == "--write-rate") {
            cfg.write_rate_mib = std::stod(need_value(a));
        } else if (a == "--bs") {
            cfg.block_size = parse_size(need_value(a));
        } else if (a == "--write-size") {
            cfg.write_size = parse_size(need_value(a));
        } else if (a == "--duration") {
            cfg.duration_sec = std::stoi(need_value(a));
        } else if (a == "--direct") {
            cfg.use_direct = (std::stoi(need_value(a)) != 0);
        } else if (a == "--allow-raw-write") {
            cfg.allow_raw_write = true;
        } else {
            throw std::runtime_error("unknown option: " + a);
        }
    }

    if (cfg.read_rate_mib >= 0.0 && cfg.read_path.empty()) {
        throw std::runtime_error("--read-rate is set but --read-path is missing");
    }

    if (cfg.write_rate_mib >= 0.0 && cfg.write_path.empty()) {
        throw std::runtime_error("--write-rate is set but --write-path is missing");
    }

    if (cfg.read_rate_mib < 0.0 && cfg.write_rate_mib < 0.0) {
        throw std::runtime_error("nothing to do: set --read-rate and/or --write-rate");
    }

    if (cfg.use_direct) {
        if (cfg.block_size % 4096 != 0) {
            throw std::runtime_error("with O_DIRECT, --bs should be multiple of 4096");
        }
    }

    return true;
}

int main(int argc, char** argv)
{
    signal(SIGINT, signal_handler);
    signal(SIGTERM, signal_handler);

    Config cfg;

    try {
        if (!parse_args(argc, argv, cfg)) {
            return 0;
        }
    } catch (const std::exception& e) {
        std::cerr << "Argument error: " << e.what() << std::endl;
        print_help(argv[0]);
        return 1;
    }

    std::cout << std::fixed << std::setprecision(2);
    std::cout << "========== M.2 PCIe Load Test ==========" << std::endl;
    std::cout << "block_size=" << static_cast<double>(cfg.block_size) / 1024.0 / 1024.0
              << " MiB"
              << ", duration=" << cfg.duration_sec
              << " sec"
              << ", direct=" << cfg.use_direct
              << std::endl;

    Counters cnt;

    std::thread tr;
    std::thread tw;

    if (cfg.read_rate_mib >= 0.0) {
        tr = std::thread(read_worker, std::cref(cfg), std::ref(cnt));
    }

    if (cfg.write_rate_mib >= 0.0) {
        tw = std::thread(write_worker, std::cref(cfg), std::ref(cnt));
    }

    uint64_t last_read = 0;
    uint64_t last_write = 0;

    auto t0 = std::chrono::steady_clock::now();

    while (g_running.load()) {
        std::this_thread::sleep_for(std::chrono::seconds(1));

        uint64_t r = cnt.read_bytes.load(std::memory_order_relaxed);
        uint64_t w = cnt.write_bytes.load(std::memory_order_relaxed);

        uint64_t dr = r - last_read;
        uint64_t dw = w - last_write;

        last_read = r;
        last_write = w;

        auto now = std::chrono::steady_clock::now();
        double elapsed = std::chrono::duration<double>(now - t0).count();

        double read_mib_s = static_cast<double>(dr) / 1024.0 / 1024.0;
        double write_mib_s = static_cast<double>(dw) / 1024.0 / 1024.0;

        double read_total_gib = static_cast<double>(r) / 1024.0 / 1024.0 / 1024.0;
        double write_total_gib = static_cast<double>(w) / 1024.0 / 1024.0 / 1024.0;

        std::cout << "[T+" << std::setw(6) << elapsed << "s] "
                  << "read=" << std::setw(8) << read_mib_s << " MiB/s, "
                  << "write=" << std::setw(8) << write_mib_s << " MiB/s, "
                  << "read_total=" << std::setw(8) << read_total_gib << " GiB, "
                  << "write_total=" << std::setw(8) << write_total_gib << " GiB, "
                  << "r_err=" << cnt.read_errors.load(std::memory_order_relaxed) << ", "
                  << "w_err=" << cnt.write_errors.load(std::memory_order_relaxed)
                  << std::endl;

        if (cfg.duration_sec > 0 && elapsed >= cfg.duration_sec) {
            g_running.store(false);
            break;
        }
    }

    if (tr.joinable()) tr.join();
    if (tw.joinable()) tw.join();

    std::cout << "========== Finished ==========" << std::endl;

    return 0;
}