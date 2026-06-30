#pragma once

#include <atomic>
#include <cstring>
#include <semaphore.h>
#include <errno.h>
#include <stdexcept>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <string>

// ---------- ----------
struct SharedControl {
    sem_t sem_empty;
    sem_t sem_full;
    std::atomic<uint64_t> write_pos;   // [0, capacity) 
    std::atomic<uint64_t> read_pos;    // [0, capacity) 
    uint64_t capacity;                 // 

    static void init(SharedControl* ctrl, size_t cap) {
        if ((cap & (cap - 1)) != 0)
            throw std::invalid_argument("capacity must be power of two");
        ctrl->capacity = cap;
        ctrl->write_pos.store(0, std::memory_order_relaxed);
        ctrl->read_pos.store(0, std::memory_order_relaxed);
        if (sem_init(&ctrl->sem_empty, 1, cap) != 0 ||
            sem_init(&ctrl->sem_full, 1, 0) != 0) {
            throw std::runtime_error("sem_init failed");
        }
    }

    static void destroy(SharedControl* ctrl) {
        sem_destroy(&ctrl->sem_empty);
        sem_destroy(&ctrl->sem_full);
    }
};

// ----------  ----------
class RingBuffer {
public:
    explicit RingBuffer(void* base_ptr)
        : ctrl_(static_cast<SharedControl*>(base_ptr)),
          data_ptr_(static_cast<char*>(base_ptr) + sizeof(SharedControl)) {}

    // 
    bool write(const void* data, size_t len) {
        if (len == 0) return true;
        if (len > ctrl_->capacity) return false;

        while (sem_wait(&ctrl_->sem_empty) != 0) {
            if (errno != EINTR) return false;
        }

        uint64_t current_write = ctrl_->write_pos.load(std::memory_order_relaxed);
        size_t offset = current_write & (ctrl_->capacity - 1);
        size_t first_chunk = std::min(len, ctrl_->capacity - offset);
        memcpy(data_ptr_ + offset, data, first_chunk);
        if (len > first_chunk) {
            memcpy(data_ptr_, static_cast<const char*>(data) + first_chunk, len - first_chunk);
        }

        // 
        uint64_t new_write = (current_write + len) & (ctrl_->capacity - 1);
        ctrl_->write_pos.store(new_write, std::memory_order_release);
        sem_post(&ctrl_->sem_full);
        return true;
    }

    // false
    bool try_write(const void* data, size_t len) {
        if (len == 0) return true;
        if (len > ctrl_->capacity) return false;
        if (sem_trywait(&ctrl_->sem_empty) != 0) return false;

        uint64_t current_write = ctrl_->write_pos.load(std::memory_order_relaxed);
        size_t offset = current_write & (ctrl_->capacity - 1);
        size_t first_chunk = std::min(len, ctrl_->capacity - offset);
        memcpy(data_ptr_ + offset, data, first_chunk);
        if (len > first_chunk) {
            memcpy(data_ptr_, static_cast<const char*>(data) + first_chunk, len - first_chunk);
        }

        uint64_t new_write = (current_write + len) & (ctrl_->capacity - 1);
        ctrl_->write_pos.store(new_write, std::memory_order_release);
        sem_post(&ctrl_->sem_full);
        return true;
    }

    // 
    bool read(void* out_buf, size_t& io_len) {
        if (out_buf == nullptr || io_len == 0) return false;

        while (sem_wait(&ctrl_->sem_full) != 0) {
            if (errno != EINTR) return false;
        }

        uint64_t current_write = ctrl_->write_pos.load(std::memory_order_acquire);
        uint64_t current_read  = ctrl_->read_pos.load(std::memory_order_relaxed);

        // capacity
        size_t available = (current_write - current_read + ctrl_->capacity) & (ctrl_->capacity - 1);
        size_t len = std::min(io_len, available);

        size_t offset = current_read & (ctrl_->capacity - 1);
        size_t first_chunk = std::min(len, ctrl_->capacity - offset);
        memcpy(out_buf, data_ptr_ + offset, first_chunk);
        if (len > first_chunk) {
            memcpy(static_cast<char*>(out_buf) + first_chunk, data_ptr_, len - first_chunk);
        }

        // 
        uint64_t new_read = (current_read + len) & (ctrl_->capacity - 1);
        ctrl_->read_pos.store(new_read, std::memory_order_release);
        sem_post(&ctrl_->sem_empty);

        io_len = len;
        return true;
    }

    // falseio_len 0
    bool try_read(void* out_buf, size_t& io_len) {
        if (out_buf == nullptr || io_len == 0) return false;
        if (sem_trywait(&ctrl_->sem_full) != 0) {
            io_len = 0;
            return false;
        }

        uint64_t current_write = ctrl_->write_pos.load(std::memory_order_acquire);
        uint64_t current_read  = ctrl_->read_pos.load(std::memory_order_relaxed);

        size_t available = (current_write - current_read + ctrl_->capacity) & (ctrl_->capacity - 1);
        size_t len = std::min(io_len, available);

        size_t offset = current_read & (ctrl_->capacity - 1);
        size_t first_chunk = std::min(len, ctrl_->capacity - offset);
        memcpy(out_buf, data_ptr_ + offset, first_chunk);
        if (len > first_chunk) {
            memcpy(static_cast<char*>(out_buf) + first_chunk, data_ptr_, len - first_chunk);
        }

        uint64_t new_read = (current_read + len) & (ctrl_->capacity - 1);
        ctrl_->read_pos.store(new_read, std::memory_order_release);
        sem_post(&ctrl_->sem_empty);

        io_len = len;
        return true;
    }

    size_t capacity() const { return ctrl_->capacity; }

private:
    SharedControl* ctrl_;
    char*          data_ptr_;
};

// ---------- / +  ----------
class RingBufferManager {
public:
    enum Role { SENDER, RECEIVER };

    //  buffer_size
    // buffer_size 
    RingBufferManager(const std::string& name, Role role, size_t buffer_size = 0)
        : name_(name), role_(role), base_(nullptr), total_size_(0), fd_(-1), ring_(nullptr) {
        if (role == SENDER) {
            if (buffer_size == 0 || (buffer_size & (buffer_size - 1)) != 0)
                throw std::invalid_argument("buffer_size must be power of two");
            create_sender(buffer_size);
        } else {
            attach_receiver();
        }
    }

    ~RingBufferManager() {
        cleanup();
    }

    // 
    RingBufferManager(const RingBufferManager&) = delete;
    RingBufferManager& operator=(const RingBufferManager&) = delete;
    RingBufferManager(RingBufferManager&& other) noexcept
        : name_(std::move(other.name_)), role_(other.role_),
          base_(other.base_), total_size_(other.total_size_),
          fd_(other.fd_), ring_(other.ring_) {
        other.base_ = nullptr;
        other.fd_ = -1;
        other.ring_ = nullptr;
    }
    RingBufferManager& operator=(RingBufferManager&& other) noexcept {
        if (this != &other) {
            cleanup();
            name_ = std::move(other.name_);
            role_ = other.role_;
            base_ = other.base_;
            total_size_ = other.total_size_;
            fd_ = other.fd_;
            ring_ = other.ring_;
            other.base_ = nullptr;
            other.fd_ = -1;
            other.ring_ = nullptr;
        }
        return *this;
    }

    RingBuffer& get_ring() {
        if (!ring_) throw std::runtime_error("RingBuffer not initialized");
        return *ring_;
    }

private:
    void create_sender(size_t buffer_size) {
        shm_unlink(name_.c_str());

        fd_ = shm_open(name_.c_str(), O_CREAT | O_RDWR | O_EXCL, 0666);
        if (fd_ == -1)
            throw std::runtime_error("shm_open failed");

        total_size_ = sizeof(SharedControl) + buffer_size;
        if (ftruncate(fd_, total_size_) == -1)
            throw std::runtime_error("ftruncate failed");

        base_ = mmap(nullptr, total_size_, PROT_READ | PROT_WRITE, MAP_SHARED, fd_, 0);
        if (base_ == MAP_FAILED)
            throw std::runtime_error("mmap failed");

        SharedControl::init(static_cast<SharedControl*>(base_), buffer_size);

        ring_ = new (std::nothrow) RingBuffer(base_);
        if (!ring_)
            throw std::runtime_error("new RingBuffer failed");
    }

    void attach_receiver() {
        fd_ = shm_open(name_.c_str(), O_RDWR, 0666);
        if (fd_ == -1)
            throw std::runtime_error("shm_open for receiver failed");

        struct stat st;
        if (fstat(fd_, &st) == -1)
            throw std::runtime_error("fstat failed");
        total_size_ = st.st_size;

        base_ = mmap(nullptr, total_size_, PROT_READ | PROT_WRITE, MAP_SHARED, fd_, 0);
        if (base_ == MAP_FAILED)
            throw std::runtime_error("mmap failed");

        ring_ = new (std::nothrow) RingBuffer(base_);
        if (!ring_)
            throw std::runtime_error("new RingBuffer failed");
    }

    void cleanup() {
        if (ring_) {
            delete ring_;
            ring_ = nullptr;
        }
        if (base_ && base_ != MAP_FAILED) {
            if (role_ == SENDER) {
                SharedControl::destroy(static_cast<SharedControl*>(base_));
            }
            munmap(base_, total_size_);
            base_ = nullptr;
        }
        if (fd_ != -1) {
            close(fd_);
            fd_ = -1;
        }
        if (role_ == SENDER) {
            shm_unlink(name_.c_str());
        }
    }

private:
    std::string name_;
    Role        role_;
    void*       base_;
    size_t      total_size_;
    int         fd_;
    RingBuffer* ring_;
};
