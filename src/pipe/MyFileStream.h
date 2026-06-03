#ifndef MYFILESTREAM_H
#define MYFILESTREAM_H

#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>


// FileOnlyStream 类定义
class FileOnlyStream {
private:
    std::ofstream file_stream;  // 文件流对象
    bool is_initialized = false;

public:

    FileOnlyStream() = default;

    void open(const std::string& filename,
              const std::string& mode = "append") {
        if (mode == "append") {
            file_stream.open(filename, std::ios::app);
        } else {
            file_stream.open(filename);
        }

        if (!file_stream.is_open()) {
            std::cerr << ":  '" << filename << "'" << std::endl;
        } else {
            is_initialized = true;
        }
    }


    /******************************
     * 重载 << 操作符
     * 支持所有可以通过 << 输出的数据类型
     ******************************/
    template<typename T>
    FileOnlyStream& operator<<(const T& value) {
        if (file_stream.is_open()) {
            file_stream << value;  // 输出到文件
        }
        return *this;  // 支持链式调用
    }

    /******************************
     * 重载 << 操作符（处理endl等操作符）
     ******************************/
    FileOnlyStream& operator<<(std::ostream& (*manip)(std::ostream&)) {
        if (file_stream.is_open()) {
            file_stream << manip;  // 处理std::endl, std::flush等
        }
        return *this;
    }

    /******************************
     * 检查文件是否成功打开
     ******************************/
    bool is_open() const {
        return file_stream.is_open();
    }

    /******************************
     * 刷新缓冲区
     * 立即将缓冲区内容写入磁盘
     ******************************/
    void flush() {
        if (file_stream.is_open()) {
            file_stream.flush();
        }
    }

    /******************************
     * 关闭文件
     ******************************/
    void close() {
        if (file_stream.is_open()) {
            file_stream.close();
        }
    }

    /******************************
     * 析构函数
     * 自动关闭文件
     ******************************/
    ~FileOnlyStream() {
        close();
    }
};

#endif // MYFILESTREAM_H
