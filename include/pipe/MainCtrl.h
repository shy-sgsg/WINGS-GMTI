#ifndef MAINCTRL_H
#define MAINCTRL_H

#include "PipeRW.h"
#include "PipeStruDef.h"
#include <queue>
#include <mqueue.h>
#include <stdio.h>
#include "MyFileStream.h"
#include "config_structs.hpp"
#include "trackModule.hpp"
#include "TrackManager.hpp"
#include "GMTIProcessor.hpp"
#include <vector>
#include <mutex>
#include <deque>
#include <cstdint>
#include <atomic>

#define UPDATEXML 0

struct GMTIResultPacket {
    std::vector<GMTIDetection> targets;
    std::vector<uint8_t> image;
    uint16_t image_rows = 0;
    uint16_t image_cols = 0;
    bool image_available = false;
    double corner_lon[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
    double corner_lat[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
    int result_file_id = 0;
};

class MainCtrl {

public:
    explicit MainCtrl(const std::string& pipeRootPath = "/home/shy/pipe_test",
                      bool enablePipes = true);
    ~MainCtrl();

    std::atomic<bool> running_{false};

    pthread_t thrRecvCmd_{};
    pthread_t thrRecvEcho_{};
    pthread_t thrProcData_{};
    pthread_t thrSendRes_{};

    // 命令接收线程
    static void* OnRecvCmdThread(void* param);
    // 回波数据接收线程
    static void* OnEchoRecvThread(void* param);
    // 回波数据处理线程
    static void* ProcessDataThread(void* param);
    // 结果发送线程
    static void* OnResSendThread(void* param);

    std::atomic<bool> IsResultReady{false};

    std::atomic<uint32_t> m_workmode{Mode_INIT};
    
    char * m_SendBuf;

    PipeRW m_RecvCmdPipe;
    PipeRW m_RecvEchoPipe;
    PipeRW m_SendResultPipe;

    void InitThread();
    void StopThreads();

    FileOnlyStream file_stream;

    void initFileStream(const std::string& filename,
                       const std::string& mode = "append") {
        file_stream.open(filename, mode);
    }

    // ===== 新增成员：GMTI 专用 =====
    Config cfg_;                           // GMTI 配置对象
    GMTIProcessor gmti_proc_;              // GMTI 处理器实例
    std::vector<std::vector<double>> all_frames_MT_;  // 所有帧的检测目标
    std::vector<int> periodList_;          // 周期列表
    std::vector<GMTIOutput> periodResults_;// 周期处理结果
    std::vector<Track> current_tracks_;    // 当前航迹结果
    TrackManager track_manager_;           // 跨周期持久航迹编号管理器
    std::vector<GMTIDetection> latest_gmti_targets_; // 最新 GMTI 目标
    std::deque<GMTIResultPacket> pending_result_packets_; // 待发送结果快照
    std::string latest_result_file_;       // 最新生成的 GMTI 航迹结果文件
    
    std::string gmti_config_xml_;          // GMTI 配置 XML 文件路径
    std::string runtime_mode_override_;    // 命令行覆盖 runtime_mode，空表示使用 XML
    bool runtime_diagnostics_override_set_ = false;
    bool runtime_diagnostics_override_ = true;
    ModeSwitchCmd last_cmd_;               // 最后接收的模式切换命令
    std::mutex result_mutex_;              // 保护 latest_gmti_targets_
    std::mutex config_mutex_;              // 保护 gmti_config_xml_ / last_cmd_
    std::string pipe_root_path_;           // 管道主路径

    // ===== 新增函数：GMTI 参数和结果处理 =====
    // 从 ModeSwitchCmd 更新 XML 文件中的参数
    bool updateXmlFromModeSwitchCmd(const ModeSwitchCmd& cmd, const std::string& xmlPath);

    // 验证 GMTI 参数的有效性
    bool validateGMTIParams();

    // 打包最新 GMTI 结果为主控协议格式
    bool packGMTIResults(const GMTIResultPacket& packet, char* buffer, size_t buffer_size, uint32_t& packed_len);

    // 本地测试模式：不依赖上位机管道，按给定回波文件顺序复用同一个 TrackManager。
    bool RunLocalTest(const std::string& xmlPath, const std::vector<std::string>& echoFiles);
    bool pipesEnabled() const { return pipes_enabled_; }

private:
    uint16_t next_result_msg_count_ = 0;
    bool pipes_enabled_ = true;
    bool recv_cmd_started_ = false;
    bool recv_echo_started_ = false;
    bool proc_data_started_ = false;
    bool send_res_started_ = false;

};

#endif // MAINCTRL_H
