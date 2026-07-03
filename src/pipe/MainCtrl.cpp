#include "pipe/MainCtrl.h"
#include "dbs/DbsFusion.hpp"
#include "trackModule.hpp"
#include "pipe/PipeStruDef.h"
#include "../auth/hardware_bind_flow/gate/security_gate.h"
#include "runtime_diagnostics.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cctype>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <mqueue.h>
#include <sstream>
#include <queue>
#include <regex>
#include <time.h>
#include "dbs/lodepng.h"

namespace {

static_assert(sizeof(ResultHeader) == 172U, "ResultHeader must match protocol V1.3 fixed fields");

std::queue<std::string> g_gmtiFileQueue;
std::mutex g_gmtiFileQueueMutex;

std::string composePipePath(const std::string& root, const char* leaf)
{
    if (root.empty()) {
        return std::string(leaf);
    }
    if (root.back() == '/') {
        return root + leaf;
    }
    return root + "/" + leaf;
}

int extractFileIdFromPath(const std::string& path)
{
    const size_t pos = path.find_last_of("/\\");
    const std::string filename = (pos == std::string::npos) ? path : path.substr(pos + 1);
    std::smatch match;
    const std::regex idPattern(R"((\d+)(?:\.(bin|dat))?)", std::regex::icase);
    if (std::regex_search(filename, match, idPattern) && match.size() >= 2) {
        return std::stoi(match[1].str());
    }
    return -1;
}

std::vector<int> buildTrackIdxRange(int id, int window)
{
    std::vector<int> values;
    if (id <= 0 || window <= 0) {
        return values;
    }
    const int start = id - window + 1;
    values.reserve(static_cast<size_t>(window));
    for (int idx = start; idx <= id; ++idx) {
        values.push_back(idx);
    }
    return values;
}

bool allocateNextResultId(GMTIProcessor& proc, const Config& cfg, int& id)
{
    std::string nextPath;
    id = -1;
    return proc.nextGMTIFileName(cfg.result_add, nextPath, id) && id > 0;
}

void clearPendingEchoQueue()
{
    std::queue<std::string> empty;
    std::lock_guard<std::mutex> lock(g_gmtiFileQueueMutex);
    g_gmtiFileQueue.swap(empty);
}

bool popPendingEchoFile(std::string& fileName)
{
    std::lock_guard<std::mutex> lock(g_gmtiFileQueueMutex);
    if (g_gmtiFileQueue.empty()) {
        return false;
    }
    fileName = std::move(g_gmtiFileQueue.front());
    g_gmtiFileQueue.pop();
    return true;
}

void pushPendingEchoFile(const std::string& fileName)
{
    std::lock_guard<std::mutex> lock(g_gmtiFileQueueMutex);
    g_gmtiFileQueue.push(fileName);
}

bool parseLocalEchoSpec(const std::string& spec, int& forcedId, std::string& echoFile)
{
    forcedId = -1;
    echoFile = spec;

    const size_t eq = spec.find('=');
    if (eq == std::string::npos || eq == 0 || eq + 1 >= spec.size()) {
        return true;
    }

    for (size_t i = 0; i < eq; ++i) {
        if (!std::isdigit(static_cast<unsigned char>(spec[i]))) {
            return true;
        }
    }

    forcedId = std::stoi(spec.substr(0, eq));
    echoFile = spec.substr(eq + 1);
    return forcedId > 0 && !echoFile.empty();
}

static inline int32_t quant_deg_to_i32(double deg)
{
    const double q = deg / LSB;
    const double rounded = (q >= 0.0) ? std::floor(q + 0.5) : std::ceil(q - 0.5);
    if (rounded > 2147483647.0) {
        return 2147483647;
    }
    if (rounded < -2147483648.0) {
        return -2147483648;
    }
    return static_cast<int32_t>(rounded);
}

static inline void put_u16_le(std::vector<uint8_t>& buf, size_t off, uint16_t value)
{
    buf[off + 0] = static_cast<uint8_t>(value & 0xFFU);
    buf[off + 1] = static_cast<uint8_t>((value >> 8) & 0xFFU);
}

static inline void put_u32_le(std::vector<uint8_t>& buf, size_t off, uint32_t value)
{
    buf[off + 0] = static_cast<uint8_t>(value & 0xFFU);
    buf[off + 1] = static_cast<uint8_t>((value >> 8) & 0xFFU);
    buf[off + 2] = static_cast<uint8_t>((value >> 16) & 0xFFU);
    buf[off + 3] = static_cast<uint8_t>((value >> 24) & 0xFFU);
}

static inline void put_i32_le(std::vector<uint8_t>& buf, size_t off, int32_t value)
{
    put_u32_le(buf, off, static_cast<uint32_t>(value));
}

static inline void put_double_le(std::vector<uint8_t>& buf, size_t off, double value)
{
    uint64_t raw = 0;
    std::memcpy(&raw, &value, sizeof(raw));
    for (int i = 0; i < 8; ++i) {
        buf[off + static_cast<size_t>(i)] = static_cast<uint8_t>((raw >> (8 * i)) & 0xFFU);
    }
}

static inline uint8_t quant_pixel_pitch(double meters)
{
    const double q = std::round(std::max(0.0, meters) / 0.01);
    if (q > 255.0) {
        return 255;
    }
    return static_cast<uint8_t>(q);
}

static inline void write_target_packet(std::vector<uint8_t>& pkt, size_t off, const GMTIDetection& det)
{
    put_u16_le(pkt, off + 0, det.id);
    put_i32_le(pkt, off + 2, quant_deg_to_i32(det.lon));
    put_i32_le(pkt, off + 6, quant_deg_to_i32(det.lat));
    put_double_le(pkt, off + 10, det.range / 1000.0);
    put_double_le(pkt, off + 18, det.speed);
    put_double_le(pkt, off + 26, det.direction);
    pkt[off + 34] = 0;
}

static std::string makeResultProductPath(const Config& cfg, const char* ext)
{
    if (cfg.result_file_id <= 0) {
        return std::string();
    }
    std::string dir = cfg.result_add;
    if (!dir.empty() && dir.back() != '/' && dir.back() != '\\') {
        dir.push_back('/');
    }
    char name[32];
    std::snprintf(name, sizeof(name), "GMTI%02d.%s", cfg.result_file_id, ext);
    return dir + name;
}

static bool loadDbsImageProduct(const Config& cfg, GMTIResultPacket& packet)
{
    const std::string pngPath = makeResultProductPath(cfg, "png");
    if (pngPath.empty()) {
        return false;
    }

    std::vector<unsigned char> image;
    unsigned width = 0;
    unsigned height = 0;
    const unsigned err = lodepng::decode(image, width, height, pngPath, LCT_GREY, 8);
    if (err != 0U) {
        std::cerr << "[GMTI][WARN] DBS image not available for protocol packet: "
                  << pngPath << " (" << lodepng_error_text(err) << ")" << std::endl;
        return false;
    }
    if (width > 65535U || height > 65535U) {
        std::cerr << "[GMTI][WARN] DBS image exceeds protocol dimensions: "
                  << height << "x" << width << std::endl;
        return false;
    }

    packet.image.assign(image.begin(), image.end());
    packet.image_rows = static_cast<uint16_t>(height);
    packet.image_cols = static_cast<uint16_t>(width);
    packet.image_available = !packet.image.empty();
    return packet.image_available;
}

static bool loadDbsCornerProduct(const Config& cfg, GMTIResultPacket& packet)
{
    const std::string txtPath = makeResultProductPath(cfg, "txt");
    if (txtPath.empty()) {
        return false;
    }

    std::ifstream in(txtPath.c_str());
    if (!in) {
        return false;
    }

    double b[4] = {0.0, 0.0, 0.0, 0.0};
    double l[4] = {0.0, 0.0, 0.0, 0.0};
    bool hasB[4] = {false, false, false, false};
    bool hasL[4] = {false, false, false, false};

    std::string key;
    char eq = '\0';
    double value = 0.0;
    while (in >> key >> eq >> value) {
        if (key.size() == 2 && (key[0] == 'B' || key[0] == 'L') &&
            key[1] >= '0' && key[1] <= '3') {
            const int idx = key[1] - '0';
            if (key[0] == 'B') {
                b[idx] = value;
                hasB[idx] = true;
            } else {
                l[idx] = value;
                hasL[idx] = true;
            }
        }
    }

    for (int i = 0; i < 4; ++i) {
        if (!hasB[i] || !hasL[i]) {
            return false;
        }
    }

    packet.corner_lat[0] = b[3];
    packet.corner_lon[0] = l[3];
    packet.corner_lat[1] = b[0];
    packet.corner_lon[1] = l[0];
    packet.corner_lat[2] = b[2];
    packet.corner_lon[2] = l[2];
    packet.corner_lat[3] = b[1];
    packet.corner_lon[3] = l[1];
    packet.corner_lat[4] = (b[0] + b[1] + b[2] + b[3]) / 4.0;
    packet.corner_lon[4] = (l[0] + l[1] + l[2] + l[3]) / 4.0;
    return true;
}

static GMTIResultPacket buildResultPacketSnapshot(const MainCtrl* host,
                                                  const std::vector<GMTIDetection>& targets)
{
    GMTIResultPacket packet;
    if (!host) {
        return packet;
    }

    packet.targets = targets;
    packet.result_file_id = host->cfg_.result_file_id;
    if (host->cfg_.enable_dbs_fusion) {
        (void)loadDbsImageProduct(host->cfg_, packet);
        (void)loadDbsCornerProduct(host->cfg_, packet);
    }
    return packet;
}

static bool reopen_result_pipe(MainCtrl* host)
{
    if (!host) {
        return false;
    }

    const std::string desiredRoot = host->cfg_.pipe_root_path.empty() ? host->pipe_root_path_ : host->cfg_.pipe_root_path;
    if (desiredRoot == host->pipe_root_path_) {
        return true;
    }

    host->m_SendResultPipe.ClosePipe();
    const std::string resultPipePath = composePipePath(desiredRoot, "pipegmtiresult");
    if (!host->m_SendResultPipe.CreatePipe(resultPipePath.c_str(), O_RDWR)) {
        std::cerr << "Failed to reopen result pipe: " << resultPipePath << std::endl;
        return false;
    }

    host->pipe_root_path_ = desiredRoot;
    return true;
}

static bool runGMTIProcessingFlow(MainCtrl* host, const std::string& echoFile, int forcedResultId = -1)
{
    TIMING_SCOPE(main_total);
    if (!host) {
        return false;
    }

    // security::reset_gate();

    // run_checkpoint(security::CheckpointId::ProgramStartup, "program_startup");
    // run_checkpoint(security::CheckpointId::SarParameterInit, "sar_parameter_init");
    // run_checkpoint(security::CheckpointId::ImagingKernel, "imaging_kernel");
    // run_checkpoint(security::CheckpointId::RangeCmp, "range_cmp");
    // run_checkpoint(security::CheckpointId::MoCo, "moco");
    // run_checkpoint(security::CheckpointId::AzComp, "az_comp");
    // run_checkpoint(security::CheckpointId::PGA, "pga");

    // if (!security::final_decision()) {
    //     std::cout << "[gate] final_decision: failed\n";
    //     return 1;
    // }

    std::string xmlPath;
    {
        std::lock_guard<std::mutex> lock(host->config_mutex_);
        xmlPath = host->gmti_config_xml_;
    }

    const auto config_load_wall_start = std::chrono::system_clock::now();
    const auto config_load_steady_start = std::chrono::high_resolution_clock::now();
    if (!host->gmti_proc_.readXmlParam(xmlPath, host->cfg_)) {
        std::cerr << "Failed to read GMTI XML config" << std::endl;
        return false;
    }

    if (host->pipesEnabled() && !reopen_result_pipe(host)) {
        return false;
    }

    if (!echoFile.empty()) {
        host->cfg_.GMTI_Data_new = echoFile;
        int trackId = forcedResultId > 0 ? forcedResultId : extractFileIdFromPath(echoFile);
        if (trackId <= 0) {
            if (!allocateNextResultId(host->gmti_proc_, host->cfg_, trackId)) {
                std::cerr << "[ERR] Cannot allocate incremental GMTI result id from result_add: "
                          << host->cfg_.result_add << std::endl;
                return false;
            }
            std::cout << "[WARN] Echo filename has no trailing id. Use next result id "
                      << trackId << " for GMTI result and track window." << std::endl;
        }
        host->cfg_.result_file_id = trackId;
        host->cfg_.track_idx_range = buildTrackIdxRange(trackId, host->cfg_.track_idx_window);
        if (forcedResultId > 0) {
            std::cout << "[LOCAL-TEST] Use forced result id " << forcedResultId
                      << " for echo file: " << echoFile << std::endl;
        }
    }

    if (host->cfg_.track_idx_range.empty()) {
        std::cerr << "track_idx_range is empty after GMTI file override" << std::endl;
        return false;
    }

    if (host->cfg_.INFO_Type) {
        host->cfg_.pkg_bytes =
            host->cfg_.info_len + host->cfg_.pulse_len * 16;
    } else if (host->cfg_.channel_mode == "separate") {
        host->cfg_.pkg_bytes = host->cfg_.info_len + host->cfg_.pulse_len * 2 * 4;
    } else if (host->cfg_.channel_mode == "interleaved") {
        host->cfg_.pkg_bytes = host->cfg_.info_len + host->cfg_.pulse_len * 4 * 4;
    } else {
        std::cerr << "[ERR] Unknown channel_mode: " << host->cfg_.channel_mode << std::endl;
        return false;
    }

    host->cfg_.lambda = C / host->cfg_.fc;
    host->cfg_.R_bin = 1.0 / host->cfg_.fs * C / 2.0;
    host->cfg_.Rg.resize(host->cfg_.rg_len);
    for (int i = 0; i < host->cfg_.rg_len; ++i) {
        host->cfg_.Rg[i] = host->cfg_.R_min + i * host->cfg_.R_bin;
    }
    const int procPulseNum = effectivePulseNum(host->cfg_);
    host->cfg_.fd_res = host->cfg_.PRF / double(procPulseNum);
    host->cfg_.az_st = 1;
    host->cfg_.az_ed = procPulseNum;
    host->cfg_.az_center = (host->cfg_.az_st + host->cfg_.az_ed) / 2;
    host->cfg_.Loc = true;

    if (!host->runtime_mode_override_.empty()) {
        host->cfg_.runtime_mode = host->runtime_mode_override_;
        std::transform(host->cfg_.runtime_mode.begin(), host->cfg_.runtime_mode.end(),
                       host->cfg_.runtime_mode.begin(),
                       [](unsigned char ch) { return static_cast<char>(std::tolower(ch)); });
        if (host->cfg_.runtime_mode == "release" || host->cfg_.runtime_mode == "formal" ||
            host->cfg_.runtime_mode == "production") {
            host->cfg_.runtime_diagnostics_enabled = false;
        } else if (host->cfg_.runtime_mode == "debug" || host->cfg_.runtime_mode == "trace") {
            host->cfg_.runtime_diagnostics_enabled = true;
        }
    }
    if (host->runtime_diagnostics_override_set_) {
        host->cfg_.runtime_diagnostics_enabled = host->runtime_diagnostics_override_;
    }
    if (!host->cfg_.runtime_diagnostics_enabled) {
        host->cfg_.track_debug_level = 0;
        host->cfg_.track_debug_dump = false;
        host->cfg_.track_debug_dump_level = 0;
    }
    const auto config_load_wall_end = std::chrono::system_clock::now();
    const auto config_load_steady_end = std::chrono::high_resolution_clock::now();
    gmti::runtime::initializeRun(host->cfg_, xmlPath, "GMTI_pipe_core");
    gmti::runtime::recordTiming(
        "config_load",
        config_load_wall_start,
        config_load_wall_end,
        std::chrono::duration_cast<std::chrono::milliseconds>(
            config_load_steady_end - config_load_steady_start).count(),
        -1,
        "readXmlParam+echoOverride+deriveRuntimeConfig");

    std::vector<std::vector<double>> posMatrix;
    int POS_num = 0;
    if (!host->cfg_.INFO_Type) {
        if (!host->gmti_proc_.POS_dataread(host->cfg_.Plane_POS_add, posMatrix, POS_num)) {
            std::cerr << "[ERR] 读取 POS 失败: " << host->cfg_.Plane_POS_add << std::endl;
            gmti::runtime::finishRun(host->cfg_, false, 1, "POS_dataread failed");
            return false;
        }
    }

    GMTIOutput::Plane plane;
    if (!host->gmti_proc_.extractPlanePV(posMatrix, host->cfg_, plane)) {
        std::cerr << "[ERR] 飞机位置提取失败" << std::endl;
        gmti::runtime::finishRun(host->cfg_, false, 1, "extractPlanePV failed");
        return false;
    }

    std::vector<int> periodList;
    if (!host->gmti_proc_.makePeriodList(plane, host->cfg_, periodList)) {
        std::cerr << "[ERR] 生成 periodList 失败" << std::endl;
        gmti::runtime::finishRun(host->cfg_, false, 1, "makePeriodList failed");
        return false;
    }

    if (host->cfg_.track_idx_range.empty()) {
        host->cfg_.track_idx_range = periodList;
    }

    std::vector<GMTIOutput> periodResults;
    const std::vector<std::vector<double>> emptyPos;
    FusionGroupContext fusionCtx;
    const std::vector<std::vector<double>> &posSource = host->cfg_.INFO_Type ? emptyPos : posMatrix;
    bool processOk = false;
    {
        TIMING_SCOPE(gmti_processing);
        processOk = host->cfg_.enable_dbs_fusion
            ? host->gmti_proc_.processPeriodsParallelFusion(periodList,
                                                            host->cfg_,
                                                            posSource,
                                                            fusionCtx,
                                                            periodResults)
            : host->gmti_proc_.processPeriodsParallel(periodList,
                                                      host->cfg_,
                                                      posSource,
                                                      periodResults);
    }
    if (!processOk) {
        std::cerr << "[ERR] GMTI period processing failed" << std::endl;
    }
    if (host->cfg_.enable_dbs_fusion && processOk) {
        TIMING_SCOPE(dbs_processing);
        if (!runDbsFusionImaging(fusionCtx, host->cfg_, true)) {
            std::cerr << "[WARN] DBS fusion imaging failed; continue with GMTI detections" << std::endl;
        }
    }

    std::vector<double> MT_acc;
    for (size_t i = 0; i < periodResults.size(); ++i) {
        auto& res = periodResults[i];
        if (!res.MT.empty()) {
            MT_acc.insert(MT_acc.end(), std::make_move_iterator(res.MT.begin()),
                          std::make_move_iterator(res.MT.end()));
        }
        // std::cout << "[OK] Period " << periodList[i] << " 完成" << std::endl;
    }

    bool wroteCurrentDetections = false;
    if (!MT_acc.empty()) {
        TIMING_SCOPE(result_write);
        if (!host->gmti_proc_.writeResult(MT_acc, host->cfg_)) {
            std::cerr << "[ERR] 整周期检测结果写盘失败" << std::endl;
        } else {
            wroteCurrentDetections = true;
        }
    } else {
        std::cerr << "[WARN] 本周期无检测目标，跳过检测结果写盘" << std::endl;
    }

    if (!wroteCurrentDetections) {
        std::cout << "[TRACK][WARN] 本周期未写入新的检测结果文件，关联窗口中的当前编号文件可能是旧数据。" << std::endl;
    }
    std::vector<GMTIDetection> currentTargets;
    {
        TIMING_SCOPE(tracking);
        currentTargets = trackModuleOnline(host->cfg_, &host->track_manager_);
    }
    GMTIResultPacket packet = buildResultPacketSnapshot(host, currentTargets);
    {
        std::lock_guard<std::mutex> lock(host->result_mutex_);
        host->latest_gmti_targets_ = currentTargets;
        host->pending_result_packets_.push_back(std::move(packet));
        host->IsResultReady.store(true);
    }
    gmti::runtime::finishRun(host->cfg_, processOk, processOk ? 0 : 1,
                             processOk ? "normal exit" : "period processing reported failure");
    return true;
}

static void printModeSwitchCmd(const ModeSwitchCmd& cmd)
{
    std::cout << std::fixed << std::setprecision(6);

    std::cout << "====== ModeSwitchCmd ======\n";
    std::cout << "head              : 0x" << std::hex << std::uppercase << cmd.head << std::dec << "\n";
    std::cout << "algoType          : 0x" << std::hex << std::uppercase << std::setw(2) << std::setfill('0')
              << static_cast<int>(cmd.algoType) << std::dec << "\n";
    std::cout << "dataLen           : " << cmd.dataLen << "\n";
    std::cout << "workMode          : " << cmd.workMode << "\n";
    std::cout << "CenterFreq        : " << cmd.CenterFreq << " Hz\n";
    std::cout << "SamplingRate      : " << cmd.SamplingRate << " Hz\n";
    std::cout << "PulseWidth        : " << cmd.PulseWidth << " s\n";
    std::cout << "BandWidth         : " << cmd.BandWidth << " Hz\n";
    std::cout << "PRF               : " << cmd.prf << " Hz\n";
    std::cout << "Kr_Sign           : " << cmd.Kr_Sign << "\n";
    std::cout << "SampleDelay       : " << cmd.SampleDelay << " s\n";
    std::cout << "reserved0         : \n";
    std::cout << "LookSide          : " << cmd.LookSide << "\n";
    std::cout << "LookDownAngle     : " << cmd.LookDownAngle << " rad\n";
    std::cout << "SquintAngle       : " << cmd.SquintAngle << " deg\n";
    std::cout << "Theta_bw          : " << cmd.Theta_bw << " rad\n";
    std::cout << "alt_scene         : " << cmd.alt_scene << " m\n";
    std::cout << "Rmin              : " << cmd.Rmin << " m\n";
    std::cout << "rfFreq            : " << cmd.rfFreq << " Hz\n";
    std::cout << "Kr                : " << cmd.Kr << " Hz/s\n";
    std::cout << "ADsamplingLen     : " << cmd.ADsamplingLen << "\n";
    std::cout << "channelSpace      : " << cmd.channelSpace << "\n";
    std::cout << "velocity          : " << cmd.velocity << " m/s\n";
    std::cout << "swath             : " << cmd.swath << " m\n";
    std::cout << "antennaAz         : " << cmd.antennaAz << " deg\n";
    std::cout << "antennaEl         : " << cmd.antennaEl << " deg\n";
    std::cout << "targetLon         : " << cmd.targetLon << " deg\n";
    std::cout << "targetLat         : " << cmd.targetLat << " deg\n";
    std::cout << "targetAlt         : " << cmd.targetAlt << " m\n";
    std::cout << "targetRange       : " << cmd.targetRange << " m\n";
    std::cout << "reserved1         : \n";
    std::cout << "tail              : 0x" << std::hex << cmd.tail << std::dec << "\n";
    std::cout << "===========================\n";
}

bool run_checkpoint(security::CheckpointId id, const char* name) {
    const security::GateCode code = security::checkpoint(id);
    std::cout << "[gate] " << name << ": "
              << security::gate_code_to_string(code) << "\n";
    return code == security::GateCode::Ok;
}

} // namespace

MainCtrl::MainCtrl(const std::string& pipeRootPath, bool enablePipes)
    : pipe_root_path_(pipeRootPath.empty() ? cfg_.pipe_root_path : pipeRootPath),
      pipes_enabled_(enablePipes)
{
    std::cout << "GMTI INIT..." << std::endl;

    m_SendBuf = static_cast<char*>(std::malloc(40 * 1024 * 1024));
    if (!m_SendBuf) {
        std::cerr << "Failed to allocate GMTI send buffer" << std::endl;
        pipes_enabled_ = false;
    }
    IsResultReady.store(false);
    m_workmode.store(Mode_INIT);
    gmti_config_xml_ = "temp_config.xml";

    if (pipes_enabled_) {
        const std::string cmdPipePath = composePipePath(pipe_root_path_, "pipegmticmd");
        const std::string echoPipePath = composePipePath(pipe_root_path_, "pipegmtiecho");
        const std::string resultPipePath = composePipePath(pipe_root_path_, "pipegmtiresult");

        const bool cmdOk = m_RecvCmdPipe.CreatePipe(cmdPipePath.c_str(), O_RDWR);
        const bool echoOk = m_RecvEchoPipe.CreatePipe(echoPipePath.c_str(), O_RDWR);
        const bool resultOk = m_SendResultPipe.CreatePipe(resultPipePath.c_str(), O_RDWR);
        if (!cmdOk || !echoOk || !resultOk) {
            std::cerr << "Failed to create one or more GMTI pipes; threads will not start" << std::endl;
            m_RecvCmdPipe.ClosePipe();
            m_RecvEchoPipe.ClosePipe();
            m_SendResultPipe.ClosePipe();
            pipes_enabled_ = false;
        }
    }

    if (pipes_enabled_) {
        InitThread();
    }
}

MainCtrl::~MainCtrl()
{
    StopThreads();
    if (m_SendBuf) {
        std::free(m_SendBuf);
        m_SendBuf = nullptr;
    }
}

void MainCtrl::InitThread()
{
    if (!pipes_enabled_) {
        return;
    }

    running_.store(true);

    int ret = 0;

    ret = pthread_create(&thrRecvCmd_, nullptr, OnRecvCmdThread, this);
    if (ret != 0) {
        std::cerr << "pthread_create OnRecvCmdThread failed: "
                  << std::strerror(ret) << std::endl;
        running_.store(false);
        StopThreads();
        return;
    }
    recv_cmd_started_ = true;

    ret = pthread_create(&thrRecvEcho_, nullptr, OnEchoRecvThread, this);
    if (ret != 0) {
        std::cerr << "pthread_create OnEchoRecvThread failed: "
                  << std::strerror(ret) << std::endl;
        running_.store(false);
        StopThreads();
        return;
    }
    recv_echo_started_ = true;

    ret = pthread_create(&thrProcData_, nullptr, ProcessDataThread, this);
    if (ret != 0) {
        std::cerr << "pthread_create ProcessDataThread failed: "
                  << std::strerror(ret) << std::endl;
        running_.store(false);
        StopThreads();
        return;
    }
    proc_data_started_ = true;

    ret = pthread_create(&thrSendRes_, nullptr, OnResSendThread, this);
    if (ret != 0) {
        std::cerr << "pthread_create OnResSendThread failed: "
                  << std::strerror(ret) << std::endl;
        running_.store(false);
        StopThreads();
        return;
    }
    send_res_started_ = true;
}

void MainCtrl::StopThreads()
{
    running_.store(false);

    if (pipes_enabled_) {
        m_RecvCmdPipe.ClosePipe();
        m_RecvEchoPipe.ClosePipe();
        m_SendResultPipe.ClosePipe();
    }

    const pthread_t self = pthread_self();
    auto joinThread = [&](pthread_t& thread, bool& started, const char* name) {
        if (!started) {
            thread = {};
            return;
        }
        if (pthread_equal(thread, self)) {
            std::cerr << "[WARN] Skip joining current thread: " << name << std::endl;
            pthread_detach(thread);
            thread = {};
            started = false;
            return;
        }
        const int ret = pthread_join(thread, nullptr);
        if (ret != 0) {
            std::cerr << "pthread_join " << name << " failed: "
                      << std::strerror(ret) << std::endl;
        }
        thread = {};
        started = false;
    };

    joinThread(thrRecvCmd_, recv_cmd_started_, "OnRecvCmdThread");
    joinThread(thrRecvEcho_, recv_echo_started_, "OnEchoRecvThread");
    joinThread(thrProcData_, proc_data_started_, "ProcessDataThread");
    joinThread(thrSendRes_, send_res_started_, "OnResSendThread");
}

bool MainCtrl::RunLocalTest(const std::string& xmlPath, const std::vector<std::string>& echoFiles)
{
    if (echoFiles.empty()) {
        std::cerr << "[LOCAL-TEST][ERR] echo file list is empty" << std::endl;
        return false;
    }

    const std::string localXmlPath = xmlPath.empty() ? "temp_config.xml" : xmlPath;
    {
        std::lock_guard<std::mutex> lock(config_mutex_);
        gmti_config_xml_ = localXmlPath;
    }
    m_workmode.store(Mode_GMTI);
    track_manager_.reset();

    bool allOk = true;
    std::cout << "[LOCAL-TEST] XML: " << localXmlPath << std::endl;
    std::cout << "[LOCAL-TEST] Echo files: " << echoFiles.size() << std::endl;
    for (size_t i = 0; i < echoFiles.size(); ++i) {
        int forcedId = -1;
        std::string echoFile;
        if (!parseLocalEchoSpec(echoFiles[i], forcedId, echoFile)) {
            allOk = false;
            std::cerr << "[LOCAL-TEST][ERR] Invalid echo spec: " << echoFiles[i] << std::endl;
            continue;
        }

        std::cout << "\n[LOCAL-TEST] ===== "
                  << (i + 1) << "/" << echoFiles.size()
                  << " echo=" << echoFile;
        if (forcedId > 0) {
            std::cout << " result_id=" << forcedId;
        }
        std::cout << " =====" << std::endl;
        if (!runGMTIProcessingFlow(this, echoFile, forcedId)) {
            allOk = false;
            std::cerr << "[LOCAL-TEST][ERR] GMTI processing failed: "
                      << echoFile << std::endl;
        }
        {
            std::lock_guard<std::mutex> lock(result_mutex_);
            pending_result_packets_.clear();
            IsResultReady.store(false);
        }
    }

    return allOk;
}

void* MainCtrl::OnRecvCmdThread(void *param)
{
    MainCtrl *pHost = static_cast<MainCtrl*>(param);
    ModeSwitchCmd cmd;

    std::cout << "sizeof ModeSwitchCmd = " << sizeof(ModeSwitchCmd) << std::endl;

    while(pHost->running_.load()) {
        const int ret = pHost->m_RecvCmdPipe.ReadData(reinterpret_cast<char*>(&cmd), sizeof(ModeSwitchCmd));
        std::cout << "Recv cmd ret = " << ret << std::endl;

        if(ret < 0) {
            if (!pHost->running_.load()) {
                break;
            }
            printf("Receive mode switch command error\n");
            usleep(100000);
            continue;
        }
        if(ret != static_cast<int>(sizeof(ModeSwitchCmd))) {
            if (!pHost->running_.load()) {
                break;
            }
            printf("Receive incomplete mode switch command\n");
            usleep(100000);
            continue;
        }

        printModeSwitchCmd(cmd);
        std::string xmlPath;
        {
            std::lock_guard<std::mutex> lock(pHost->config_mutex_);
            pHost->last_cmd_ = cmd;
            if(cmd.workMode == Mode_GMTI)
            {
                pHost->gmti_config_xml_ = "temp_config.xml";
            }
            xmlPath = pHost->gmti_config_xml_;
        }

        if(pHost->m_workmode.load() == cmd.workMode)
        {
            printf("Already in cmd workmode\n");
        }
        else if(cmd.workMode == Mode_INIT)
        {
            printf("Switch to MODE_INIT\n");
        }
        else if(cmd.workMode == Mode_GMTI)
        {
            printf("Switch to Mode_GMTI (mode=%d)\n", cmd.workMode);
        }
        else
        {
            printf("Switch to mode %d\n", cmd.workMode);
        }

        if (UPDATEXML) {
            if (!pHost->updateXmlFromModeSwitchCmd(cmd, xmlPath)) {
                std::cerr << "Failed to update XML config from ModeSwitchCmd: "
                          << xmlPath << std::endl;
            }
        }

        pHost->m_workmode.store(cmd.workMode);
    }
    return nullptr;
}

void* MainCtrl::OnEchoRecvThread(void *param)
{
    MainCtrl *pHost = static_cast<MainCtrl*>(param);

    std::cout << "OnEchoRecvThread start" << std::endl;

    char buf[512] = {0};
    while(pHost->running_.load())
    {
        const int ret = pHost->m_RecvEchoPipe.ReadData(buf, 512);
        if(ret < 1) {
            if (!pHost->running_.load()) {
                break;
            }
            printf("Receive echo data error\n");
            usleep(100000);
            continue;
        }

        std::string fileName(buf, static_cast<size_t>(ret));
        while (!fileName.empty() && (fileName.back() == '\0' || fileName.back() == '\n' || fileName.back() == '\r')) {
            fileName.pop_back();
        }
        if (!fileName.empty()) {
            pushPendingEchoFile(fileName);
        }
    }
    return nullptr;
}

void* MainCtrl::ProcessDataThread(void* param)
{
    MainCtrl *pHost = static_cast<MainCtrl*>(param);

    std::cout << "ProcessDataThread start" << std::endl;

    while(pHost->running_.load()) {
        std::string fileName;
        if (!popPendingEchoFile(fileName))
        {
            usleep(10000);
            continue;
        }

        std::cout << "open file: " << fileName << std::endl;

        const bool ok = runGMTIProcessingFlow(pHost, fileName);
        {
            std::lock_guard<std::mutex> lock(pHost->result_mutex_);
            pHost->IsResultReady.store(!pHost->pending_result_packets_.empty());
        }
        if (!ok) {
            std::cerr << "[ERR] GMTI processing failed for echo file: " << fileName << std::endl;
        }

        if (pHost->m_workmode.load() != Mode_GMTI)
        {
            clearPendingEchoQueue();
            std::cout << "Standby mode: cleared pending echo queue after current file" << std::endl;
        }
    }
    return nullptr;
}

void* MainCtrl::OnResSendThread(void* param)
{
    MainCtrl *pHost = static_cast<MainCtrl*>(param);

    std::cout << "OnResSendThread start" << std::endl;

    while(pHost->running_.load())
    {
        GMTIResultPacket packet;
        bool hasPacket = false;
        {
            std::lock_guard<std::mutex> lock(pHost->result_mutex_);
            if (pHost->pending_result_packets_.empty()) {
                pHost->IsResultReady.store(false);
            } else {
                packet = pHost->pending_result_packets_.front();
                pHost->IsResultReady.store(true);
                hasPacket = true;
            }
        }

        if (!hasPacket)
        {
            usleep(100000);
            continue;
        }

        std::cout << "GMTI mode: Sending protocol packet..." << std::endl;

        uint32_t packed_len = 0;
        if (!pHost->packGMTIResults(packet, pHost->m_SendBuf, 40U * 1024U * 1024U, packed_len))
        {
            std::cerr << "Failed to pack GMTI protocol packet" << std::endl;
            usleep(100000);
            continue;
        }

        const ssize_t written = pHost->m_SendResultPipe.WriteData(pHost->m_SendBuf, static_cast<int>(packed_len));
        if (!pHost->running_.load()) {
            break;
        }
        if (written == static_cast<ssize_t>(packed_len))
        {
            std::cout << "GMTI protocol packet sent: " << packed_len << " bytes" << std::endl;
            std::lock_guard<std::mutex> lock(pHost->result_mutex_);
            if (!pHost->pending_result_packets_.empty()) {
                pHost->pending_result_packets_.pop_front();
            }
            pHost->latest_gmti_targets_.clear();
            pHost->latest_gmti_targets_.shrink_to_fit();
            pHost->IsResultReady.store(!pHost->pending_result_packets_.empty());
        }
        else
        {
            std::cerr << "Failed to send GMTI protocol packet" << std::endl;
            usleep(100000);
        }
        std::cout << "end send res" << std::endl;
    }
    return nullptr;
}

bool MainCtrl::updateXmlFromModeSwitchCmd(const ModeSwitchCmd& cmd, const std::string& xmlPath)
{
    TiXmlDocument doc(xmlPath.c_str());
    if (!doc.LoadFile())
    {
        std::cerr << "Failed to load XML file: " << xmlPath << std::endl;
        return false;
    }

    TiXmlElement* root = doc.FirstChildElement("GMTI");
    if (!root)
    {
        std::cerr << "Invalid XML structure: root <GMTI> not found" << std::endl;
        return false;
    }

    TiXmlElement* param = root->FirstChildElement("GMTI_parameter");
    if (!param)
    {
        std::cerr << "Invalid XML structure: <GMTI_parameter> not found" << std::endl;
        return false;
    }

    auto setText = [&](const char* name, double value) {
        TiXmlElement* elem = param->FirstChildElement(name);
        if (!elem) {
            elem = new TiXmlElement(name);
            param->LinkEndChild(elem);
        }
        std::ostringstream oss;
        oss << std::setprecision(15) << value;
        elem->Clear();
        elem->LinkEndChild(new TiXmlText(oss.str().c_str()));
    };

    auto setTextInt = [&](const char* name, int value) {
        TiXmlElement* elem = param->FirstChildElement(name);
        if (!elem) {
            elem = new TiXmlElement(name);
            param->LinkEndChild(elem);
        }
        std::ostringstream oss;
        oss << value;
        elem->Clear();
        elem->LinkEndChild(new TiXmlText(oss.str().c_str()));
    };

    if (cmd.CenterFreq > 0)
    {
        setText("fc", cmd.CenterFreq / 1e9);
    }

    if (cmd.SamplingRate > 0)
    {
        setText("fs", cmd.SamplingRate / 1e6);
    }

    if (cmd.prf > 0)
    {
        setText("PRF", cmd.prf);
    }

    if (cmd.BandWidth > 0)
    {
        setText("Br", cmd.BandWidth / 1e6);
    }

    if (cmd.PulseWidth > 0)
    {
        setText("Tr", cmd.PulseWidth * 1e6);
    }

    if (cmd.alt_scene > 0)
    {
        setText("alt_scene", cmd.alt_scene);
    }

    if (cmd.Rmin > 0)
    {
        setText("Rmin", cmd.Rmin);
    }

    if (cmd.velocity > 0)
    {
        setText("v_platform", cmd.velocity);
    }

    if (cmd.LookSide == 0 || cmd.LookSide == 1)
    {
        setTextInt("squint_side", static_cast<int>(cmd.LookSide));
    }

    if (cmd.SquintAngle != 0)
    {
        setText("squint_angle", cmd.SquintAngle);
    }

    if (cmd.targetLon != 0)
    {
        setText("targetLon", cmd.targetLon);
    }

    if (cmd.targetLat != 0)
    {
        setText("targetLat", cmd.targetLat);
    }

    if (!doc.SaveFile())
    {
        std::cerr << "Failed to save XML file: " << xmlPath << std::endl;
        return false;
    }

    std::cout << "GMTI parameters updated from ModeSwitchCmd" << std::endl;
    return true;
}

bool MainCtrl::validateGMTIParams()
{
    if (cfg_.pulse_num <= 0 || cfg_.rg_len <= 0)
    {
        std::cerr << "Invalid pulse_num or rg_len" << std::endl;
        return false;
    }

    if (cfg_.fc <= 0 || cfg_.fs <= 0 || cfg_.PRF <= 0)
    {
        std::cerr << "Invalid frequency parameters" << std::endl;
        return false;
    }

    if (cfg_.GMTI_Data_new.empty())
    {
        std::cerr << "GMTI data path not specified" << std::endl;
        return false;
    }

    return true;
}

bool MainCtrl::packGMTIResults(const GMTIResultPacket& packet, char* buffer, size_t buffer_size, uint32_t& packed_len)
{
    constexpr size_t kTargetPacketSize = 35U;

    if (!buffer || buffer_size < sizeof(ResultHeader))
    {
        std::cerr << "Invalid buffer for packing results" << std::endl;
        return false;
    }

    const size_t target_count = std::min<size_t>(packet.targets.size(), 65535U);
    const size_t image_bytes = packet.image_available ? packet.image.size() : 0U;
    const size_t total_bytes = sizeof(ResultHeader) + target_count * kTargetPacketSize + image_bytes;
    if (total_bytes > buffer_size) {
        std::cerr << "Buffer overflow: GMTI packet too large" << std::endl;
        return false;
    }
    if (total_bytes > 0xFFFFFFFFULL) {
        std::cerr << "GMTI packet length exceeds protocol U32 range" << std::endl;
        return false;
    }

    std::vector<uint8_t> pkt(total_bytes, 0U);
    ResultHeader header{};
    header.head = 0xAA55;
    header.msgLen = static_cast<uint32_t>(total_bytes - sizeof(uint16_t) - sizeof(uint8_t));
    header.msgAddr = 0;
    header.msgType = 0;
    header.msgCount = next_result_msg_count_++;
    header.srcId = 0;
    header.dstId = 0;
    header.cmdType = static_cast<uint8_t>(Mode_GMTI);
    header.cmdCount = static_cast<uint8_t>(header.msgCount & 0xFFU);
    header.height = packet.image_rows;
    header.width = packet.image_cols;
    header.availFlag = packet.image_available ? 0xFFFF : 0x0000;
    header.roll = 0;
    header.heading = 0;
    header.pitch = 0;
    header.navLon = 0;
    header.navLat = 0;
    header.navHeight = 0;
    header.velNorth = 0;
    header.velUp = 0;
    header.velEast = 0;
    header.hour = 0;
    header.minute = 0;
    header.second = 0;
    header.millisec = 0;
    std::memset(header.posReserve, 0, sizeof(header.posReserve));
    header.hLeftTop = 0;
    header.hLeftDown = 0;
    header.hRightDown = 0;
    header.hRightTop = 0;
    header.hCenter = 0;
    header.lonLeftTop = quant_deg_to_i32(packet.corner_lon[0]);
    header.lonLeftDown = quant_deg_to_i32(packet.corner_lon[1]);
    header.lonRightDown = quant_deg_to_i32(packet.corner_lon[2]);
    header.lonRightTop = quant_deg_to_i32(packet.corner_lon[3]);
    header.lonCenter = quant_deg_to_i32(packet.corner_lon[4]);
    header.latLeftTop = quant_deg_to_i32(packet.corner_lat[0]);
    header.latLeftDown = quant_deg_to_i32(packet.corner_lat[1]);
    header.latRightDown = quant_deg_to_i32(packet.corner_lat[2]);
    header.latRightTop = quant_deg_to_i32(packet.corner_lat[3]);
    header.latCenter = quant_deg_to_i32(packet.corner_lat[4]);
    header.rLeftTop = 0;
    header.rLeftDown = 0;
    header.rRightDown = 0;
    header.rRightTop = 0;
    header.rCenter = 0;
    header.reserve1 = 0;
    header.pixelPitch = quant_pixel_pitch(cfg_.dbs_out_res_m);
    header.LookDownAngle = 0;
    header.SquintAngle = static_cast<uint16_t>(std::lround(std::fabs(cfg_.squint_angle) * 100.0));
    header.LookSide = (cfg_.squint_side == 0) ? 0x00 : 0xFF;
    std::memset(header.reserve2, 0, sizeof(header.reserve2));
    header.checksum = 0;
    header.targetNum = static_cast<uint16_t>(target_count);

    std::memcpy(pkt.data(), &header, sizeof(header));

    size_t off = sizeof(ResultHeader);
    for (size_t i = 0; i < target_count; ++i) {
        write_target_packet(pkt, off, packet.targets[i]);
        off += kTargetPacketSize;
    }
    if (image_bytes > 0U) {
        std::memcpy(pkt.data() + off, packet.image.data(), image_bytes);
        off += image_bytes;
    }

    uint8_t checksum = 0;
    for (size_t i = 2; i < 169; ++i) {
        checksum = static_cast<uint8_t>(checksum + pkt[i]);
    }
    pkt[169] = checksum;

    std::memcpy(buffer, pkt.data(), total_bytes);
    packed_len = static_cast<uint32_t>(total_bytes);
    return true;
}
