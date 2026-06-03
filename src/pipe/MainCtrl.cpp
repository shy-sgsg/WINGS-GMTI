#include "MainCtrl.h"
#include "trackModule.hpp"
#include "PipeStruDef.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <mqueue.h>
#include <queue>
#include <regex>
#include <time.h>

namespace {

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
    const std::regex idPattern(R"((\d+)(?:\.bin)?$)", std::regex::icase);
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

static inline uint16_t quant_speed_to_u16(double speed)
{
    const double q = std::round(std::max(0.0, speed) / 0.01);
    if (q > 65535.0) {
        return 65535;
    }
    return static_cast<uint16_t>(q);
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

static inline void write_target_packet(std::vector<uint8_t>& pkt, size_t off, const GMTIDetection& det)
{
    put_u16_le(pkt, off + 0, det.id);
    put_i32_le(pkt, off + 2, quant_deg_to_i32(det.lon));
    put_i32_le(pkt, off + 6, quant_deg_to_i32(det.lat));
    put_u16_le(pkt, off + 10, quant_speed_to_u16(det.speed));

    uint64_t dir_raw = 0;
    uint64_t range_raw = 0;
    std::memcpy(&dir_raw, &det.direction, sizeof(dir_raw));
    std::memcpy(&range_raw, &det.range, sizeof(range_raw));
    for (int i = 0; i < 8; ++i) {
        pkt[off + 12 + static_cast<size_t>(i)] = static_cast<uint8_t>((dir_raw >> (8 * i)) & 0xFFU);
        pkt[off + 20 + static_cast<size_t>(i)] = static_cast<uint8_t>((range_raw >> (8 * i)) & 0xFFU);
    }
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

static bool runGMTIProcessingFlow(MainCtrl* host, const std::string& echoFile)
{
    TIMING_SCOPE(main_total);
    if (!host) {
        return false;
    }

    if (!host->gmti_proc_.readXmlParam(host->gmti_config_xml_, host->cfg_)) {
        std::cerr << "Failed to read GMTI XML config" << std::endl;
        return false;
    }

    if (!reopen_result_pipe(host)) {
        return false;
    }

    if (!echoFile.empty()) {
        host->cfg_.GMTI_Data_new = echoFile;
        const int trackId = extractFileIdFromPath(echoFile);
        if (trackId > 0) {
            host->cfg_.track_idx_range = buildTrackIdxRange(trackId, host->cfg_.track_idx_window);
        }
    }

    if (host->cfg_.track_idx_range.empty()) {
        std::cerr << "track_idx_range is empty after GMTI file override" << std::endl;
        return false;
    }

    if (host->cfg_.INFO_Type) {
        host->cfg_.pkg_bytes = 256 + 4096 * 16;
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
    host->cfg_.fd_res = host->cfg_.PRF / double(host->cfg_.pulse_num);
    host->cfg_.az_st = 1;
    host->cfg_.az_ed = host->cfg_.pulse_num;
    host->cfg_.az_center = (host->cfg_.az_st + host->cfg_.az_ed) / 2;
    host->cfg_.Loc = true;

    std::vector<std::vector<double>> posMatrix;
    int POS_num = 0;
    if (!host->cfg_.INFO_Type) {
        if (!host->gmti_proc_.POS_dataread(host->cfg_.Plane_POS_add, posMatrix, POS_num)) {
            std::cerr << "[ERR] 读取 POS 失败: " << host->cfg_.Plane_POS_add << std::endl;
            return false;
        }
    }

    GMTIOutput::Plane plane;
    if (!host->gmti_proc_.extractPlanePV(posMatrix, host->cfg_, plane)) {
        std::cerr << "[ERR] 飞机位置提取失败" << std::endl;
        return false;
    }

    std::vector<int> periodList;
    if (!host->gmti_proc_.makePeriodList(plane, host->cfg_, periodList)) {
        std::cerr << "[ERR] 生成 periodList 失败" << std::endl;
        return false;
    }

    if (host->cfg_.track_idx_range.empty()) {
        host->cfg_.track_idx_range = periodList;
    }

    std::vector<GMTIOutput> periodResults;
    const std::vector<std::vector<double>> emptyPos;
    if (!host->gmti_proc_.processPeriodsParallel(periodList,
                                                 host->cfg_,
                                                 host->cfg_.INFO_Type ? emptyPos : posMatrix,
                                                 periodResults)) {
        std::cerr << "[ERR] GMTI period processing failed" << std::endl;
    }

    std::vector<double> MT_acc;
    for (size_t i = 0; i < periodResults.size(); ++i) {
        auto& res = periodResults[i];
        if (!res.MT.empty()) {
            MT_acc.insert(MT_acc.end(), std::make_move_iterator(res.MT.begin()),
                          std::make_move_iterator(res.MT.end()));
        }
        std::cout << "[OK] Period " << periodList[i] << " 完成" << std::endl;
    }

    if (!MT_acc.empty()) {
        if (!host->gmti_proc_.writeResult(MT_acc, host->cfg_)) {
            std::cerr << "[ERR] 整周期检测结果写盘失败" << std::endl;
        }
    } else {
        std::cerr << "[WARN] 本周期无检测目标，跳过检测结果写盘" << std::endl;
    }

    host->latest_gmti_targets_ = trackModule(host->cfg_);
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

} // namespace

MainCtrl::MainCtrl(const std::string& pipeRootPath)
    : pipe_root_path_(pipeRootPath.empty() ? cfg_.pipe_root_path : pipeRootPath)
{
    std::cout << "GMTI INIT..." << std::endl;

    m_SendBuf = static_cast<char*>(std::malloc(40 * 1024 * 1024));
    IsResultReady = false;

    const std::string cmdPipePath = composePipePath(pipe_root_path_, "pipegmticmd");
    const std::string echoPipePath = composePipePath(pipe_root_path_, "pipegmtiecho");
    const std::string resultPipePath = composePipePath(pipe_root_path_, "pipegmtiresult");

    m_RecvCmdPipe.CreatePipe(cmdPipePath.c_str(), O_RDWR);
    m_RecvEchoPipe.CreatePipe(echoPipePath.c_str(), O_RDWR);
    m_SendResultPipe.CreatePipe(resultPipePath.c_str(), O_RDWR);

    m_workmode = Mode_INIT;
    InitThread();
}

MainCtrl::~MainCtrl()
{
    if (m_SendBuf) {
        std::free(m_SendBuf);
        m_SendBuf = nullptr;
    }
}

void MainCtrl::InitThread()
{
    std::cout << "GMTI mode..." << std::endl;

    pthread_t thrRecvCmd;
    pthread_create(&thrRecvCmd, NULL, OnRecvCmdThread, this);

    pthread_t thrRecvEcho;
    pthread_create(&thrRecvEcho, NULL, OnEchoRecvThread, this);

    pthread_t thrProcData;
    pthread_create(&thrProcData, NULL, ProcessDataThread, this);

    pthread_t thrSendRes;
    pthread_create(&thrSendRes, NULL, OnResSendThread, this);
}

void* MainCtrl::OnRecvCmdThread(void *param)
{
    MainCtrl *pHost = static_cast<MainCtrl*>(param);
    ModeSwitchCmd cmd;

    std::cout << "sizeof ModeSwitchCmd = " << sizeof(ModeSwitchCmd) << std::endl;

    while(true) {
        const int ret = pHost->m_RecvCmdPipe.ReadData(reinterpret_cast<char*>(&cmd), sizeof(ModeSwitchCmd));
        std::cout << "Recv cmd ret = " << ret << std::endl;

        if(ret < 0) {
            printf("Receive mode switch command error\n");
            continue;
        }

        printModeSwitchCmd(cmd);
        pHost->last_cmd_ = cmd;

        if(pHost->m_workmode == cmd.workMode)
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
            pHost->gmti_config_xml_ = "temp_config.xml";
        }
        else
        {
            printf("Switch to mode %d\n", cmd.workMode);
        }

        if (!pHost->gmti_config_xml_.empty()) {
            if (!pHost->updateXmlFromModeSwitchCmd(cmd, pHost->gmti_config_xml_)) {
                std::cerr << "Failed to update XML config from ModeSwitchCmd: "
                          << pHost->gmti_config_xml_ << std::endl;
            }
        }

        pHost->m_workmode = cmd.workMode;
    }
}

void* MainCtrl::OnEchoRecvThread(void *param)
{
    MainCtrl *pHost = static_cast<MainCtrl*>(param);

    std::cout << "OnEchoRecvThread start" << std::endl;

    char buf[512] = {0};
    while(true)
    {
        const int ret = pHost->m_RecvEchoPipe.ReadData(buf, 512);
        if(ret < 1) {
            printf("Receive echo data error\n");
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
}

void* MainCtrl::ProcessDataThread(void* param)
{
    MainCtrl *pHost = static_cast<MainCtrl*>(param);

    std::cout << "ProcessDataThread start" << std::endl;

    while(true) {
        std::string fileName;
        if (!popPendingEchoFile(fileName))
        {
            usleep(100);
            continue;
        }

        std::cout << "open file: " << fileName << std::endl;

        const bool ok = runGMTIProcessingFlow(pHost, fileName);
        pHost->IsResultReady = ok;

        if (pHost->m_workmode != Mode_GMTI)
        {
            clearPendingEchoQueue();
            std::cout << "Standby mode: cleared pending echo queue after current file" << std::endl;
        }
    }
}

void* MainCtrl::OnResSendThread(void* param)
{
    MainCtrl *pHost = static_cast<MainCtrl*>(param);

    std::cout << "OnResSendThread start" << std::endl;

    while(true)
    {
        if(pHost->IsResultReady == false)
        {
            sleep(1);
            continue;
        }

        std::cout << "GMTI mode: Sending protocol packet..." << std::endl;

        uint32_t packed_len = 0;
        if (!pHost->packGMTIResults(pHost->m_SendBuf, 40U * 1024U * 1024U, packed_len))
        {
            std::cerr << "Failed to pack GMTI protocol packet" << std::endl;
            pHost->IsResultReady = false;
            continue;
        }

        const ssize_t written = pHost->m_SendResultPipe.WriteData(pHost->m_SendBuf, static_cast<int>(packed_len));
        if (written == static_cast<ssize_t>(packed_len))
        {
            std::cout << "GMTI protocol packet sent: " << packed_len << " bytes" << std::endl;
            std::lock_guard<std::mutex> lock(pHost->result_mutex_);
            pHost->latest_gmti_targets_.clear();
            pHost->latest_gmti_targets_.shrink_to_fit();
            pHost->IsResultReady = false;
        }
        else
        {
            std::cerr << "Failed to send GMTI protocol packet" << std::endl;
        }
        std::cout << "end send res" << std::endl;
    }
}

bool MainCtrl::updateXmlFromModeSwitchCmd(const ModeSwitchCmd& cmd, const std::string& xmlPath)
{
    TiXmlDocument doc(xmlPath.c_str());
    if (!doc.LoadFile())
    {
        std::cerr << "Failed to load XML file: " << xmlPath << std::endl;
        return false;
    }

    TiXmlElement* root = doc.RootElement();
    if (!root)
    {
        std::cerr << "Invalid XML structure: no root element" << std::endl;
        return false;
    }

    if (cmd.CenterFreq > 0)
    {
        if (TiXmlElement* fc_elem = root->FirstChildElement("fc")) {
            fc_elem->SetDoubleAttribute("value", cmd.CenterFreq / 1e9);
        }
    }

    if (cmd.SamplingRate > 0)
    {
        if (TiXmlElement* fs_elem = root->FirstChildElement("fs")) {
            fs_elem->SetDoubleAttribute("value", cmd.SamplingRate / 1e6);
        }
    }

    if (cmd.prf > 0)
    {
        if (TiXmlElement* prf_elem = root->FirstChildElement("PRF")) {
            prf_elem->SetDoubleAttribute("value", cmd.prf);
        }
    }

    if (cmd.BandWidth > 0)
    {
        const double pulse_width = 1.0 / cmd.BandWidth;
        if (TiXmlElement* tr_elem = root->FirstChildElement("Tr")) {
            tr_elem->SetDoubleAttribute("value", pulse_width);
        }
    }

    if (cmd.alt_scene > 0)
    {
        if (TiXmlElement* alt_elem = root->FirstChildElement("alt_scene")) {
            alt_elem->SetDoubleAttribute("value", cmd.alt_scene);
        }
    }

    if (cmd.Rmin > 0)
    {
        if (TiXmlElement* rmin_elem = root->FirstChildElement("R_min")) {
            rmin_elem->SetDoubleAttribute("value", cmd.Rmin);
        }
    }

    if (cmd.velocity > 0)
    {
        if (TiXmlElement* v_elem = root->FirstChildElement("v_platform")) {
            v_elem->SetDoubleAttribute("value", cmd.velocity);
        }
    }

    if (cmd.targetLon != 0)
    {
        if (TiXmlElement* lon_elem = root->FirstChildElement("targetLon")) {
            lon_elem->SetDoubleAttribute("value", cmd.targetLon);
        }
    }

    if (cmd.targetLat != 0)
    {
        if (TiXmlElement* lat_elem = root->FirstChildElement("targetLat")) {
            lat_elem->SetDoubleAttribute("value", cmd.targetLat);
        }
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

bool MainCtrl::packGMTIResults(char* buffer, size_t buffer_size, uint32_t& packed_len)
{
    constexpr size_t kTargetPacketSize = 28U;

    if (!buffer || buffer_size < sizeof(ResultHeader))
    {
        std::cerr << "Invalid buffer for packing results" << std::endl;
        return false;
    }

    std::vector<GMTIDetection> targets;
    {
        std::lock_guard<std::mutex> lock(result_mutex_);
        targets = latest_gmti_targets_;
    }

    const size_t target_count = std::min<size_t>(targets.size(), 65535U);
    const size_t total_bytes = sizeof(ResultHeader) + target_count * kTargetPacketSize;
    if (total_bytes > buffer_size) {
        std::cerr << "Buffer overflow: GMTI packet too large" << std::endl;
        return false;
    }

    std::vector<uint8_t> pkt(total_bytes, 0U);
    ResultHeader header{};
    header.head = 0xAA55;
    header.msgLen = static_cast<uint32_t>(total_bytes - sizeof(uint16_t) - sizeof(uint32_t));
    header.msgAddr = 0;
    header.msgType = 0;
    header.msgCount = next_result_msg_count_++;
    header.srcId = 0;
    header.dstId = 0;
    header.cmdType = static_cast<uint8_t>(Mode_GMTI);
    header.cmdCount = static_cast<uint8_t>(header.msgCount & 0xFFU);
    header.height = 0;
    header.width = 0;
    header.availFlag = 0xFFFF;
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
    header.lonLeftTop = 0;
    header.lonLeftDown = 0;
    header.lonRightDown = 0;
    header.lonRightTop = 0;
    header.lonCenter = 0;
    header.latLeftTop = 0;
    header.latLeftDown = 0;
    header.latRightDown = 0;
    header.latRightTop = 0;
    header.latCenter = 0;
    header.rLeftTop = 0;
    header.rLeftDown = 0;
    header.rRightDown = 0;
    header.rRightTop = 0;
    header.rCenter = 0;
    header.reserve1 = 0;
    header.pixelPitch = 0;
    header.LookDownAngle = 0;
    header.SquintAngle = static_cast<uint16_t>(std::lround(std::fabs(cfg_.squint_angle) * 100.0));
    header.LookSide = (cfg_.squint_side == 0) ? 0x00 : 0xFF;
    std::memset(header.reserve2, 0, sizeof(header.reserve2));
    header.checksum = 0;
    header.targetNum = static_cast<uint16_t>(target_count);

    std::memcpy(pkt.data(), &header, sizeof(header));

    size_t off = sizeof(ResultHeader);
    for (size_t i = 0; i < target_count; ++i) {
        write_target_packet(pkt, off, targets[i]);
        off += kTargetPacketSize;
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
