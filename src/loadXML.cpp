#include "GMTIProcessor.hpp"
#include <iostream>
#include <string>
#include <fstream>
#include <cassert>
#include <sstream>
#include <cctype>
#include <regex>
#include "rangeCompress.hpp"

// 读取 XML 配置文件并填充到 Config 结构体中
bool GMTIProcessor::readXmlParam(const std::string &xmlFile, Config &cfg)
{
    TiXmlDocument doc;
    if (!doc.LoadFile(xmlFile.c_str()))
    {
        std::cerr << "Error loading XML file: " << xmlFile;
        if (doc.Error())
        {
            std::cerr << " (" << doc.ErrorDesc() << ")";
        }
        std::cerr << std::endl;
        return false; // 如果加载失败，返回 false
    }

    // 访问根元素
    TiXmlElement *root = doc.FirstChildElement("GMTI");
    if (!root)
    {
        std::cerr << "Root element <GMTI> not found." << std::endl;
        return false; // 如果根元素不存在，返回 false
    }

    // 访问 GMTI_parameter 元素
    TiXmlElement *param = root->FirstChildElement("GMTI_parameter");
    if (!param)
    {
        std::cerr << "GMTI_parameter element not found." << std::endl;
        return false; // 如果没有找到 GMTI_parameter 元素，返回 false
    }

    // 提取文本内容的辅助函数
    auto getTextContent = [](TiXmlElement *element) -> std::string
    {
        if (element && element->GetText())
        {
            return element->GetText();
        }
        return "";
    };

    auto parseIntList = [](const std::string& text) -> std::vector<int>
    {
        std::vector<int> values;
        std::string token;
        std::stringstream ss(text);
        while (std::getline(ss, token, ',')) {
            std::string trimmed;
            for (char ch : token) {
                if (!std::isspace(static_cast<unsigned char>(ch))) {
                    trimmed.push_back(ch);
                }
            }
            if (!trimmed.empty()) {
                values.push_back(std::stoi(trimmed));
            }
        }
        return values;
    };

    auto extractFileId = [](const std::string& path) -> int
    {
        const size_t pos = path.find_last_of("/\\");
        const std::string filename = (pos == std::string::npos) ? path : path.substr(pos + 1);
        std::smatch match;
        const std::regex idPattern(R"((\d+)(?:\.bin)?$)", std::regex::icase);
        if (std::regex_search(filename, match, idPattern) && match.size() >= 2) {
            return std::stoi(match[1].str());
        }
        return -1;
    };

    auto buildTrackIdxRange = [](int id, int window) -> std::vector<int>
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
    };

    // 从 XML 读取并填充 cfg 结构体
    cfg.GMTI_Data_add = getTextContent(param->FirstChildElement("GMTI_data"));
    cfg.GMTI_Data_add2 = getTextContent(param->FirstChildElement("GMTI_data2"));
    cfg.GMTI_Data_new = getTextContent(param->FirstChildElement("GMTI_data_new"));
    cfg.result_add = getTextContent(param->FirstChildElement("result_add"));
    {
        const std::string pipeRootPath = getTextContent(param->FirstChildElement("pipe_root_path"));
        if (!pipeRootPath.empty()) {
            cfg.pipe_root_path = pipeRootPath;
        }
    }
    cfg.Plane_POS_add = getTextContent(param->FirstChildElement("Plane_POS"));
    cfg.CAR_POS_add = getTextContent(param->FirstChildElement("CAR_POS"));
    cfg.reffunc_add = getTextContent(param->FirstChildElement("reffunc_add"));
    cfg.channel_mode = getTextContent(param->FirstChildElement("isSeparated"));
    cfg.INFO_Type = std::stoi(getTextContent(param->FirstChildElement("INFO_Type")));
    cfg.iq_compose = getTextContent(param->FirstChildElement("iq_compose"));

    // 解析数字字段并填充到结构体中
    cfg.isPC = std::stoi(getTextContent(param->FirstChildElement("isPC")));
    cfg.hasRefFunc = std::stoi(getTextContent(param->FirstChildElement("hasRefFunc")));
    cfg.info_len = std::stoi(getTextContent(param->FirstChildElement("info_len")));
    cfg.pulse_len = std::stoi(getTextContent(param->FirstChildElement("pulse_len")));
    cfg.rg_len = std::stoi(getTextContent(param->FirstChildElement("rg_len")));
    cfg.pulse_num = std::stoi(getTextContent(param->FirstChildElement("pulse_num")));
    cfg.pulse_dec = std::stoi(getTextContent(param->FirstChildElement("pulse_dec")));
    cfg.fc = std::stod(getTextContent(param->FirstChildElement("fc"))) * 1e9;  // 转换为赫兹
    cfg.Br = std::stod(getTextContent(param->FirstChildElement("Br"))) * 1e6;  // 转换为赫兹
    cfg.fs = std::stod(getTextContent(param->FirstChildElement("fs"))) * 1e6;  // 转换为赫兹
    cfg.Tr = std::stod(getTextContent(param->FirstChildElement("Tr"))) * 1e-6; // 转换为秒
    cfg.PRF = std::stod(getTextContent(param->FirstChildElement("PRF")));
    cfg.az_count = std::stoi(getTextContent(param->FirstChildElement("az_count")));
    cfg.week = std::stoi(getTextContent(param->FirstChildElement("week_offset")));
    cfg.d_channel = std::stod(getTextContent(param->FirstChildElement("d_chan")));
    cfg.pf = std::stod(getTextContent(param->FirstChildElement("pf")));
    cfg.R_min = std::stod(getTextContent(param->FirstChildElement("Rmin")));
    cfg.L0 = std::stod(getTextContent(param->FirstChildElement("ref_lon")));
    cfg.MT_nowz = std::stod(getTextContent(param->FirstChildElement("ref_H")));
    cfg.secBias = std::stoi(getTextContent(param->FirstChildElement("secBias")));
    cfg.skip_az_num = std::stoi(getTextContent(param->FirstChildElement("skip_pulses")));
    cfg.calib_coef = std::stod(getTextContent(param->FirstChildElement("calib_coef")));

    if (cfg.GMTI_Data_new.empty()) {
        cfg.GMTI_Data_new = cfg.GMTI_Data_add;
    }

    cfg.wavepos_st = std::stoi(getTextContent(param->FirstChildElement("wavepos_st")));
    cfg.wavepos_ed = std::stoi(getTextContent(param->FirstChildElement("wavepos_ed")));
    cfg.wavepos_skip = std::stoi(getTextContent(param->FirstChildElement("wavepos_skip")));
    cfg.min_points = std::stoi(getTextContent(param->FirstChildElement("min_points")));
    cfg.min_len = std::stoi(getTextContent(param->FirstChildElement("min_len")));

    cfg.rg_st = std::stoi(getTextContent(param->FirstChildElement("rg_st")));
    cfg.rg_ed = std::stoi(getTextContent(param->FirstChildElement("rg_ed")));
    cfg.squint_side = std::stoi(getTextContent(param->FirstChildElement("squint_side")));

    cfg.roi_ll_deg[0] = std::stod(getTextContent(param->FirstChildElement("roi_lat1")));
    cfg.roi_ll_deg[1] = std::stod(getTextContent(param->FirstChildElement("roi_lng1")));
    cfg.roi_ll_deg[2] = std::stod(getTextContent(param->FirstChildElement("roi_Ry1")));
    cfg.roi_ll_deg[3] = std::stod(getTextContent(param->FirstChildElement("roi_Rx1")));

    cfg.wavepos_use_roi = (getTextContent(param->FirstChildElement("wavepos_use_roi")) == "true");

    const std::string waveposParallelStr = getTextContent(param->FirstChildElement("wavepos_parallel"));
    if (!waveposParallelStr.empty()) {
        cfg.wavepos_parallel = (waveposParallelStr == "true" || waveposParallelStr == "1");
    }

    const std::string enableDbsFusionStr = getTextContent(param->FirstChildElement("enable_dbs_fusion"));
    if (!enableDbsFusionStr.empty()) {
        cfg.enable_dbs_fusion = (enableDbsFusionStr == "true" || enableDbsFusionStr == "1");
    }

    const std::string estimateErrorAngleStr = getTextContent(param->FirstChildElement("estimate_error_angle"));
    if (!estimateErrorAngleStr.empty()) {
        cfg.estimate_error_angle = (estimateErrorAngleStr == "true" || estimateErrorAngleStr == "1");
    }

    cfg.lat_st = std::stod(getTextContent(param->FirstChildElement("lat_st")));
    cfg.lat_ed = std::stod(getTextContent(param->FirstChildElement("lat_ed")));
    cfg.lon_st = std::stod(getTextContent(param->FirstChildElement("lon_st")));
    cfg.lon_ed = std::stod(getTextContent(param->FirstChildElement("lon_ed")));
    {
        const std::string squintAngleStr = getTextContent(param->FirstChildElement("squint_angle"));
        if (!squintAngleStr.empty()) {
            cfg.squint_angle = std::stod(squintAngleStr);
        }
    }

    // 可选调试配置：节点缺失时使用结构体默认值
    const std::string mtSortByUtcStr = getTextContent(param->FirstChildElement("mt_sort_by_utc"));
    if (!mtSortByUtcStr.empty()) {
        cfg.mt_sort_by_utc = (mtSortByUtcStr == "true" || mtSortByUtcStr == "1");
    }

    const std::string trackDebugLevelStr = getTextContent(param->FirstChildElement("track_debug_level"));
    if (!trackDebugLevelStr.empty()) {
        cfg.track_debug_level = std::stoi(trackDebugLevelStr);
    }

    const std::string trackDebugFramesStr = getTextContent(param->FirstChildElement("track_debug_frames"));
    if (!trackDebugFramesStr.empty()) {
        cfg.track_debug_frames = std::stoi(trackDebugFramesStr);
    }

    const std::string trackDebugPointsStr = getTextContent(param->FirstChildElement("track_debug_points"));
    if (!trackDebugPointsStr.empty()) {
        cfg.track_debug_points = std::stoi(trackDebugPointsStr);
    }

    const std::string trackIdxWindowStr = getTextContent(param->FirstChildElement("track_idx_window"));
    if (!trackIdxWindowStr.empty()) {
        cfg.track_idx_window = std::stoi(trackIdxWindowStr);
    }

    const std::string trackTruthThresholdStr = getTextContent(param->FirstChildElement("track_truth_threshold"));
    if (!trackTruthThresholdStr.empty()) {
        cfg.track_truth_threshold = std::stoi(trackTruthThresholdStr);
    }

    const std::string trackIdxRangeStr = getTextContent(param->FirstChildElement("track_idx_range"));
    if (!trackIdxRangeStr.empty()) {
        cfg.track_idx_range = parseIntList(trackIdxRangeStr);
    }

    const std::string trackGateStr = getTextContent(param->FirstChildElement("track_gate_m"));
    if (!trackGateStr.empty()) {
        cfg.track_gate_m = std::stod(trackGateStr);
    }

    const std::string trackVmaxStr = getTextContent(param->FirstChildElement("track_v_max"));
    if (!trackVmaxStr.empty()) {
        cfg.track_v_max = std::stod(trackVmaxStr);
    }

    if (cfg.track_idx_range.empty()) {
        const int trackId = extractFileId(cfg.GMTI_Data_new.empty() ? cfg.GMTI_Data_add : cfg.GMTI_Data_new);
        cfg.track_idx_range = buildTrackIdxRange(trackId, cfg.track_idx_window);
    }

    return true; // 如果一切顺利，返回 true
}

using namespace std;

// 从 POS 文件读取数据并返回是否成功
bool GMTIProcessor::POS_dataread(const std::string &posFile, std::vector<std::vector<double>> &POS_data, int &POS_num)
{
    // 打开 POS 文件
    const int pos_len = 7;
    std::ifstream posFileStream(posFile, std::ios::binary);
    if (!posFileStream)
    {
        std::cerr << "无法打开 POS 文件: " << posFile << std::endl;
        return false; // 文件打开失败，返回 false
    }

    // 获取文件大小并计算帧数
    posFileStream.seekg(0, std::ios::end);            // 移动到文件末尾
    std::streamsize fileSize = posFileStream.tellg(); // 获取文件大小
    posFileStream.seekg(0, std::ios::beg);            // 移动回文件开头

    // 每一帧数据的字节数：pos_len 个 double，每个 double 占 8 字节
    const int POS_pkg_len = pos_len * sizeof(double);
    // 计算文件中包含的帧数
    POS_num = static_cast<int>(fileSize / POS_pkg_len);

    // 检查文件大小与帧数是否匹配
    if (fileSize % POS_pkg_len != 0)
    {
        std::cerr << "文件大小与帧数不匹配，可能包含损坏数据。" << std::endl;
        return false; // 如果文件大小无法整除每帧字节数，返回 false
    }

    // 读取 POS 数据
    POS_data.resize(POS_num, std::vector<double>(pos_len)); // 根据帧数调整二维 vector 的大小
    for (int i = 0; i < POS_num; ++i)
    {
        for (int j = 0; j < pos_len; ++j)
        {
            posFileStream.read(reinterpret_cast<char *>(&POS_data[i][j]), sizeof(double));
        }
    }

    posFileStream.close();

    DBG("成功读取 POS 数据");

    return true; // 数据读取成功，返回 true
}

// 打开回波文件（兼容单/双文件）
bool GMTIProcessor::openEchoFiles(const Config &cfg, std::ifstream &fid1, std::ifstream &fid2)
{
    if (cfg.channel_mode == "separate")
    { // 双文件模式
        fid1.open(cfg.GMTI_Data_add, std::ios::binary);
        fid2.open(cfg.GMTI_Data_add2, std::ios::binary);
        if (!fid1.is_open() || !fid2.is_open())
        {
            std::cerr << "无法打开回波数据双文件" << std::endl;
            return false;
        }
    }
    else if (cfg.channel_mode == "interleaved")
    { // 单文件模式
        fid1.open(cfg.GMTI_Data_add, std::ios::binary);
        fid2.close(); // 无第二个文件
        if (!fid1.is_open())
        {
            std::cerr << "无法打开交织回波数据文件" << std::endl;
            return false;
        }
    }
    else
    {
        std::cerr << "未知 channel_mode: " << cfg.channel_mode << std::endl;
        return false;
    }

    return true;
}

bool GMTIProcessor::readPulseBlock(const Config& cfg,
                                   int beamskip,
                                   std::vector<std::complex<double>>& data1,
                                   std::vector<std::complex<double>>& data2,
                                   std::vector<double>& utc,
                                   double& theta_sq)
{
    TIMING_SCOPE(readPulseBlock);
    const char* filepath = cfg.GMTI_Data_add2.c_str();
    std::vector<double> fw_angle_deg;

    if(!readBeamRaw(cfg, filepath, beamskip, data2, fw_angle_deg, utc)) {
        std::cerr << "读取回波数据失败" << std::endl;
        return false;
    }

    const char* filepath2 = cfg.GMTI_Data_add.c_str();
    if(!readBeamRaw(cfg, filepath2, beamskip, data1, fw_angle_deg, utc)) {
        std::cerr << "读取回波数据失败" << std::endl;
        return false;
    }

    // 计算方位角平方的平均值
    theta_sq = fw_angle_deg[fw_angle_deg.size() / 2];

    return true;
}




// 提取 UTC 时间（秒）
bool GMTIProcessor::extractUTC(const std::vector<uint8_t>& headerBytes, const Config& cfg,
                                std::vector<double>& t) {
    size_t pulse_num = headerBytes.size() / cfg.info_len;  // 计算脉冲数量（假设headerBytes为一维数组）

    // 校验headerBytes大小是否足够
    assert(headerBytes.size() >= cfg.info_len * pulse_num && "headerBytes 大小不正确");

    double TickHz = 1e8;  // 固定为 1e8
    double DayBias = cfg.week * 24 * 3600;  // 从 cfg 获取 DayBias
    double SecBias = cfg.secBias;  // 从 cfg 获取 SecBias

    // BCD 解码：时分秒
    auto bcd2dec = [](uint8_t u8) -> double {
        return static_cast<double>(u8) - 6.0 * std::floor(static_cast<double>(u8) / 16.0);
    };

    t.resize(pulse_num);  // 分配时间结果向量

    if (cfg.INFO_Type) {  // 如果 INFO_Type 为 true
        for (size_t i = 0; i < pulse_num; ++i) {
            // 从 headerBytes 中提取 BCD 编码的时、分、秒
            double hh = bcd2dec(headerBytes[36 + i * cfg.info_len]);  // 小时 0..23
            double mm = bcd2dec(headerBytes[37 + i * cfg.info_len]);  // 分钟 0..59
            double ss = bcd2dec(headerBytes[38 + i * cfg.info_len]);  // 秒 0..59

            // 32-bit 计数（小端序）→ 秒的小数部分
            double ds = static_cast<double>(headerBytes[42 + i * cfg.info_len]) * 256.0 * 256.0 * 256.0 +
                        static_cast<double>(headerBytes[41 + i * cfg.info_len]) * 256.0 * 256.0 +
                        static_cast<double>(headerBytes[40 + i * cfg.info_len]) * 256.0 +
                        static_cast<double>(headerBytes[39 + i * cfg.info_len]);

            double fracSec = ds / TickHz;

            // 组合为秒
            t[i] = DayBias + SecBias + hh * 3600.0 + mm * 60.0 + ss + fracSec;
        }
    } else {  // 如果 INFO_Type 为 false
        for (size_t i = 0; i < pulse_num; ++i) {
            // 从 headerBytes 中提取 BCD 编码的时、分、秒
            double hh = bcd2dec(headerBytes[241 + i * cfg.info_len]);  // 小时 0..23
            double mm = bcd2dec(headerBytes[242 + i * cfg.info_len]);  // 分钟 0..59
            double ss = bcd2dec(headerBytes[243 + i * cfg.info_len]);  // 秒 0..59

            // 32-bit 计数（小端序）→ 秒的小数部分
            double ds = static_cast<double>(headerBytes[247 + i * cfg.info_len]) * 256.0 * 256.0 * 256.0 +
                        static_cast<double>(headerBytes[246 + i * cfg.info_len]) * 256.0 * 256.0 +
                        static_cast<double>(headerBytes[245 + i * cfg.info_len]) * 256.0 +
                        static_cast<double>(headerBytes[244 + i * cfg.info_len]);

            double fracSec = ds / TickHz;

            // 组合为秒
            t[i] = DayBias + SecBias + hh * 3600.0 + mm * 60.0 + ss + fracSec;
        }
    }

    return true;  // 成功提取 UTC 时间
}

// 读取脉冲块（float版本, GPU路径优化使用）
bool GMTIProcessor::readPulseBlockFloat(const Config& cfg,
                                        int beamskip,
                                        std::vector<std::complex<float>>& data1,
                                        std::vector<std::complex<float>>& data2,
                                        std::vector<double>& utc,
                                        double& theta_sq)
{
    TIMING_SCOPE(readPulseBlockFloat);
    const char* filepath = cfg.GMTI_Data_add2.c_str();
    std::vector<double> fw_angle_deg;

    if(!readBeamRawFloat(cfg, filepath, beamskip, data2, fw_angle_deg, utc)) {
        std::cerr << "读取回波数据失败 (float)" << std::endl;
        return false;
    }

    const char* filepath2 = cfg.GMTI_Data_add.c_str();
    if(!readBeamRawFloat(cfg, filepath2, beamskip, data1, fw_angle_deg, utc)) {
        std::cerr << "读取回波数据失败 (float)" << std::endl;
        return false;
    }

    // 计算方位角平方的平均值
    theta_sq = fw_angle_deg[fw_angle_deg.size() / 2];

    return true;
}
