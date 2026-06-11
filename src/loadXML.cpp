#include "GMTIProcessor.hpp"
#include <iostream>
#include <string>
#include <fstream>
#include <cassert>
#include <sstream>
#include <cctype>
#include <regex>
#include <algorithm>
#include <stdexcept>
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

    auto getTextContent = [](TiXmlElement *element) -> std::string
    {
        if (element && element->GetText())
        {
            return element->GetText();
        }
        return "";
    };

    auto trimText = [](const std::string& text) -> std::string
    {
        size_t first = 0;
        while (first < text.size() && std::isspace(static_cast<unsigned char>(text[first]))) {
            ++first;
        }
        size_t last = text.size();
        while (last > first && std::isspace(static_cast<unsigned char>(text[last - 1]))) {
            --last;
        }
        return text.substr(first, last - first);
    };

    auto getOptionalText = [&](const char *name) -> std::string
    {
        return trimText(getTextContent(param->FirstChildElement(name)));
    };

    auto getRequiredText = [&](const char *name) -> std::string
    {
        const std::string text = getOptionalText(name);
        if (text.empty()) {
            std::ostringstream oss;
            oss << "[XML][ERR] Missing required field <" << name << "> in " << xmlFile;
            throw std::runtime_error(oss.str());
        }
        return text;
    };

    auto parseRequiredInt = [&](const char *name) -> int
    {
        const std::string text = getRequiredText(name);
        try {
            size_t pos = 0;
            const int value = std::stoi(text, &pos);
            if (pos != text.size()) {
                throw std::invalid_argument("trailing characters");
            }
            return value;
        } catch (const std::exception& e) {
            std::ostringstream oss;
            oss << "[XML][ERR] Invalid integer field <" << name << ">=\"" << text
                << "\" in " << xmlFile << ": " << e.what();
            throw std::runtime_error(oss.str());
        }
    };

    auto parseRequiredDouble = [&](const char *name) -> double
    {
        const std::string text = getRequiredText(name);
        try {
            size_t pos = 0;
            const double value = std::stod(text, &pos);
            if (pos != text.size()) {
                throw std::invalid_argument("trailing characters");
            }
            return value;
        } catch (const std::exception& e) {
            std::ostringstream oss;
            oss << "[XML][ERR] Invalid numeric field <" << name << ">=\"" << text
                << "\" in " << xmlFile << ": " << e.what();
            throw std::runtime_error(oss.str());
        }
    };

    auto parseOptionalInt = [&](const char *name, int current) -> int
    {
        const std::string text = getOptionalText(name);
        if (text.empty()) {
            return current;
        }
        try {
            size_t pos = 0;
            const int value = std::stoi(text, &pos);
            if (pos != text.size()) {
                throw std::invalid_argument("trailing characters");
            }
            return value;
        } catch (const std::exception& e) {
            std::ostringstream oss;
            oss << "[XML][ERR] Invalid optional integer field <" << name << ">=\"" << text
                << "\" in " << xmlFile << ": " << e.what();
            throw std::runtime_error(oss.str());
        }
    };

    auto parseOptionalDouble = [&](const char *name, double current) -> double
    {
        const std::string text = getOptionalText(name);
        if (text.empty()) {
            return current;
        }
        try {
            size_t pos = 0;
            const double value = std::stod(text, &pos);
            if (pos != text.size()) {
                throw std::invalid_argument("trailing characters");
            }
            return value;
        } catch (const std::exception& e) {
            std::ostringstream oss;
            oss << "[XML][ERR] Invalid optional numeric field <" << name << ">=\"" << text
                << "\" in " << xmlFile << ": " << e.what();
            throw std::runtime_error(oss.str());
        }
    };

    auto parseOptionalSizeT = [&](const char *name, size_t current) -> size_t
    {
        const std::string text = getOptionalText(name);
        if (text.empty()) {
            return current;
        }
        try {
            size_t pos = 0;
            const unsigned long long value = std::stoull(text, &pos);
            if (pos != text.size()) {
                throw std::invalid_argument("trailing characters");
            }
            return value > 0 ? static_cast<size_t>(value) : current;
        } catch (const std::exception& e) {
            std::ostringstream oss;
            oss << "[XML][ERR] Invalid optional size field <" << name << ">=\"" << text
                << "\" in " << xmlFile << ": " << e.what();
            throw std::runtime_error(oss.str());
        }
    };

    auto parseOptionalBool = [&](const char *name, bool current) -> bool
    {
        const std::string text = getOptionalText(name);
        if (text.empty()) {
            return current;
        }
        if (text == "true" || text == "1") {
            return true;
        }
        if (text == "false" || text == "0") {
            return false;
        }
        std::ostringstream oss;
        oss << "[XML][ERR] Invalid boolean field <" << name << ">=\"" << text
            << "\" in " << xmlFile << ". Expected true/false/1/0.";
        throw std::runtime_error(oss.str());
    };

    auto parseOptionalString = [&](const char *name, const std::string& current) -> std::string
    {
        const std::string text = getOptionalText(name);
        return text.empty() ? current : text;
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
                try {
                    values.push_back(std::stoi(trimmed));
                } catch (const std::exception& e) {
                    std::ostringstream oss;
                    oss << "[XML][ERR] Invalid integer in <track_idx_range>: \""
                        << trimmed << "\": " << e.what();
                    throw std::runtime_error(oss.str());
                }
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

    try {
    cfg.GMTI_Data_add = getOptionalText("GMTI_data");
    cfg.GMTI_Data_add2 = getOptionalText("GMTI_data2");
    cfg.GMTI_Data_new = getOptionalText("GMTI_data_new");
    cfg.result_add = getRequiredText("result_add");
    {
        const std::string pipeRootPath = getOptionalText("pipe_root_path");
        if (!pipeRootPath.empty()) {
            cfg.pipe_root_path = pipeRootPath;
        }
    }
    cfg.Plane_POS_add = getOptionalText("Plane_POS");
    cfg.CAR_POS_add = getOptionalText("CAR_POS");
    cfg.reffunc_add = getOptionalText("reffunc_add");
    cfg.channel_mode = getRequiredText("isSeparated");
    cfg.INFO_Type = parseRequiredInt("INFO_Type");
    cfg.iq_compose = getOptionalText("iq_compose");

    // 解析数字字段并填充到结构体中
    cfg.isPC = parseRequiredInt("isPC");
    cfg.hasRefFunc = parseRequiredInt("hasRefFunc");
    cfg.info_len = parseRequiredInt("info_len");
    cfg.pulse_len = parseRequiredInt("pulse_len");
    cfg.rg_len = parseRequiredInt("rg_len");
    cfg.pulse_num = parseRequiredInt("pulse_num");
    cfg.read_pulse_num = parseOptionalInt("read_pulse_num", cfg.read_pulse_num);
    cfg.read_pulse_offset = parseOptionalInt("read_pulse_offset", cfg.read_pulse_offset);
    cfg.pulse_dec = parseRequiredInt("pulse_dec");
    cfg.fc = parseRequiredDouble("fc") * 1e9;  // 转换为赫兹
    cfg.Br = parseRequiredDouble("Br") * 1e6;  // 转换为赫兹
    cfg.fs = parseRequiredDouble("fs") * 1e6;  // 转换为赫兹
    cfg.Tr = parseRequiredDouble("Tr") * 1e-6; // 转换为秒
    cfg.PRF = parseRequiredDouble("PRF");
    cfg.az_count = parseRequiredInt("az_count");
    cfg.beamwidth_deg = parseOptionalDouble("boshu", cfg.beamwidth_deg);
    cfg.loc_beam_gate_deg = parseOptionalDouble("loc_beam_gate_deg", cfg.loc_beam_gate_deg);
    cfg.week = parseRequiredInt("week_offset");
    cfg.d_channel = parseRequiredDouble("d_chan");
    cfg.pf = parseRequiredDouble("pf");
    cfg.R_min = parseRequiredDouble("Rmin");
    cfg.L0 = parseRequiredDouble("ref_lon");
    cfg.MT_nowz = parseRequiredDouble("ref_H");
    cfg.secBias = parseRequiredInt("secBias");
    cfg.skip_az_num = parseRequiredInt("skip_pulses");
    cfg.calib_coef = parseRequiredDouble("calib_coef");

    if (cfg.GMTI_Data_new.empty()) {
        cfg.GMTI_Data_new = cfg.GMTI_Data_add;
    }

    cfg.wavepos_st = parseRequiredInt("wavepos_st");
    cfg.wavepos_ed = parseRequiredInt("wavepos_ed");
    cfg.wavepos_skip = parseRequiredInt("wavepos_skip");
    cfg.scan_min_deg = parseOptionalDouble("scan_min_deg", cfg.scan_min_deg);
    cfg.scan_max_deg = parseOptionalDouble("scan_max_deg", cfg.scan_max_deg);
    cfg.min_points = parseRequiredInt("min_points");
    cfg.min_len = parseRequiredInt("min_len");

    cfg.rg_st = parseRequiredInt("rg_st");
    cfg.rg_ed = parseRequiredInt("rg_ed");
    cfg.squint_side = parseRequiredInt("squint_side");

    cfg.roi_ll_deg[0] = parseRequiredDouble("roi_lat1");
    cfg.roi_ll_deg[1] = parseRequiredDouble("roi_lng1");
    cfg.roi_ll_deg[2] = parseRequiredDouble("roi_Ry1");
    cfg.roi_ll_deg[3] = parseRequiredDouble("roi_Rx1");

    cfg.wavepos_use_roi = parseOptionalBool("wavepos_use_roi", cfg.wavepos_use_roi);

    cfg.wavepos_parallel = parseOptionalBool("wavepos_parallel", cfg.wavepos_parallel);

    cfg.enable_dbs_fusion = parseOptionalBool("enable_dbs_fusion", cfg.enable_dbs_fusion);

    {
        cfg.dbs_out_res_m = parseOptionalDouble("raw_fenbianlv", cfg.dbs_out_res_m);
        cfg.dbs_beam_skip = std::max(1, parseOptionalInt("n_tiaoguo", cfg.dbs_beam_skip));
        cfg.dbs_range_skip = std::max(1, parseOptionalInt("len_tiaoguo", cfg.dbs_range_skip));
        cfg.dbs_interp_mode = parseOptionalInt("dbs_interp_mode", cfg.dbs_interp_mode);
        cfg.dbs_max_mosaic_pixels = parseOptionalSizeT("dbs_max_mosaic_pixels",
                                                       cfg.dbs_max_mosaic_pixels);
    }

    cfg.estimate_error_angle = parseOptionalBool("estimate_error_angle", cfg.estimate_error_angle);

    cfg.lat_st = parseRequiredDouble("lat_st");
    cfg.lat_ed = parseRequiredDouble("lat_ed");
    cfg.lon_st = parseRequiredDouble("lon_st");
    cfg.lon_ed = parseRequiredDouble("lon_ed");
    cfg.squint_angle = parseOptionalDouble("squint_angle", cfg.squint_angle);

    // 可选调试配置：节点缺失时使用结构体默认值
    cfg.mt_sort_by_utc = parseOptionalBool("mt_sort_by_utc", cfg.mt_sort_by_utc);
    cfg.track_debug_level = parseOptionalInt("track_debug_level", cfg.track_debug_level);
    cfg.track_debug_frames = parseOptionalInt("track_debug_frames", cfg.track_debug_frames);
    cfg.track_debug_points = parseOptionalInt("track_debug_points", cfg.track_debug_points);
    cfg.track_idx_window = parseOptionalInt("track_idx_window", cfg.track_idx_window);
    cfg.track_truth_threshold = parseOptionalInt("track_truth_threshold", cfg.track_truth_threshold);
    cfg.track_confirm_window = parseOptionalInt("track_confirm_window", cfg.track_idx_window);
    cfg.track_confirm_hits = parseOptionalInt("track_confirm_hits", cfg.track_truth_threshold);
    cfg.track_max_missed = parseOptionalInt("track_max_missed", cfg.track_max_missed);
    cfg.track_tentative_max_missed = parseOptionalInt("track_tentative_max_missed",
                                                       cfg.track_tentative_max_missed);
    cfg.track_default_dt = parseOptionalDouble("track_default_dt", cfg.track_default_dt);
    cfg.track_chi2_gate = parseOptionalDouble("track_chi2_gate", cfg.track_chi2_gate);
    cfg.track_tentative_gate_scale = parseOptionalDouble("track_tentative_gate_scale",
                                                         cfg.track_tentative_gate_scale);
    cfg.track_tentative_chi2_scale = parseOptionalDouble("track_tentative_chi2_scale",
                                                         cfg.track_tentative_chi2_scale);
    cfg.track_dummy_cost = parseOptionalDouble("track_dummy_cost", cfg.track_dummy_cost);
    cfg.track_linearity_window = parseOptionalInt("track_linearity_window", cfg.track_linearity_window);
    cfg.track_min_linearity_confirm = parseOptionalDouble("track_min_linearity_confirm",
                                                          cfg.track_min_linearity_confirm);
    cfg.track_speed_smooth_weight = parseOptionalDouble("track_speed_smooth_weight",
                                                        cfg.track_speed_smooth_weight);
    cfg.track_heading_weight = parseOptionalDouble("track_heading_weight", cfg.track_heading_weight);
    cfg.track_process_noise_pos = parseOptionalDouble("track_process_noise_pos",
                                                      cfg.track_process_noise_pos);
    cfg.track_process_noise_vel = parseOptionalDouble("track_process_noise_vel",
                                                      cfg.track_process_noise_vel);
    cfg.track_measurement_noise_pos = parseOptionalDouble("track_measurement_noise_pos",
                                                          cfg.track_measurement_noise_pos);
    cfg.track_debug_dump = parseOptionalBool("track_debug_dump", cfg.track_debug_dump);
    cfg.track_debug_dir = parseOptionalString("track_debug_dir", cfg.track_debug_dir);
    cfg.track_debug_dump_level = parseOptionalInt("track_debug_dump_level",
                                                  cfg.track_debug_dump_level);

    const std::string trackIdxRangeStr = getOptionalText("track_idx_range");
    if (!trackIdxRangeStr.empty()) {
        cfg.track_idx_range = parseIntList(trackIdxRangeStr);
    }

    cfg.track_gate_m = parseOptionalDouble("track_gate_m", cfg.track_gate_m);
    cfg.track_v_max = parseOptionalDouble("track_v_max", cfg.track_v_max);

    if (cfg.track_idx_range.empty()) {
        const int trackId = extractFileId(cfg.GMTI_Data_new.empty() ? cfg.GMTI_Data_add : cfg.GMTI_Data_new);
        if (trackId > 0) {
            cfg.result_file_id = trackId;
        } else {
            std::string nextPath;
            int nextId = -1;
            if (!nextGMTIFileName(cfg.result_add, nextPath, nextId)) {
                std::cerr << "[XML][ERR] Cannot allocate incremental GMTI result id from result_add: "
                          << cfg.result_add << std::endl;
                return false;
            }
            cfg.result_file_id = nextId;
            std::cout << "[XML][WARN] Cannot extract file id from echo path. Use next result id "
                      << nextId << " for GMTI result and track window." << std::endl;
        }
        cfg.track_idx_range = buildTrackIdxRange(cfg.result_file_id, cfg.track_idx_window);
    }
    if (cfg.result_file_id <= 0) {
        const int trackId = extractFileId(cfg.GMTI_Data_new.empty() ? cfg.GMTI_Data_add : cfg.GMTI_Data_new);
        if (trackId > 0) {
            cfg.result_file_id = trackId;
        } else {
            std::string nextPath;
            int nextId = -1;
            if (!nextGMTIFileName(cfg.result_add, nextPath, nextId)) {
                std::cerr << "[XML][ERR] Cannot allocate incremental GMTI result id from result_add: "
                          << cfg.result_add << std::endl;
                return false;
            }
            cfg.result_file_id = nextId;
            std::cout << "[XML][WARN] Cannot extract file id from echo path. Use next result id "
                      << nextId << " for GMTI result." << std::endl;
        }
    }
    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        return false;
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
                                   std::vector<std::complex<float>>& data1,
                                   std::vector<std::complex<float>>& data2,
                                   std::vector<double>& utc,
                                   double& theta_sq)
{
    TIMING_SCOPE(readPulseBlock);
    const char* filepath = cfg.GMTI_Data_add2.c_str();
    std::vector<double> fw_angle_deg;

    if(!readBeamRawFloat(cfg, filepath, beamskip, data2, fw_angle_deg, utc)) {
        std::cerr << "读取回波数据失败" << std::endl;
        return false;
    }

    const char* filepath2 = cfg.GMTI_Data_add.c_str();
    if(!readBeamRawFloat(cfg, filepath2, beamskip, data1, fw_angle_deg, utc)) {
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
