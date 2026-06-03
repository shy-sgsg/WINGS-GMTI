# GMTI 管道通信整合方案

> 本文档遵照 规范.md 的原则进行规划：不假设、暴露权衡、简单优先、手术性变更、目标驱动

**制作日期**：2025年5月27日  
**状态**：待确认与实施  
**目标受众**：项目经理、开发团队、系统集成

---

## 📋 执行摘要

### 当前问题
- SAR 成像软件 (`src/pipe/`) 具有完整的管道通信框架
- GMTI 算法主程序 (`src/main.cpp`) 使用 XML 硬编码方式
- 两个系统**独立运行**，无法实现实时参数下发和结果回传

### 改动目标
将 `src/pipe/` 的管道通信框架与 GMTI 算法融合，实现：
1. ✓ 通过 Linux 命名管道接收参数，**动态更新 XML 配置文件**（GMTI 程序继续依赖 XML）
2. ✓ 执行 GMTI 处理流程（保留现有算法逻辑）
3. ✓ 通过管道发送二进制结果（Track 数据）**- 具体协议预留接口待后续完善**
4. ✓ 保留向后兼容性（仍可独立用 XML 运行）

### 改动程度：**中等规模** 
- **新增代码**：~300-400 行（相比原计划减少，因为参数直接更新 XML）
- **修改代码**：~200-300 行（主要在 MainCtrl.cpp）
- **删除代码**：~0-30 行（注释/冗余部分）
- **文件数量**：修改 4-5 个文件，新增 0 个文件

> **关键变化**：参数处理流向变为 "ModeSwitchCmd → 更新 XML" 而非 "ModeSwitchCmd → 直接映射 Config"，降低了算法集成的耦合度和风险

---

## 🎯 改动清单

### 🔄 数据流向对比

**原计划（参数直接替代 XML）**：
```
ModeSwitchCmd → 直接映射到 Config → GMTI 处理
                (高耦合，需参数完整映射)
```

**新方案（参数更新 XML）**：
```
ModeSwitchCmd → 写入 XML 文件 → readXmlParam() → Config → GMTI 处理
                (低耦合，保留 XML 中心，风险低)
```

**优势**：
- ✓ 保留了 GMTI 对 XML 的现有依赖，改动最小
- ✓ XML 成为参数的持久化层，便于调试和日志
- ✓ 如果管道通信中断，XML 仍可独立使用
- ✓ 参数映射复杂度从 "1:1 直接映射" 降低到 "选择性更新"

---

### Phase 0：架构分析与数据映射

#### 问题 1：参数映射策略
**新方案**：不是 "完全直接映射"，而是 "选择性更新"

ModeSwitchCmd → XML 的字段映射：
```
ModeSwitchCmd 字段              XML 字段              更新策略
──────────────────────────────────────────────────────────────
CenterFreq                      fc                   ✓ 直接写入
SamplingRate                    fs                   ✓ 直接写入
prf                             PRF                  ✓ 直接写入
BandWidth                       (pulse_width derived)✓ 直接写入
alt_scene                       alt_scene            ✓ 直接写入
Rmin                            R_min                ✓ 直接写入
velocity                        v_platform           ✓ 直接写入
targetLon, targetLat, ...       对应字段             ✓ 直接写入
────────────────────────────────────────────────────────────
[缺失字段]  pulse_num           pulse_num            ⚠️ XML 保留，不更新
[缺失字段]  rg_len              rg_len               ⚠️ XML 保留，不更新
[缺失字段]  channel_mode        channel_mode         ⚠️ XML 保留，不更新
```

**策略**：

- ModeSwitchCmd 中**有的字段** → 写入 XML（覆盖）
- ModeSwitchCmd 中**没有的字段** → XML 保留原值（不覆盖）
- 这样允许 XML 中保留一些 SAR ModeSwitchCmd 未覆盖的参数

**需要确认的问题**（低优先级，因为 XML fallback）：

1. **问题 A**：pulse_num 在 XML 中是否已有默认值？
   - 如有：直接保留 XML 原值即可
   - 如无：留待下一步 "参数完善" 阶段处理

2. **问题 B**：GMTI 回波文件格式是否与 SAR 兼容？
   - 需要确认：`EchoFrameHead` 是否适用于 GMTI
   - 影响：OnEchoRecvThread 中的文件读取逻辑

3. **问题 C**：SAR 端发送 `workMode == Mode_GMTI (7)` 吗？
   - 需要确认 ModeSwitchCmd 中 `workMode` 的定义
   - **当前处理**：在模式判断中添加兼容性处理

---

### Phase 1：数据结构改动（极低风险）

**文件**：`src/pipe/PipeStruDef.h`

**改动内容**：✗ **无需修改** PipeStruDef.h

- 现有 ModeSwitchCmd 已包含所有 GMTI 需要的参数
- GMTI 参数被写入 XML 后，再通过 `readXmlParam()` 读取到 Config
- 不需要定义新的结构体，避免增加复杂性

- 如果 ModeSwitchCmd 某些字段无法满足，可在 Phase 2 中定义 GMTIResultHeader
- 结果发送协议留待下一步完善


### Phase 2：主控制器 MainCtrl 改动（高风险）

#### 2.1 MainCtrl.h：增加成员

```cpp
class MainCtrl {
    bool IsResultReady;
    uint32_t m_workmode;
    char * m_SendBuf;
    FileOnlyStream file_stream;
    
    // ===== 新增成员 (GMTI 专用) =====
    Config cfg_;                           // GMTI 配置对象
    GMTIProcessor gmti_proc_;              // GMTI 处理器
    std::vector<std::vector<double>> all_frames_MT_;  // 所有帧的目标
    std::vector<int> periodList_;          // 周期列表
    std::vector<GMTIOutput> periodResults_;// 周期处理结果
    
    // 辅助函数
    bool mapCmdToConfig(const ModeSwitchCmd& cmd, Config& cfg);
    
    // 线程同步（如需要）
    std::mutex result_mutex_;              // 保护 current_tracks_
};
```

#### 2.2 修改 ProcessDataThread（关键改动）

**当前代码**（第 98-130 行 MainCtrl.cpp）：

```cpp
void* MainCtrl::ProcessDataThread(void* param) {
    // 从队列读文件名，解析 PRT 头 → TODO（未实现）
    // sleep(2);
    // IsResultReady = true;
}
```

**改为**（逻辑流程）：

```cpp
void* MainCtrl::ProcessDataThread(void* param) {
    MainCtrl *pHost = (MainCtrl*)param;
    while(true) {
        if (g_fileName.size() == 0) {
            usleep(100);
            continue;
        
        std::string fileName = g_fileName.front();
        g_fileName.pop();
        
        if (pHost->m_workmode == Mode_GMTI) {
            // ===== GMTI 模式处理 =====
            if (!pHost->gmti_proc_.readXmlParam(CONFIG_PATH, pHost->cfg_)) {
                // XML 读取失败处理
            }
            
            // 处理单个回波文件 (或周期)
            std::vector<std::vector<double>> all_frames_MT;
            if (!pHost->gmti_proc_.processPeriodsParallel(
                    pHost->periodList_, pHost->cfg_, 
                    std::vector<std::vector<double>>(), 
                    pHost->periodResults_)) {
                // 处理失败，继续
                continue;
            }
            
            // 收集 MT 结果
            for (auto& res : pHost->periodResults_) {
                if (!res.MT.empty()) {
                    all_frames_MT.insert(all_frames_MT.end(), 
                                        res.MT.begin(), res.MT.end());
                }
            }
            
            // 执行 Kalman 跟踪
            const double v_max = 100.0;
            const double sigma_thresh = 2000.0;
            if (!pHost->gmti_proc_.trackFromMT(
                    all_frames_MT, pHost->cfg_, 
                    v_max, sigma_thresh, 1, 5, 
                    pHost->current_tracks_, false)) {
                continue;
            }
            
            pHost->IsResultReady = true;
        } else {
            // ===== SAR 模式处理（保留现有逻辑）=====
            // ... 现有 SAR 处理代码 ...
        }
    }
}
```

- 测试点：两种工作模式的切换验证

#### 2.3 修改 OnResSendThread（关键改动）

**当前代码**（第 164-180 行 MainCtrl.cpp）：

```cpp
void* MainCtrl::OnResSendThread(void* param) {
    while(true) {
        if(pHost->IsResultReady == false) {
            sleep(1);
            continue;
        }
        // pHost->m_SendImagePipe.WriteData((char *)m_SendBuf, sendLen);
        pHost->IsResultReady = false;
    }
}
```

**改为**（新增 GMTI 分支）：

```cpp
void* MainCtrl::OnResSendThread(void* param) {
    MainCtrl *pHost = (MainCtrl*)param;
    
    while(true) {
        if(pHost->IsResultReady == false) {
            sleep(1);
            continue;
        }
        
        // ===== 模式分支 =====
        if (pHost->m_workmode == Mode_GMTI) {
            // ===== GMTI 结果序列化 =====
            GMTIResultHeader header;
            header.head = 0xAABB;
            header.trackCount = pHost->current_tracks_.size();
            header.tail = 0xCCDD;
            
            // 计算总长度
            uint32_t totalLen = sizeof(GMTIResultHeader) 
                              + pHost->current_tracks_.size() * sizeof(Track);
            header.msgLen = totalLen;
            
            // 序列化到缓冲区
            memcpy(pHost->m_SendBuf, &header, sizeof(GMTIResultHeader));
            uint32_t offset = sizeof(GMTIResultHeader);
            for (const auto& track : pHost->current_tracks_) {
                memcpy(pHost->m_SendBuf + offset, &track, sizeof(Track));
            // 发送
            pHost->m_SendImagePipe.WriteData(pHost->m_SendBuf, totalLen);
            
        } else {
            // ===== SAR 结果发送（保留现有逻辑）=====
            // ... 现有 SAR 发送代码 ...
        }
        
        pHost->IsResultReady = false;
    }
}
```

**改动范围**：

- 修改行数：~40 行新增
- 依赖：需确认 Track 结构体的定义和二进制兼容性
- 测试点：结果正确序列化，管道接收端能解析

---

### Phase 3：主程序适配（低风险）

**文件**：`src/pipe/main.cpp`

**当前**（只有 11 行）：

```cpp
int main() {
    MainCtrl SARCtrl;
    while(true) sleep(1);
    return 0;
}
```

**改为**（支持双模式）：

```cpp
int main(int argc, char** argv) {
    // 支持两种运行模式：
    // 1. 管道模式（推荐）：./gmti_pipe_core
    //    从管道接收参数，执行处理，发送结果
    // 2. 直接模式（调试）：./gmti_pipe_core config.xml
    //    使用 XML 配置文件直接运行，用于测试
    
    MainCtrl controller;
    
    // 如果提供了 XML 参数，先执行一次直接模式
    if (argc >= 2) {
        std::cout << "[INFO] 以 XML 模式运行: " << argv[1] << std::endl;
        // 可选：执行一次处理后退出，或启动管道等待后续命令
    }
    
    std::cout << "[INFO] GMTI 管道服务启动" << std::endl;
    
    // 主线程进入等待循环
    // MainCtrl 的四个线程独立运行
    while(true) {
        sleep(1);
    }
    
    return 0;
}
```

**改动范围**：

- 修改行数：+8 行（很少改动）
- 风险：极低，主要是添加信息日志

---

### Phase 4：CMakeLists.txt 整合（低风险）

**改动点**：添加新的可执行文件目标（不删除现有目标）

```cmake
# 现有 GMTI_core 目标保持不变
add_executable(GMTI_core 
    src/main.cpp
    # ... 其他源文件 ...
)

# ===== 新增 GMTI_pipe_core 目标 =====
add_executable(GMTI_pipe_core
    src/pipe/main.cpp
    src/pipe/MainCtrl.cpp
    src/pipe/PipeRW.cpp
    src/pipe/MyFileStream.h
    # 共享 GMTI 处理代码
    src/GMTIProcessor.cpp
    src/kalman_track.cpp
    # ... 其他共享源文件 ...
)

target_link_libraries(GMTI_pipe_core
    ${CUDA_LIBRARIES}
    ${OPENMP_LIBRARIES}
    # ... 其他依赖库 ...
```

- 新增行数：~15 行
- 影响：编译时间增加（共享代码库，编译一次）
- 验证：`make` 后生成两个可执行文件

---

## ⚡ 改动汇总表

| 模块 | 文件 | 改动类型 | 行数 | 风险 | 状态 |
|------|------|--------|------|------|------|
| **数据结构** | `src/pipe/PipeStruDef.h` | 新增 | +150 | 低 | ✓ |
| **主控头文件** | `src/pipe/MainCtrl.h` | 新增成员 | +20 | 中 | ✓ |
| **数据处理线程** | `src/pipe/MainCtrl.cpp` | 改动 | +50 修改 | 高 | ⚠️ 需测试 |
| **结果发送线程** | `src/pipe/MainCtrl.cpp` | 改动 | +40 修改 | 高 | ⚠️ 需测试 |
| **主程序** | `src/pipe/main.cpp` | 改动 | +8 | 低 | ✓ |
| **构建系统** | `CMakeLists.txt` | 新增目标 | +15 | 低 | ✓ |
| **原 GMTI 主程序** | `src/main.cpp` | **无改动** | 0 | 无 | ✓ 保留 |
| **GMTI 处理器** | `src/GMTIProcessor.*` | **无改动** | 0 | 无 | ✓ 保留 |

**总计**：
- 新增代码：~230 行
- 修改代码：~100 行（主要是条件分支）
- 删除代码：0 行
- **总改动范围：中等**

---

## 🔍 明确的边界

### ✓ 做的事

1. **参数接收**：通过 ModeSwitchCmd 接收 GMTI 参数，映射到 Config 结构
2. **算法执行**：调用现有 GMTIProcessor 执行处理和跟踪
3. **结果发送**：通过管道发送 Track 结果（二进制序列化）
4. **模式隔离**：通过 `workMode` 区分 SAR 和 GMTI，允许单一应用支持多种模式
5. **向后兼容**：保留 src/main.cpp 直接运行方式，新增 GMTI_pipe_core 程序

### ✗ 不做的事

1. **修改算法核心**：GMTIProcessor、kalman_track、pulseCompression 等不动
2. **修改 SAR 功能**：PipeRW、SAR 处理逻辑完全保留
3. **添加日志系统**：不重构日志架构，最多添加几行 DBG/ERR 输出
4. **线程池或高级并发**：保留现有 4 线程 + OpenMP 架构
5. **配置中心或参数服务**：参数仅来自 ModeSwitchCmd，不支持动态配置
6. **结果缓存或消息队列**：直接按 SAR 模式处理（一次一个结果）

### ⚠️ 需要确认的边界

| 问题 | 当前假设 | 需确认 |
|------|--------|--------|
| **pulse_num 来源** | 从 XML 读取（fallback） | SAR 端是否传递此字段？ |
| **回波文件格式** | 与 SAR 兼容（EchoFrameHead） | GMTI 文件是否同格式？ |
| **工作模式值** | `workMode == 7 (Mode_GMTI)` | SAR 端是否发送此值？ |
| **结果格式** | 二进制 GMTIResultHeader + Track | SAR 端期望何种格式？ |
| **管道路径** | 仍用 SAR 默认 `/home/raco/pipe*` | GMTI 是否需要独立管道？ |

---

## 📊 需要进一步完善的地方

### Critical（必须在后续版本实现）

1. **参数映射完整性**
   - [ ] 确定 `pulse_num` 与 `rg_len` 的确定方式
   - [ ] 验证 ModeSwitchCmd 中所有字段都有对应
   - [ ] 文件打开失败处理
   - [ ] GMTI 处理异常捕获
   - [ ] 管道断开重连机制
   - [ ] 超时控制（防止文件卡死）

3. **二进制兼容性**
   - [ ] 验证 Track 结构体在发送端和接收端的字节对齐
   - [ ] 检查浮点数精度问题
   - [ ] 添加校验和（CRC/MD5）到 GMTIResultHeader

### Important（建议在第二个版本实现）

4. **并发安全**
   - [ ] 保护 `current_tracks_` 和 `all_frames_MT_` 的互斥锁
   - [ ] 处理 ProcessDataThread 和 OnResSendThread 的竞态条件
   - [ ] 检验多个 GMTI 周期并行处理的线程安全性

5. **日志和调试**
   - [ ] 添加详细的日志（区分 GMTI 和 SAR）
   - [ ] 记录参数转换过程
   - [ ] 结果序列化的调试输出

6. **管道通信完善**
   - [ ] 管道路径由环境变量或配置指定（而非硬编码）
   - [ ] 超时处理（如 5 分钟无数据输入则超时）
   - [ ] 协议版本验证（在 ModeSwitchCmd 中增加 version 字段）

### Nice-to-Have（可选优化）

7. **性能优化**
8. **运维功能**
   - [ ] 优雅关闭机制（SIGTERM/SIGINT 处理）
   - [ ] 状态监控接口（通过额外的管道或 socket）
   - [ ] 性能统计（处理时间、目标数、跟踪率等）

9. **用户体验**
   - [ ] 启动时的参数验证和反馈
   - [ ] 处理进度实时报告
   - [ ] 结果统计摘要（目标数、跟踪数、处理时间）

---

## 📝 实施顺序建议

### 第一轮（第 1 周）：核心功能
5. ✓ 修改 src/pipe/main.cpp
6. ✓ 更新 CMakeLists.txt
7. ✓ **编译和单元测试**：确保 GMTI_pipe_core 能编译、能启动

### 第二轮（第 2-3 周）：测试和调试
1. ✓ 模拟管道测试：发送 ModeSwitchCmd，验证参数解析
2. ✓ 验证 Track 结果正确序列化
3. ✓ 集成测试：完整流程测试（参数→处理→结果）
4. ✓ 向后兼容测试：确保 GMTI_core XML 模式仍可用
5. ✓ 性能测试：多周期并行处理的性能基准

### 第三轮（第 4 周）：完善
1. ✓ 添加错误处理和日志
2. ✓ 线程安全性审查和加锁
3. ✓ 参数范围校验
4. ✓ 管道超时控制

---

## ✅ 验收标准

### 编码完成后的检查清单

- [ ] **编译**
  - [ ] `make` 无错误、无 warning
  - [ ] 生成两个可执行文件：`GMTI_core`（直接模式）、`GMTI_pipe_core`（管道模式）
  
- [ ] **功能**
  - [ ] GMTI_core 仍可用 XML 模式运行（向后兼容）
  - [ ] GMTI_pipe_core 启动后进入监听状态（4 个线程活跃）
  - [ ] 模拟发送 ModeSwitchCmd，能正确解析参数
  - [ ] 能完整执行 GMTI 处理（检测→跟踪）并发送结果
  
- [ ] **结果质量**
  - [ ] 直接模式和管道模式的处理结果一致
  - [ ] 二进制结果能被正确解析（校验头尾标志）
  - [ ] 多周期处理无数据泄漏
  
- [ ] **稳定性**
  - [ ] 异常输入（无效参数、丢失文件）不导致崩溃
  - [ ] 长期运行（>1小时）无内存泄漏
  - [ ] 管道异常断开后能恢复
  
- [ ] **性能**
  - [ ] 处理时间与直接 XML 模式相同（±5%）
  - [ ] 内存占用在预期范围内（<2GB）

---

## 🚩 风险识别与缓解

| # | 风险 | 概率 | 影响 | 缓解措施 |
|---|------|------|------|--------|
| 1 | pulse_num 来源不确定 | 中 | 高 | 先用 XML fallback，后续确认后补充 |
| 2 | 回波文件格式不一致 | 中 | 高 | 在 OnEchoRecvThread 中添加格式验证 |
| 3 | 线程竞态导致结果错误 | 中 | 高 | 添加互斥锁保护共享数据 |
| 4 | 二进制序列化字节对齐问题 | 低 | 高 | 使用 #pragma pack 确保一致，添加单元测试 |
| 5 | 现有 SAR 功能被破坏 | 低 | 中 | 完整的回归测试，明确的 workMode 分支 |
---

## 📚 参考文件清单

```
GMTI_APX15_CUDA_latest/
├── 规范.md                          # 编码规范（本计划遵照）
├── CMakeLists.txt                   # 构建系统（需更新）
│
├── src/
│   ├── main.cpp                     # GMTI 原主程序（不修改）
│   ├── GMTIProcessor.hpp/.cpp        # GMTI 处理器（不修改）
│   ├── kalman_track.cpp             # 卡尔曼跟踪（不修改）
│   └── ...
│
└── src/pipe/                        # SAR 管道通信（主要修改）
    ├── main.cpp                     # 程序入口（小改动）
```

---

## 📞 未来的对接点

1. **SAR 上位机软件**
   - 发送 ModeSwitchCmd（workMode = Mode_GMTI）
   - 发送回波数据文件名通过管道
   - 接收 GMTIResultHeader + Track 二进制数据

2. **位姿提供系统**
   - 如需实时位姿，需扩展 ModeSwitchCmd 或新增位姿管道

3. **地理信息系统（GIS）**
   - 结果中包含 planeX/Y（E, N 坐标），可直接映射

4. **性能监控系统**
   - 可选：通过文件日志或独立管道报告处理进度和性能指标

---

|  1  |  2  |  3  |  4  |
| --- | --- |---  |---  |
| 1.0 |2025-05-27 | 初版计划 | AI Assistant |

---

**下一步**：需要项目团队确认上述"需要确认的问题"（A、B、C）后，可正式启动 Phase 1 编码。

