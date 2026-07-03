#ifndef RUNTIME_DIAGNOSTICS_HPP
#define RUNTIME_DIAGNOSTICS_HPP

#include <chrono>
#include <string>

struct Config;

namespace gmti {
namespace runtime {

struct RunPaths {
    std::string runtime_config_json;
    std::string runtime_config_txt;
    std::string timing_metrics_csv;
    std::string run_manifest_json;
};

void initializeRun(const Config& cfg,
                   const std::string& config_path,
                   const std::string& executable_path);

void finishRun(const Config& cfg, bool normal_exit, int exit_code,
               const std::string& notes = "");

void writeRuntimeConfigDump(const Config& cfg,
                            const std::string& config_path,
                            const std::string& executable_path);

void recordTiming(const char* scope_name,
                  std::chrono::system_clock::time_point start_time,
                  std::chrono::system_clock::time_point end_time,
                  long long elapsed_ms,
                  int period_id = -1,
                  const std::string& extra = "");

void flushTimingMetrics();

bool diagnosticsEnabled();
const std::string& runId();
const std::string& caseId();
const std::string& resultId();
RunPaths paths();
std::string currentIsoTime();

class TimingScope {
public:
    explicit TimingScope(const char* name, int period_id = -1,
                         const std::string& extra = "");
    ~TimingScope();

private:
    const char* name_;
    int period_id_;
    std::string extra_;
    std::chrono::high_resolution_clock::time_point steady_start_;
    std::chrono::system_clock::time_point wall_start_;
};

} // namespace runtime
} // namespace gmti

#endif // RUNTIME_DIAGNOSTICS_HPP
