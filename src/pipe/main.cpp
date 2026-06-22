#include <iostream>
#include <string>
#include <vector>
#include <pthread.h>
#include "pipe/MainCtrl.h"
#include <csignal>
#include <atomic>
#include <unistd.h>
#include "trig_lut.hpp"


using namespace std;

static atomic<bool> g_running{true};

void handle_signal(int)
{
    g_running.store(false);
}

namespace {

bool applyTrigArg(const string& arg, const char* next, bool& consumedNext)
{
    consumedNext = false;
    const string prefix = "--trig-mode=";
    if (arg.compare(0, prefix.size(), prefix) == 0) {
        return gmti::trig_lut::setModeFromString(arg.c_str() + prefix.size(), false);
    }
    if (arg == "--trig-mode") {
        if (!next) {
            cerr << "[ERR] --trig-mode 需要参数: lut|math|compare\n";
            return false;
        }
        consumedNext = true;
        return gmti::trig_lut::setModeFromString(next, false);
    }
    return true;
}

void printUsage(const char* prog)
{
    cerr << "Usage:\n"
         << "  " << prog << " [--trig-mode lut|math|compare]\n"
         << "  " << prog << " [--trig-mode lut|math|compare] --local-test <xml> <echo_file1> [echo_file2 ...]\n"
         << "  " << prog << " [--trig-mode lut|math|compare] --local-test <xml> <id=echo_file1> [id=echo_file2 ...]\n";
}

} // namespace

int main(int argc, char** argv)
{
    gmti::trig_lut::configureFromEnv();

    vector<string> args;
    args.reserve(static_cast<size_t>(argc > 0 ? argc - 1 : 0));
    for (int i = 1; i < argc; ++i) {
        const string arg = argv[i];
        bool consumedNext = false;
        if (arg == "--trig-mode" || arg.compare(0, string("--trig-mode=").size(), "--trig-mode=") == 0) {
            if (!applyTrigArg(arg, (i + 1 < argc) ? argv[i + 1] : nullptr, consumedNext)) {
                return 1;
            }
            if (consumedNext) ++i;
        } else {
            args.push_back(arg);
        }
    }

    if (!gmti::trig_lut::initialize(true)) {
        std::cerr << "[ERR] 三角函数路径初始化失败\n";
        return 1;
    }
    gmti::trig_lut::benchmark(1u << 22);

    if (!args.empty() && args[0] == "--local-test") {
        if (args.size() < 3) {
            printUsage(argv[0]);
            return 1;
        }

        const string xmlPath = args[1];
        vector<string> echoFiles;
        echoFiles.reserve(args.size() - 2);
        for (size_t i = 2; i < args.size(); ++i) {
            echoFiles.push_back(args[i]);
        }

        MainCtrl gmtiCtrl("", false);
        return gmtiCtrl.RunLocalTest(xmlPath, echoFiles) ? 0 : 1;
    }

    if (!args.empty()) {
        printUsage(argv[0]);
        return 1;
    }

    std::signal(SIGINT, handle_signal);
    std::signal(SIGTERM, handle_signal);

    MainCtrl gmtiCtrl;

    while (g_running.load()) {
        sleep(1);
    }

    std::cout << "GMTI exiting..." << std::endl;
    gmtiCtrl.StopThreads();
    return 0;
}
