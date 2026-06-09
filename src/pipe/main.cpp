#include <iostream>
#include <string>
#include <vector>
#include <pthread.h>
#include "pipe/MainCtrl.h"

using namespace std;

namespace {

void printUsage(const char* prog)
{
    cerr << "Usage:\n"
         << "  " << prog << "\n"
         << "  " << prog << " --local-test <xml> <echo_file1> [echo_file2 ...]\n"
         << "  " << prog << " --local-test <xml> <id=echo_file1> [id=echo_file2 ...]\n";
}

} // namespace

int main(int argc, char** argv)
{
    if (argc >= 2 && string(argv[1]) == "--local-test") {
        if (argc < 4) {
            printUsage(argv[0]);
            return 1;
        }

        const string xmlPath = argv[2];
        vector<string> echoFiles;
        echoFiles.reserve(static_cast<size_t>(argc - 3));
        for (int i = 3; i < argc; ++i) {
            echoFiles.push_back(argv[i]);
        }

        MainCtrl gmtiCtrl("", false);
        return gmtiCtrl.RunLocalTest(xmlPath, echoFiles) ? 0 : 1;
    }

    if (argc != 1) {
        printUsage(argv[0]);
        return 1;
    }

    GMTIProcessor proc;
    Config cfg;
    if (!proc.readXmlParam("temp_config.xml", cfg)) {
        cfg.pipe_root_path = "/home/raco";
    }

    MainCtrl gmtiCtrl(cfg.pipe_root_path);

    pthread_exit(nullptr);
}
