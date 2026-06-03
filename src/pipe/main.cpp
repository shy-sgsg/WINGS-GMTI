#include <iostream>
#include <string>
#include <pthread.h>
#include "MainCtrl.h"
#include "GMTIProcessor.hpp"

using namespace std;

int main()
{
    Config cfg;
    GMTIProcessor proc;
    if (!proc.readXmlParam("../temp_config.xml", cfg)) {
        cfg.pipe_root_path = "/home/shy/pipe_test";
    }

    MainCtrl gmtiCtrl(cfg.pipe_root_path);

    pthread_exit(nullptr);
}
