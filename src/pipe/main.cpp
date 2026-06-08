#include <iostream>
#include <string>
#include <pthread.h>
#include "pipe/MainCtrl.h"

using namespace std;

int main()
{
    GMTIProcessor proc;
    Config cfg;
    if (!proc.readXmlParam("temp_config.xml", cfg)) {
        cfg.pipe_root_path = "/home/raco";
    }

    MainCtrl gmtiCtrl(cfg.pipe_root_path);

    pthread_exit(nullptr);
}
