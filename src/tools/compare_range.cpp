#include <iostream>
#include <vector>
#include <complex>
#include "GMTIProcessor.hpp"

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: compare_range <xml-config-path>" << std::endl;
        return 1;
    }
    std::string xml = argv[1];
    Config cfg;
    GMTIProcessor proc;
    if (!proc.readXmlParam(xml, cfg)) {
        std::cerr << "Failed to read XML" << std::endl;
        return 1;
    }

    // delegate to member debug helper which has access to private methods
    if (!proc.debug_compare_range(cfg, 1)) {
        std::cerr << "debug_compare_range failed" << std::endl;
        return 1;
    }
    return 0;
}
