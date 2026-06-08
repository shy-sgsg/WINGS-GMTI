#ifndef NEW_PROTOCOL_READER_HPP
#define NEW_PROTOCOL_READER_HPP

#include "config_structs.hpp"

#include <complex>
#include <vector>

bool readPulseBlockNewProtocol(const Config &cfg,
                               int beamskip,
                               std::vector<std::complex<float>> &data1,
                               std::vector<std::complex<float>> &data2,
                               std::vector<double> &utc,
                               double &theta_sq,
                               std::vector<std::vector<double>> &posRaw);

#endif // NEW_PROTOCOL_READER_HPP
