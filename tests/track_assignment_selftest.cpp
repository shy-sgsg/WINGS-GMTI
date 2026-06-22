#include "TrackManager.hpp"

#include <iostream>

int main()
{
    return runTrackAssociationSelfTests(std::cout) ? 0 : 1;
}
