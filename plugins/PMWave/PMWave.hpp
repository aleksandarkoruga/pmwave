// PluginPMWave.hpp
// Aleksandar Koruga (aleksandar.koruga@gmail.com)

#pragma once

#include "SC_PlugIn.hpp"

namespace PMWave {

class PMWave : public SCUnit {
public:
    PMWave();

    // Destructor
    // ~PMWave();

private:
    // Calc function
    void next(int nSamples);

    // Member variables
};

} // namespace PMWave
