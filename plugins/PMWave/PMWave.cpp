// PluginPMWave.cpp
// Aleksandar Koruga (aleksandar.koruga@gmail.com)

#include "SC_PlugIn.hpp"
#include "PMWave.hpp"

static InterfaceTable* ft;

namespace PMWave {

PMWave::PMWave() {
    mCalcFunc = make_calc_function<PMWave, &PMWave::next>();
    next(1);
}

void PMWave::next(int nSamples) {
    const float* input = in(0);
    const float* gain = in(1);
    float* outbuf = out(0);

    // simple gain function
    for (int i = 0; i < nSamples; ++i) {
        outbuf[i] = input[i] * gain[i];
    }
}

} // namespace PMWave

PluginLoad(PMWaveUGens) {
    // Plugin magic
    ft = inTable;
    registerUnit<PMWave::PMWave>(ft, "PMWave", false);
}
