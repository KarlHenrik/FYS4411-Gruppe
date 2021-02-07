#include "wavefunction.h"

WaveFunction::WaveFunction(System* system) {
    m_system = system;
}

void WaveFunction::setParameters(std::vector<double> parameters) {
    m_parameters = parameters;
}
