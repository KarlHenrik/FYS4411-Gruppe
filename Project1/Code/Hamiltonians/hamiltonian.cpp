#include "hamiltonian.h"
#include "../system.h"
#include "../particle.h"
#include "../WaveFunctions/wavefunction.h"

Hamiltonian::Hamiltonian(System* system) {
    m_system = system;
}

double Hamiltonian::computeKinetic(std::vector<Particle*> particles) {
    //return -0.5 * m_system->getWaveFunction()->ComputeLocalFullDer(particles);
    return -0.5 * m_system->getWaveFunction()->computeDoubleDerivative(particles);
    //computeLocalDoubleDerivative(particles);

    //
}
