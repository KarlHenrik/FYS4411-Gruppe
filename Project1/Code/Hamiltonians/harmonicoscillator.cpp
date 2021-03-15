#include "harmonicoscillator.h"
#include "../system.h"
#include "../particle.h"
#include "../WaveFunctions/wavefunction.h"


HarmonicOscillator::HarmonicOscillator(System* system, double omega) : Hamiltonian(system) {
    m_omega  = omega;
    m_omegaSq = omega * omega;
}

double HarmonicOscillator::computeEnergy(std::vector<Particle*> particles) {
    double potentialEnergy = 0;
    for (int i = 0; i < (int) particles.size(); i++) {
        potentialEnergy += 0.5 * m_omegaSq * particles[i]->getLengthSq();
    }
    double kineticEnergy = computeKinetic(particles);
    return kineticEnergy + potentialEnergy;
}
