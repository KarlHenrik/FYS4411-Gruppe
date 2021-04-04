#include "harmonicoscillator.h"
#include "../system.h"
#include "../particle.h"
#include "../WaveFunctions/wavefunction.h"
#include <vector>

using namespace std;

HarmonicOscillator::HarmonicOscillator(System* system, double omega, double beta) : Hamiltonian(system) {
    m_omega  = omega;
    m_omegaSq = omega * omega;
    ell = {1, 1, beta};
}

double HarmonicOscillator::computeEnergy(std::vector<Particle*> particles, int thread) {
    double potentialEnergy = 0;
    for (int i = 0; i < (int) particles.size(); i++) {
        for (int d = 0; d < particles[0]->getDims(); d++) {
            potentialEnergy += 0.5 * m_omegaSq * particles[i]->getPosition().at(d) * particles[i]->getPosition().at(d) * ell.at(d);
        }
    }
    double kineticEnergy = computeKinetic(particles, thread);
    return kineticEnergy + potentialEnergy;
}
