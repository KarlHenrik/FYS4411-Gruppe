#pragma once
#include "hamiltonian.h"
#include <vector>

class HarmonicOscillator : public Hamiltonian {
public:
    HarmonicOscillator(System* system, double omega);
    double computeEnergy(std::vector<Particle*> particles, int thread);

private:
    double m_omega = 0;
    double m_omegaSq = 0;
};
