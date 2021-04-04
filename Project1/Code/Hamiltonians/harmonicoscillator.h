#pragma once
#include "hamiltonian.h"

using namespace std;

class HarmonicOscillator : public Hamiltonian {
public:
    HarmonicOscillator(System* system, double omega, double beta);
    double computeEnergy(std::vector<Particle*> particles, int thread);

private:
    vector<double> ell;
    double m_omega = 0;
    double m_omegaSq = 0;
};
