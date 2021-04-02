#pragma once
#include "wavefunction.h"

using namespace std;

class SimpleGaussian : public WaveFunction {
public:
    SimpleGaussian(class System* system, double alpha, double a, double beta);
    double computeDoubleDerivative(vector<class Particle*> particles, int thread);
    vector<double> computeQF(vector<Particle*> particles, int particle_idx, Particle* randParticle, vector<double> oldPos, int thread);
    double computeRatio(vector<Particle*> particles, int particle_idx, Particle* randParticle, vector<double> oldPos, int thread_num);
    double computeParamDer(vector<Particle*> particles);
private:
    vector<double> ell;
};
