#pragma once
#include "wavefunction.h"

using namespace std;

class Correlated : public WaveFunction {
public:
    Correlated(class System* system, double alpha, double a, double beta);
    double computeDoubleDerivative(vector<class Particle*> particles, int thread);
    vector<double> computeQF(vector<Particle*> particles, int particle_idx, Particle* randParticle, vector<double> oldPos, int thread);
    double computeRatio(vector<Particle*> particles, int particle_idx, Particle* randParticle, vector<double> oldPos, int thread_num);
    double computeParamDer(vector<Particle*> particles);
    // Only implemented in this class
    void setup(vector<Particle*> particles, int thread_num);
    void clear(vector<Particle*> particles, int thread);
    void revertDists(vector<Particle*> particles, int particle_idx, int thread);
    // Only declare in this class
    void computeDists(vector<Particle*> particles, int particle_idx, Particle* randParticle, int thread);

private:
    vector<double> ell;
    double m_a;
    double*** dists;
    double*** old_dists;
    double**** unit_vecs;
    double**** old_unit_vecs;
};
