#pragma once
#include "wavefunction.h"

using namespace std;

class Correlated : public WaveFunction {
public:
    Correlated(class System* system, double alpha);
    double computeDoubleDerivative(vector<class Particle*> particles);
    vector<double> computeQF(vector<Particle*> particles, int particle_idx, Particle* randParticle, vector<double> oldPos, int thread);
    double computeRatio(vector<Particle*> particles, int particle_idx, Particle* randParticle, vector<double> oldPos, double oldLengthSq, int thread_num);
    double computeParamDer(vector<Particle*> particles);
    void setup(vector<Particle*> particles, int thread_num);
    void updateNewDist(vector<Particle*> particles, int particle_idx, Particle* randParticle, vector<double> oldPos, int thread_num);

private:

    double m_a = 1; // what should this be!!??!?
    double*** old_Dist;
    double*** new_Dist;
    double**** unit_Vec;
};
