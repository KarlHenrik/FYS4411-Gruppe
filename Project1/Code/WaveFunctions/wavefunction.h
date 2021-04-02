#pragma once
#include <vector>

using namespace std;

class WaveFunction {
public:
    WaveFunction(class System* system);
    int     getNumberOfParameters() { return m_numberOfParameters; }
    vector<double> getParameters() { return m_parameters; }
    virtual void setParameters(vector<double>);

    virtual double computeDoubleDerivative(vector<class Particle*> particles, int thread) = 0;
    virtual vector<double> computeQF(vector<Particle*> particles, int particle_idx, Particle* randParticle, vector<double> oldPos, int thread) = 0;
    virtual double computeRatio(vector<Particle*> particles, int particle_idx, Particle* randParticle, vector<double> oldPos, int thread_num) = 0;
    virtual double computeParamDer(vector<Particle*> particles) = 0;

    virtual void clear(vector<Particle*> /*particles*/, int /*thread*/) { return; }
    virtual void setup(vector<Particle*> /*particles*/, int /*thread*/) { return; }
    virtual void revertDists(vector<Particle*> /*particles*/, int /*particle_idx*/, int /*thread*/) { return; }

protected:
    int     m_numberOfParameters = 0;
    vector<double> m_parameters = vector<double>();
    class System* m_system = nullptr;
};
