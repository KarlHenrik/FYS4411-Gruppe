#pragma once
#include "particle.h"

using namespace std;

class Sampler {
public:
    Sampler(class System* system, int num_threads);
    void setNumberOfMetropolisSteps(int steps);
    void sample(bool acceptedStep, std::vector<Particle*> particles, int thread_num);
    void updateVals(std::vector<Particle*> particles, int thread_num);
    void printOutputToTerminal();
    void computeAverages();
    double getEnergy()          { return m_energy; }

private:
    int     m_numberOfMetropolisSteps = 0;
    int     m_num_threads;
    class System* m_system = nullptr;
    // Calculated each loop
    vector<double> m_localEnergy;
    vector<double> m_localEnergy2;
    // Cumulative
    vector<double> m_cumulativeEnergy;
    vector<double> m_cumulativeEnergy2;
    vector<int> m_stepNumber;
    // End values
    double  m_energy = 0;
    double  m_energy2 = 0;
};
