#pragma once
#include "particle.h"

#include <string>

using namespace std;

class Sampler {
public:
    Sampler(class System* system, int num_threads);
    void sample(bool acceptedStep, vector<Particle*> particles, int thread_num);
    void updateVals(vector<Particle*> particles, int thread_num);
    string outputText();
    void computeAverages();
    double getEnergy() { return m_energy_EV; }
    double getParamDer() { return m_dEnergy; }
private:
    int m_num_threads;
    class System* m_system = nullptr;
    // Calculated each loop
    vector<double> m_energy;
    vector<double> m_energy2;
    vector<double> m_dPsi;
    // Totals
    vector<int> m_counter;
    vector<double> m_total_energy;
    vector<double> m_total_energy2;
    vector<double> m_total_dPsi;
    vector<double> m_total_dPsiEL;
    // End values
    double m_counts = 0;
    double m_energy_EV = 0;
    double m_energy2_EV = 0;
    double m_dPsi_EV = 0;
    double m_dPsiEL_EV = 0;
    double m_dEnergy = 0;

};
/*
class SavingSampler : public Sampler {
public:
    void sample(bool acceptedStep, vector<Particle*> particles, int thread_num);
    void updateVals(vector<Particle*> particles, int thread_num);
    string outputText();
    void computeAverages();
private:

};
*/
