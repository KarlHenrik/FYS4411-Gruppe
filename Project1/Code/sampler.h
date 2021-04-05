#pragma once
#include "particle.h"

#include <string>

using namespace std;

class Sampler {
public:
    Sampler(class System* system, int num_threads);
    virtual void sample(bool acceptedStep, vector<Particle*> particles, int thread_num);
    virtual void updateVals(vector<Particle*> particles, int thread_num);
    virtual string outputText();
    virtual void computeAverages();
    double getEnergy() { return m_energy_EV; }
    double getParamDer() { return m_dEnergy; }
protected:
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

class SavingSampler : public Sampler {
public:
    SavingSampler(class System* system, int num_threads, string fileName);
    void sample(bool acceptedStep, vector<Particle*> particles, int thread_num);
    void updateVals(vector<Particle*> particles, int thread_num);
    string outputText();
    void computeAverages() { return; }
    void oneBody(vector<Particle*> particles, int thread_num);
private:
    double *m_arr_energy;
    double **m_oneBody;
    int m_paralellSize;
    vector<int> m_step;
    double m_max_rad = 4;
    int m_bins = 50;
    string m_fileName;
};
