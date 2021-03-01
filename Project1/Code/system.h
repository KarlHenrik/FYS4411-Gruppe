#pragma once
#include <vector>
#include <Math/random.h>
#include "particle.h"
#include <string>
#include <fstream>
#include <iomanip>

using namespace std;

class System {
public:
    System();
    System(int num_threads);
    System(int num_threads, int seed);
    bool ImmetropolisStep           (std::vector<Particle*> particles, double& waveFuncValue);
    bool metropolisStep             (std::vector<Particle*> particles, double& waveFuncValue);
    void runMetropolisSteps         (bool m_choice);
    void setNumberOfParticles       (int numberOfParticles);
    void setNumberOfDimensions      (int numberOfDimensions);
    void setNumberOfSteps           (int numberOfSteps);
    void setStepLength              (double stepLength);
    void setTimeStep                (double timestep);
    void setEquilibrationFraction   (double equilibrationFraction);
    void setChoice                  (bool ImpSampling);
    void setHamiltonian             (class Hamiltonian* hamiltonian);
    void setWaveFunction            (class WaveFunction* waveFunction);
    void setInitialState            (class InitialState* initialState);
    void addOutput                  (string text);
    class WaveFunction*             getWaveFunction()   { return m_waveFunction; }
    class Hamiltonian*              getHamiltonian()    { return m_hamiltonian; }
    class Sampler*                  getSampler()        { return m_sampler; }
    class Random*                   getRandomEngine();
    class InitialState*             getInitialState()   { return m_initialState; }
    int getNumberOfParticles()          { return m_numberOfParticles; }
    int getNumberOfDimensions()         { return m_numberOfDimensions; }
    int getNumberOfMetropolisSteps()    { return m_numberOfMetropolisSteps; }
    double getEquilibrationFraction()   { return m_equilibrationFraction; }
    string getOutput()                  { return m_output; }
    bool getChoice()                    { return m_choice; }

private:
    int                             m_num_threads;
    int                             m_numberOfParticles = 0;
    int                             m_numberOfDimensions = 0;
    int                             m_numberOfMetropolisSteps = 0;
    double                          m_equilibrationFraction = 0.0;
    double                          m_stepLength = 0.1;
    double                          m_timestep = 0.1;
    bool                            m_choice;
    class WaveFunction*             m_waveFunction = nullptr;
    class Hamiltonian*              m_hamiltonian = nullptr;
    class InitialState*             m_initialState = nullptr;
    class Sampler*                  m_sampler = nullptr;
    vector<class Random*>           m_randoms;
    string                          m_output = "";
};
