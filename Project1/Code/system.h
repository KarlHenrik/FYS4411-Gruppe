#pragma once
#include "particle.h"
#include <string>

using namespace std;

class System {
public:
    System(int num_threads = -1, int seed = -1);
    bool imMetropolisStep           (vector<Particle*> particles);
    bool metropolisStep             (vector<Particle*> particles);
    void runMetropolisSteps         (bool m_choice = false);
    void setNumberOfParticles       (int numberOfParticles);
    void setNumberOfDimensions      (int numberOfDimensions);
    void setNumberOfSteps           (int metroSteps, int equiSteps);
    void setStepLength              (double stepLength);
    void setTimeStep                (double timestep);
    void setChoice                  (bool ImpSampling);
    void setHamiltonian             (class Hamiltonian* hamiltonian);
    void setWaveFunction            (class WaveFunction* waveFunction);
    void setInitialState            (class InitialState* initialState);
    void setSampler                 (class Sampler* sampler);
    void addOutput                  (string text);
    void clearOutput                ();
    class WaveFunction*             getWaveFunction()   { return m_waveFunction; }
    class Hamiltonian*              getHamiltonian()    { return m_hamiltonian; }
    class Sampler*                  getSampler()        { return m_sampler; }
    class Random*                   getRandomEngine();
    class InitialState*             getInitialState()   { return m_initialState; }
    int getNumberOfParticles()          { return m_numberOfParticles; }
    int getNumberOfDimensions()         { return m_numberOfDimensions; }
    int getMetroSteps()                 { return m_metroSteps; }
    int getEquiSteps()                  { return m_equiSteps; }
    string getOutput()                  { return m_output; }
    bool getChoice()                    { return m_choice; }
    int getNumerOfThreads()             { return m_num_threads; }

private:
    int                             m_num_threads;
    int                             m_numberOfParticles = 0;
    int                             m_numberOfDimensions = 0;
    int                             m_metroSteps = 0;
    int                             m_equiSteps  = 0;
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
