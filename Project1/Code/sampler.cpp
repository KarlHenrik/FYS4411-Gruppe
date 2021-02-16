#include <iostream>
#include <cmath>
#include <vector>
#include "sampler.h"
#include "system.h"
#include "particle.h"
#include "Hamiltonians/hamiltonian.h"
#include "WaveFunctions/wavefunction.h"

#include <string>
#include <fstream>
#include <iomanip>
#include <cstdlib>

using namespace std;


Sampler::Sampler(System* system, int num_threads) {
    m_system = system;
    m_num_threads = num_threads;

    // Setting up vectors with values for all threads
    for (int i = 0; i < num_threads; i++) {
        m_localEnergy.push_back(0);
        m_localEnergy2.push_back(0);
        m_cumulativeEnergy.push_back(0);
        m_cumulativeEnergy2.push_back(0);
        m_stepNumber.push_back(0);
    }
}

void Sampler::setNumberOfMetropolisSteps(int steps) {
    m_numberOfMetropolisSteps = steps;
}

void Sampler::sample(bool acceptedStep, std::vector<Particle*> particles, int thread_num) {
    // Make sure the sampling variable(s) are initialized at the first step.
    if (m_stepNumber.at(thread_num) == 0) {
        m_cumulativeEnergy.at(thread_num) = 0;
        m_cumulativeEnergy2.at(thread_num) = 0;
    }
    // Only calculate and update values when step is accepted, and after equilibration!
    if (acceptedStep) {

        m_localEnergy.at(thread_num) = m_system->getHamiltonian()->computeLocalEnergy(particles);
        m_localEnergy2.at(thread_num) = m_localEnergy.at(thread_num)*m_localEnergy.at(thread_num);
    }
    m_cumulativeEnergy.at(thread_num) += m_localEnergy.at(thread_num);
    m_cumulativeEnergy2.at(thread_num) += m_localEnergy2.at(thread_num);
    m_stepNumber.at(thread_num)++;
}

void Sampler::updateVals(std::vector<Particle*> particles, int thread_num) {
    // Updates values after equilibration
    m_localEnergy.at(thread_num) = m_system->getHamiltonian()->computeLocalEnergy(particles);
    m_localEnergy2.at(thread_num) = m_localEnergy.at(thread_num)*m_localEnergy.at(thread_num);
}

void Sampler::computeAverages() {
    // Not all steps are sampled, and we need to divide my the number of threads
    double averageFac = 1.0 / ( m_system->getNumberOfMetropolisSteps() * (1.0 - m_system->getEquilibrationFraction()) * m_num_threads);
    // Summing values from all threads
    for (int i = 0; i < m_num_threads; i++) {
        m_energy += m_cumulativeEnergy.at(i);
        m_energy2 += m_cumulativeEnergy2.at(i);
    }
    // Scaling values to be correct
    m_energy *= averageFac;
    m_energy2 *= averageFac;
}

void Sampler::printOutputToTerminal() {
    int p  = m_system->getWaveFunction()->getNumberOfParameters();
    std::vector<double> pa = m_system->getWaveFunction()->getParameters();

    for (int i=0; i < p; i++) {
        cout << left << " " << setw(15) << pa.at(i);
    }
    cout << setw(15) << m_energy;
    cout << setw(15) << m_energy2;
    cout << setw(15) << m_energy2 - (m_energy*m_energy);
    cout << endl;
}

/*
void Sampler::printOutputToTerminal() {
    int     np = m_system->getNumberOfParticles();
    int     nd = m_system->getNumberOfDimensions();
    int     ms = m_system->getNumberOfMetropolisSteps();
    int     p  = m_system->getWaveFunction()->getNumberOfParameters();
    double  ef = m_system->getEquilibrationFraction();
    std::vector<double> pa = m_system->getWaveFunction()->getParameters();

    cout << endl;
    cout << "  -- System info -- " << endl;
    cout << " Number of particles  : " << np << endl;
    cout << " Number of dimensions : " << nd << endl;
    cout << " Number of Metropolis steps run : 10^" << std::log10(ms) << endl;
    cout << " Number of equilibration steps  : 10^" << std::log10(std::round(ms*ef)) << endl;
    cout << endl;
    cout << "  -- Wave function parameters -- " << endl;
    cout << " Number of parameters : " << p << endl;
    for (int i=0; i < p; i++) {
        cout << " Parameter " << i+1 << " : " << pa.at(i) << endl;
    }
    cout << endl;
    cout << "  -- Results -- " << endl;
    cout << " Energy : " << m_energy << endl;
    cout << endl;
}
*/
