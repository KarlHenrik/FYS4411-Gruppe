#include <iostream>
#include <cmath>
#include <vector>
#include "sampler.h"
#include "system.h"
#include "particle.h"
#include "Hamiltonians/hamiltonian.h"
#include "WaveFunctions/wavefunction.h"

#include <string>
#include <sstream>
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
        m_counter.push_back(0);
    }
}

void Sampler::setNumberOfSamples(int samples) {
    m_numberOfSamples = samples;
}

void Sampler::sample(bool acceptedStep, std::vector<Particle*> particles, int thread_num) {
    // Make sure the sampling variable(s) are initialized at the first step.
    if (m_stepNumber.at(thread_num) == 0) {
        m_cumulativeEnergy.at(thread_num) = 0;
        m_cumulativeEnergy2.at(thread_num) = 0;
    }
    // Only calculate and update values when step is accepted, and after equilibration!
    if (acceptedStep) {
        m_counter.at(thread_num)++;
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
    double averageFac = 1.0 / m_numberOfSamples;
    // Summing values from all threads
    for (int i = 0; i < m_num_threads; i++) {
        m_energy += m_cumulativeEnergy.at(i);
        m_energy2 += m_cumulativeEnergy2.at(i);
        m_counts += m_counter.at(i);
    }
    // Scaling values to be correct
    m_energy *= averageFac;
    m_energy2 *= averageFac;
}

string Sampler::outputText() {
    stringstream buffer;

    int p  = m_system->getWaveFunction()->getNumberOfParameters();
    std::vector<double> pa = m_system->getWaveFunction()->getParameters();
    double var = m_energy2 - (m_energy*m_energy);
    for (int i=0; i < p; i++) {
        buffer << left << " " << setw(15) << pa.at(i);
    }
    buffer << setw(15) << m_energy;
    buffer << setw(15) << m_energy2;
    buffer << setw(15) << var;
    buffer << setw(15) << m_counts/m_numberOfSamples;
    buffer << endl;

    return buffer.str();
}
