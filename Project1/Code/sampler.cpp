#include "sampler.h"
#include "system.h"
#include "particle.h"
#include "Hamiltonians/hamiltonian.h"
#include "WaveFunctions/wavefunction.h"

#include <cmath>
#include <vector>

#include <string>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <iostream>

using namespace std;


Sampler::Sampler(System* system, int num_threads) {
    m_system = system;
    m_num_threads = num_threads;

    // Setting up vectors with values for all threads
    for (int i = 0; i < num_threads; i++) {
        m_counter.push_back(0);

        m_energy.push_back(0);
        m_energy2.push_back(0);
        m_dPsi.push_back(0);

        m_total_energy.push_back(0);
        m_total_energy2.push_back(0);
        m_total_dPsi.push_back(0);
        m_total_dPsiEL.push_back(0);
    }
}

void Sampler::sample(bool acceptedStep, std::vector<Particle*> particles, int thread_num) {
    // Only calculate and update values when step is accepted, and after equilibration!
    if (acceptedStep) {
        m_counter.at(thread_num)++;

        m_energy.at(thread_num) = m_system->getHamiltonian()->computeEnergy(particles);
        m_energy2.at(thread_num) = m_energy.at(thread_num) * m_energy.at(thread_num);
        m_dPsi.at(thread_num) = m_system->getWaveFunction()->computeParamDer(particles);
    }
    m_total_energy.at(thread_num) += m_energy.at(thread_num);
    m_total_energy2.at(thread_num) += m_energy2.at(thread_num);
    m_total_dPsi.at(thread_num) += m_dPsi.at(thread_num);
    m_total_dPsiEL.at(thread_num) += m_dPsi.at(thread_num) * m_energy.at(thread_num);
}

void Sampler::updateVals(std::vector<Particle*> particles, int thread_num) {
    // Updates values after equilibration
    m_energy.at(thread_num) = m_system->getHamiltonian()->computeEnergy(particles);
    m_energy2.at(thread_num) = m_energy.at(thread_num) * m_energy.at(thread_num);
    m_dPsi.at(thread_num) = m_system->getWaveFunction()->computeParamDer(particles);
}

void Sampler::computeAverages() {
    // Not all steps are sampled, and we need to divide my the number of threads
    double averageFac = 1.0 / m_system->getMetroSteps();
    // Summing values from all threads
    for (int i = 0; i < m_num_threads; i++) {
        m_counts += m_counter.at(i);

        m_energy_EV += m_total_energy.at(i);
        m_energy2_EV += m_total_energy2.at(i);
        m_dPsi_EV += m_total_dPsi.at(i);
        m_dPsiEL_EV += m_total_dPsiEL.at(i);
    }
    // Scaling values to be correct
    m_energy_EV *= averageFac;
    m_energy2_EV *= averageFac;
    m_dPsi_EV *= averageFac;
    m_dPsiEL_EV *= averageFac;
    // Calculating derivative of wave function wrt. variational parameter
    m_dEnergy = 2 * (m_dPsiEL_EV - m_dPsi_EV * m_energy_EV);
}

string Sampler::outputText() {
    stringstream buffer;

    int p  = m_system->getWaveFunction()->getNumberOfParameters();
    std::vector<double> pa = m_system->getWaveFunction()->getParameters();
    double var = m_energy2_EV - (m_energy_EV * m_energy_EV);
    for (int i=0; i < p; i++) {
        buffer << left << " " << setw(15) << pa.at(i);
    }
    buffer << setw(15) << m_energy_EV;
    buffer << setw(15) << m_energy2_EV;
    buffer << setw(15) << var;
    buffer << setw(15) << m_counts / m_system->getMetroSteps();
    buffer << endl;

    return buffer.str();
}

// ------------- SavingSampler functions -----------------

SavingSampler::SavingSampler(System* system, int num_threads) : Sampler(system, num_threads) {
    m_arr_energy = new double[system->getMetroSteps()];
    m_paralellSize = (int) system->getMetroSteps() / num_threads;

    // Setting up vectors with values for all threads
    for (int i = 0; i < num_threads; i++) {
        m_step.push_back(0);
    }
}

void SavingSampler::sample(bool acceptedStep, vector<Particle*> particles, int thread_num) {
    // Only calculate and update values when step is accepted, and after equilibration!
    if (acceptedStep) {
        m_counter.at(thread_num)++;
        m_energy.at(thread_num) = m_system->getHamiltonian()->computeEnergy(particles);
    }
    m_total_energy.at(thread_num) += m_energy.at(thread_num);

    m_step.at(thread_num)++;
    m_arr_energy[m_paralellSize * thread_num + m_step.at(thread_num) - 1] = m_total_energy.at(thread_num) / m_step.at(thread_num);
}

void SavingSampler::updateVals(vector<Particle*> particles, int thread_num) {
    m_energy.at(thread_num) = m_system->getHamiltonian()->computeEnergy(particles);
}

string SavingSampler::outputText() {
    stringstream buffer;
    for (int i = 0; i < m_system->getMetroSteps(); i++) {
        buffer << setw(15) << m_arr_energy[i] << endl;
    }
    return buffer.str();
}
void SavingSampler::computeAverages() {
    return;
}
