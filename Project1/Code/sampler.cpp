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

        m_energy.at(thread_num) = m_system->getHamiltonian()->computeEnergy(particles, thread_num);
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
    m_energy.at(thread_num) = m_system->getHamiltonian()->computeEnergy(particles, thread_num);
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

SavingSampler::SavingSampler(System* system, int num_threads, string fileName) : Sampler(system, num_threads) {
    m_arr_energy = new double [system->getMetroSteps()];
    m_oneBody = new double *[num_threads];
    m_fileName = fileName;
    for (int t = 0; t < num_threads; t++) {
        m_oneBody[t] = new double [m_bins];
        for (int i = 0; i < m_bins; i++) {
            m_oneBody[t][i] = 0;
        }
    }
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
        m_energy.at(thread_num) = m_system->getHamiltonian()->computeEnergy(particles, thread_num);
    }
    m_step.at(thread_num)++;
    m_arr_energy[m_paralellSize * thread_num + m_step.at(thread_num) - 1] = m_energy.at(thread_num);
    oneBody(particles, thread_num);
}

void SavingSampler::updateVals(vector<Particle*> particles, int thread_num) {
    m_energy.at(thread_num) = m_system->getHamiltonian()->computeEnergy(particles, thread_num);
}

string SavingSampler::outputText() {
    ofstream ofile;
    ofile.open("Output/" + m_fileName + "_oB" + ".txt");
    ofile << setiosflags(ios::showpoint | ios::uppercase);

    ofile << setw(15) << "r";
    ofile << setw(15) << "Density" << endl;

    double step = m_max_rad / m_bins;
    for (int i = 0; i < m_bins; i++) {
        ofile << setw(15) << (i + 1) * step;

        double density = 0;
        for (int t = 0; t < m_num_threads; t++) {
            density += m_oneBody[t][i];
        }
        ofile << setw(15) << density << endl;
    }
    ofile.close();

    // Normal output
    stringstream buffer;
    for (int i = 0; i < m_system->getMetroSteps(); i++) {
        buffer << setw(15) << m_arr_energy[i] << endl;
    }
    double accepted = 0;
    for (int t = 0; t < m_num_threads; t++) {
        accepted += m_counter.at(t);
    }
    buffer << setw(15) << accepted / m_system->getMetroSteps() << endl;

    return buffer.str();
}

void SavingSampler::oneBody(vector<Particle*> particles, int thread_num) {
    double step = m_max_rad / m_bins;
    double rad;
    for (int p1 = 0; p1 < (int) particles.size(); p1++) {
        rad = sqrt(particles[p1]->getLengthSq());

        for (int bin = 0; bin < m_bins; bin++) {
            if (rad < (bin + 1) * step) {
                m_oneBody[thread_num][bin]++;
                break;
            }
        }
    }
}
