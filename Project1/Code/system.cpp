#include "system.h"
#include <cassert>
#include "sampler.h"
#include "particle.h"
#include "WaveFunctions/wavefunction.h"
#include "Hamiltonians/hamiltonian.h"
#include "InitialStates/initialstate.h"
#include "Math/random.h"
#include "omp.h"
#include <string>
#include <iostream>
using namespace std;

System::System() {
    m_num_threads = omp_get_max_threads();

    for (int i = 0; i < m_num_threads; i++) {
        m_randoms.push_back(new Random());
    }
}

System::System(int num_threads) {
    m_num_threads = num_threads;
    omp_set_num_threads(num_threads);

    for (int i = 0; i < m_num_threads; i++) {
        m_randoms.push_back(new Random());
    }
}

System::System(int num_threads, int seed) {
    m_num_threads = num_threads;
    omp_set_num_threads(num_threads);

    for (int i = 0; i < m_num_threads; i++) {
        m_randoms.push_back(new Random(seed + i));
    }
}

bool System::metropolisStep(std::vector<Particle*> particles, double& waveFuncValue) {
    Random* private_random = m_randoms.at(omp_get_thread_num());

    Particle* randParticle = particles[private_random->nextInt(0, m_numberOfParticles - 1)];
    std::vector<double> oldPos = randParticle->getPosition();
    double oldLengthSq = randParticle->getLengthSq();

    //double dir = (private_random->nextInt(1) - 0.5) * 2;

    // Obtaining the vector contaning drift force entries from SimpleGaussian
    std::vector<double> oldQuantumForce = m_waveFunction->ComputeQF(randParticle,oldPos);
    // Declaring a vector to store the change of position in each spatial direction
    std::vector<double> move(m_numberOfDimensions, 0);

    // Looping over each spatial component to update the positions
    for (int i = 0; i < m_numberOfDimensions; i++) {
      move[i] = oldPos[i] + 0.5*m_timestep*oldQuantumForce[i]
              + private_random->nextGaussian(0.,1.)*sqrt(m_timestep);
    }

    // Changing each component of randParticles position
    randParticle->adjustLangevin(move);
    std::vector<double> newPos = randParticle->getPosition();

    //randParticle->adjustPosition(m_stepLength * dir, private_random->nextInt(0, m_numberOfDimensions - 1));
    double newWaveFuncValue = m_waveFunction->evaluateChange(randParticle, waveFuncValue, oldLengthSq);
    std::vector<double> newQuantumForce = m_waveFunction->ComputeQF(randParticle, newPos);
    double GreensFunctionRatio = 0.0;

    for (int i = 0; i < m_numberOfDimensions; i++) {
      GreensFunctionRatio += (oldQuantumForce[i] + newQuantumForce[i])
      *(oldPos[i] - newPos[i] + m_timestep*(oldQuantumForce[i] - newQuantumForce[i]));
    }

    GreensFunctionRatio = exp(0.5*GreensFunctionRatio);

    if (private_random->nextDouble() < GreensFunctionRatio*std::pow(newWaveFuncValue / waveFuncValue, 2)) {
        waveFuncValue = newWaveFuncValue;
        oldQuantumForce = newQuantumForce;
        return true;
    } else {
        randParticle->setPosition(oldPos);
        return false;
    }
}

void System::runMetropolisSteps() {
    int equilibrationSteps = m_numberOfMetropolisSteps * m_equilibrationFraction;
    int parallellSteps = (m_numberOfMetropolisSteps - equilibrationSteps) / m_num_threads + equilibrationSteps;

    m_sampler = new Sampler(this, m_num_threads);
    m_sampler->setNumberOfSamples((parallellSteps - equilibrationSteps) * m_num_threads);
    # pragma omp parallel
    {
        // Private variables setup for each thread
        std::vector<Particle*> private_particles = m_initialState->newParticles();
        double private_waveFuncValue = m_waveFunction->evaluate(private_particles);
        int thread_num = omp_get_thread_num();
        // Equilibration
        for (int i = 0; i < equilibrationSteps; i++) {
            metropolisStep(private_particles, private_waveFuncValue);
        }
        m_sampler->updateVals(private_particles, thread_num);
        // Steps with sampling
        for (int i = equilibrationSteps; i < parallellSteps; i++) {
            bool acceptedStep = metropolisStep(private_particles, private_waveFuncValue);
            m_sampler->sample(acceptedStep, private_particles, thread_num);
        }
    }
    m_sampler->computeAverages();
    addOutput(m_sampler->outputText());
}

void System::addOutput(string text) {
    m_output += text;
}

void System::setNumberOfSteps(int numberOfSteps) {
    m_numberOfMetropolisSteps = numberOfSteps;
}

void System::setNumberOfParticles(int numberOfParticles) {
    m_numberOfParticles = numberOfParticles;
}

void System::setNumberOfDimensions(int numberOfDimensions) {
    m_numberOfDimensions = numberOfDimensions;
}

void System::setStepLength(double stepLength) {
    assert(stepLength >= 0);
    m_stepLength = stepLength;
}

void System::setTimeStep(double timestep) {
  m_timestep = timestep;
}

void System::setEquilibrationFraction(double equilibrationFraction) {
    assert(equilibrationFraction >= 0);
    m_equilibrationFraction = equilibrationFraction;
}

void System::setHamiltonian(Hamiltonian* hamiltonian) {
    m_hamiltonian = hamiltonian;
}

void System::setWaveFunction(WaveFunction* waveFunction) {
    m_waveFunction = waveFunction;
}

void System::setInitialState(InitialState* initialState) {
    m_initialState = initialState;
}

Random* System::getRandomEngine() {
    return m_randoms.at(omp_get_thread_num());
}
