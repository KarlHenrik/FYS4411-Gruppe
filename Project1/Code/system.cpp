#include "system.h"
#include <cassert>
#include "sampler.h"
#include "particle.h"
#include "WaveFunctions/wavefunction.h"
#include "Hamiltonians/hamiltonian.h"
#include "InitialStates/initialstate.h"
#include "Math/random.h"
#include "omp.h"

#include <iostream>
using namespace std;

System::System() {
    m_random = new Random();
}

System::System(int seed) {
    m_random = new Random(seed);
}

bool System::metropolisStep(std::vector<Particle*> particles, double& waveFuncValue) {
    Particle* randParticle = particles[m_random->nextInt(0, m_numberOfParticles - 1)];
    std::vector<double> oldPos = randParticle->getPosition();
    double oldLengthSq = randParticle->getLengthSq();

    double dir = (m_random->nextInt(1) - 0.5) * 2;

    randParticle->adjustPosition(m_stepLength * dir, m_random->nextInt(0, m_numberOfDimensions - 1));
    double newWaveFuncValue = m_waveFunction->evaluateChange(randParticle, waveFuncValue, oldLengthSq);

    if (m_random->nextDouble() < std::pow(newWaveFuncValue / waveFuncValue, 2)) {
        waveFuncValue = newWaveFuncValue;
        return true;
    } else {
        randParticle->setPosition(oldPos);
        return false;
    }
}

void System::runMetropolisSteps() {
    int num_threads = 6;
    omp_set_num_threads(num_threads);

    m_sampler = new Sampler(this, num_threads);
    m_sampler->setNumberOfMetropolisSteps(m_numberOfMetropolisSteps);
    int equilibrationSteps = m_numberOfMetropolisSteps * m_equilibrationFraction;
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
        for (int i = equilibrationSteps; i < m_numberOfMetropolisSteps; i++) {
            bool acceptedStep = metropolisStep(private_particles, private_waveFuncValue);
            m_sampler->sample(acceptedStep, private_particles, thread_num);
        }
    }
    m_sampler->computeAverages();
    m_sampler->printOutputToTerminal();
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
