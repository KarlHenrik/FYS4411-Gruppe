#include "system.h"
#include <cassert>
#include "sampler.h"
#include "particle.h"
#include "WaveFunctions/wavefunction.h"
#include "Hamiltonians/hamiltonian.h"
#include "InitialStates/initialstate.h"
#include "Math/random.h"


System::System() {
    m_random = new Random();
}

System::System(int seed) {
    m_random = new Random(seed);
}

bool System::metropolisStep() {
    /* Perform the actual Metropolis step: Choose a particle at random and
     * change it's position by a random amount, and check if the step is
     * accepted by the Metropolis test (compare the wave function evaluated
     * at this new position with the one at the old position).
     */
    Particle* randParticle = m_particles[m_random->nextInt(0, m_numberOfParticles - 1)];
    std::vector<double> oldPos = randParticle->getPosition();
    double oldLengthSq = randParticle->getLengthSq();

    randParticle->adjustPosition(m_stepLength, m_random->nextInt(0, m_numberOfDimensions - 1));
    double newWaveFuncValue = m_waveFunction->evaluateChange(randParticle, m_waveFunctionValue, oldLengthSq);

    if (m_random->nextDouble() > std::pow(newWaveFuncValue / m_waveFunctionValue, 2)) {
        m_waveFunctionValue = newWaveFuncValue;
        return true;
    } else {
        randParticle->setPosition(oldPos);
        return false;
    }
}

void System::runMetropolisSteps(int numberOfMetropolisSteps) {
    m_particles                 = m_initialState->getParticles();
    m_sampler                   = new Sampler(this);
    m_numberOfMetropolisSteps   = numberOfMetropolisSteps;
    m_sampler->setNumberOfMetropolisSteps(numberOfMetropolisSteps);
    m_waveFunctionValue = m_waveFunction->evaluate(m_particles);
    /* Here you should sample the energy (and maybe other things using
     * the m_sampler instance of the Sampler class. Make sure, though,
     * to only begin sampling after you have let the system equilibrate
     * for a while. You may handle this using the fraction of steps which
     * are equilibration steps; m_equilibrationFraction.
     */
    int equilibrationSteps = numberOfMetropolisSteps * m_equilibrationFraction;
    for (int i = 0; i < equilibrationSteps; i++) {
        metropolisStep();
    }
    for (int i = equilibrationSteps; i < numberOfMetropolisSteps; i++) {
        bool acceptedStep = metropolisStep();

        m_sampler->sample(acceptedStep);
    }

    m_sampler->computeAverages();
    m_sampler->printOutputToTerminal();
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
