#include "system.h"
#include "sampler.h"
#include "particle.h"
#include "WaveFunctions/wavefunction.h"
#include "Hamiltonians/hamiltonian.h"
#include "InitialStates/initialstate.h"
#include "Math/random.h"
#include "omp.h"

#include <string>

using namespace std;

System::System(int num_threads, int seed) {
    if (num_threads == -1) {
        m_num_threads = omp_get_max_threads();
    } else {
        m_num_threads = num_threads;
        omp_set_num_threads(num_threads);
    }

    if (seed == -1) {
        for (int i = 0; i < m_num_threads; i++) {
            m_randoms.push_back(new Random());
        }
    } else {
        for (int i = 0; i < m_num_threads; i++) {
            m_randoms.push_back(new Random(seed + i));
        }
    }
}

bool System::ImmetropolisStep(vector<Particle*> particles, double& waveFuncValue) {
    Random* private_random = m_randoms.at(omp_get_thread_num());

    Particle* randParticle = particles[private_random->nextInt(0, m_numberOfParticles - 1)];
    vector<double> oldPos = randParticle->getPosition();
    double oldLengthSq = randParticle->getLengthSq();

    // Obtaining the vector contaning drift force entries from SimpleGaussian
    vector<double> oldQuantumForce = m_waveFunction->ComputeQF(randParticle,oldPos);
    // Declaring a vector to store the change of position in each spatial direction
    vector<double> move(m_numberOfDimensions, 0);

    // Looping over each spatial component to update the positions
    for (int i = 0; i < m_numberOfDimensions; i++) {
        move[i] = 0.5 * m_timestep * oldQuantumForce[i] + private_random->nextGaussian(0.,1.) * sqrt(m_timestep);
    }

    // Changing randParticle position
    randParticle->adjustLangevin(move);
    vector<double> newPos = randParticle->getPosition();

    double newWaveFuncValue = m_waveFunction->evaluateChange(randParticle, waveFuncValue, oldLengthSq);
    vector<double> newQuantumForce = m_waveFunction->ComputeQF(randParticle, newPos);
    double GreensFunctionRatio = 0.0;

    for (int i = 0; i < m_numberOfDimensions; i++) {
        GreensFunctionRatio += (oldQuantumForce[i] + newQuantumForce[i]) * (oldPos[i] - newPos[i] +
                                    m_timestep * (oldQuantumForce[i] - newQuantumForce[i]));
    }
    GreensFunctionRatio = exp(0.5 * GreensFunctionRatio);

    if (private_random->nextDouble() < GreensFunctionRatio * pow(newWaveFuncValue / waveFuncValue, 2)) {
        waveFuncValue = newWaveFuncValue;
        oldQuantumForce = newQuantumForce;
        return true;
    } else {
        randParticle->setPosition(oldPos);
        return false;
    }
}

bool System::metropolisStep(vector<Particle*> particles, double& waveFuncValue) {
    Random* private_random = m_randoms.at(omp_get_thread_num());

    Particle* randParticle = particles[private_random->nextInt(0, m_numberOfParticles - 1)];
    vector<double> oldPos = randParticle->getPosition();
    double oldLengthSq = randParticle->getLengthSq();

    double dir = (private_random->nextInt(1) - 0.5) * 2;

    randParticle->adjustPosition(m_stepLength * dir, private_random->nextInt(0, m_numberOfDimensions - 1));
    double newWaveFuncValue = m_waveFunction->evaluateChange(randParticle, waveFuncValue, oldLengthSq);

    if (private_random->nextDouble() < pow(newWaveFuncValue / waveFuncValue, 2)) {
        waveFuncValue = newWaveFuncValue;
        return true;
    } else {
        randParticle->setPosition(oldPos);
        return false;
    }
}

void System::runMetropolisSteps(bool m_choice) {
    int parallellSteps = m_metroSteps / m_num_threads + m_equiSteps;
    # pragma omp parallel
    {
        // Private variables setup for each thread
        vector<Particle*> private_particles = m_initialState->newParticles();
        double private_waveFuncValue = m_waveFunction->evaluate(private_particles);
        int thread_num = omp_get_thread_num();

        if (m_choice == 1) { // Importance sampling
            for (int i = 0; i < m_equiSteps; i++) { // Equilibration
                ImmetropolisStep(private_particles, private_waveFuncValue);
            }
            m_sampler->updateVals(private_particles, thread_num);
            for (int i = m_equiSteps; i < parallellSteps; i++) { // Steps which are sampled
                bool acceptedStep = ImmetropolisStep(private_particles, private_waveFuncValue);
                m_sampler->sample(acceptedStep, private_particles, thread_num);
            }
        } else { // Standard metropolis
            for (int i = 0; i < m_equiSteps; i++) { // Equilibration
                metropolisStep(private_particles, private_waveFuncValue);
            }
            m_sampler->updateVals(private_particles, thread_num);
            for (int i = m_equiSteps; i < parallellSteps; i++) { // Steps which are sampled
                bool acceptedStep = metropolisStep(private_particles, private_waveFuncValue);
                m_sampler->sample(acceptedStep, private_particles, thread_num);
            }
        }
    }
    m_sampler->computeAverages();
    addOutput(m_sampler->outputText());
}

void System::addOutput(string text) {
    m_output += text;
}

void System::clearOutput() {
    m_output = "";
}

// ----------------- Setters and getters ----------------------

void System::setNumberOfSteps(int metroSteps, int equiSteps) {
    m_metroSteps = metroSteps;
    m_equiSteps = equiSteps;
}

void System::setNumberOfParticles(int numberOfParticles) {
    m_numberOfParticles = numberOfParticles;
}

void System::setNumberOfDimensions(int numberOfDimensions) {
    m_numberOfDimensions = numberOfDimensions;
}

void System::setStepLength(double stepLength) {
    m_stepLength = stepLength;
}

void System::setTimeStep(double timestep) {
    m_timestep = timestep;
}

void System::setChoice(bool ImpSampling) {
    m_choice = ImpSampling;
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

void System::setSampler(Sampler* sampler) {
    m_sampler = sampler;
}

Random* System::getRandomEngine() {
    return m_randoms.at(omp_get_thread_num());
}
