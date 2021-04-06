#include "catch.hpp"

#include "../system.h"
#include "../WaveFunctions/wavefunction.h"
#include "../WaveFunctions/simplegaussian.h"
#include "../WaveFunctions/correlated.h"
#include "../Hamiltonians/hamiltonian.h"
#include "../Hamiltonians/harmonicoscillator.h"
#include "../InitialStates/initialstate.h"
#include "../InitialStates/randomuniform.h"
#include "../paramTester.h"
#include "../sampler.h"

#include <string>

#include <iostream>
#include <cmath> // log10
#include <string> // output
#include <iomanip>
#include <fstream>
#include <sstream>

using namespace std;

TEST_CASE("param derivative is zero at 0.5 for simple gaussion", "[sampler->getParamDer()]") {
    // Physical system parameters
    int numberOfDimensions  = 1;
    int numberOfParticles   = 2;
    double omega            = 1.0;       // Oscillator frequency
    double alpha            = 0.5;       // Variational parameter, initial value
    double a                = 0;     // Interaction radius
    double beta             = 1;   // Shift of trap in the z-direction
    // Metropolis parameters
    int metroSteps = 10000;
    int equiSteps = 10000;           // Steps are specified down with the type of calculation
    bool ImpSampling        = 1;         // cout << "Perform importance sampling? [0,1] : "; cin >> ImpSampling;
    double stepLength       = 0.1;       // Metropolis step length used without importance sampling
    double timestep         = 0.1;      // Time step used in importance sampling movement
    // Other parameters
    int seed                = 42;        // Seed for the random number generator. -1 means random seed
    int num_threads         = -1;        // Number of threads for calculation. -1 means max (automatic)
    string systemName = "temp/temp";

    // ---------------SYSTEM SETUP-----------------------
    System* system = new System(num_threads, seed);
    system->setHamiltonian    (new HarmonicOscillator(system, omega, beta));
    system->setWaveFunction   (new SimpleGaussian(system, alpha, a, beta));
    system->setInitialState   (new RandomUniform(system, numberOfDimensions, numberOfParticles, a));
    system->setStepLength     (stepLength);
    system->setTimeStep       (timestep);
    system->setChoice         (ImpSampling);

    // Testing sampler
    system->setNumberOfSteps(metroSteps, equiSteps);
    Sampler* sampler = new Sampler(system, system->getNumerOfThreads());
    system->setSampler(sampler);

    system->runMetropolisSteps(system->getChoice());

    REQUIRE( sampler->getParamDer() == Approx(0).epsilon(0.00001) );
}

TEST_CASE("param derivative is non-zero at 0.5 for simple gaussion with beta != 1", "[sampler->getParamDer()]") {
    // Physical system parameters
    int numberOfDimensions  = 3;
    int numberOfParticles   = 2;
    double omega            = 1.0;       // Oscillator frequency
    double alpha            = 0.5;       // Variational parameter, initial value
    double a                = 0;     // Interaction radius
    double beta             = sqrt(8);   // Shift of trap in the z-direction
    // Metropolis parameters
    int metroSteps = 100000;
    int equiSteps = 10000;           // Steps are specified down with the type of calculation
    bool ImpSampling        = 1;         // cout << "Perform importance sampling? [0,1] : "; cin >> ImpSampling;
    double stepLength       = 0.1;       // Metropolis step length used without importance sampling
    double timestep         = 0.01;      // Time step used in importance sampling movement
    // Other parameters
    int seed                = 42;        // Seed for the random number generator. -1 means random seed
    int num_threads         = -1;        // Number of threads for calculation. -1 means max (automatic)
    string systemName = "temp/temp";

    // ---------------SYSTEM SETUP-----------------------
    System* system = new System(num_threads, seed);
    system->setHamiltonian    (new HarmonicOscillator(system, omega, beta));
    system->setWaveFunction   (new SimpleGaussian(system, alpha, a, beta));
    system->setInitialState   (new RandomUniform(system, numberOfDimensions, numberOfParticles, a));
    system->setStepLength     (stepLength);
    system->setTimeStep       (timestep);
    system->setChoice         (ImpSampling);

    // Testing sampler
    system->setNumberOfSteps(metroSteps, equiSteps);
    Sampler* sampler = new Sampler(system, system->getNumerOfThreads());
    system->setSampler(sampler);

    system->runMetropolisSteps(system->getChoice());

    REQUIRE( sampler->getParamDer() != Approx(0).epsilon(0.00001) );
}

TEST_CASE("param derivative points toward 0.5 for simple gaussion", "[sampler->getParamDer()]") {
    // Physical system parameters
    int numberOfDimensions  = 3;
    int numberOfParticles   = 2;
    double omega            = 1.0;       // Oscillator frequency
    double alpha            = 0.4;       // Variational parameter, initial value
    double a                = 0;     // Interaction radius
    double beta             = 1;   // Shift of trap in the z-direction
    // Metropolis parameters
    int metroSteps = 100000;
    int equiSteps = 10000;           // Steps are specified down with the type of calculation
    bool ImpSampling        = 1;         // cout << "Perform importance sampling? [0,1] : "; cin >> ImpSampling;
    double stepLength       = 0.1;       // Metropolis step length used without importance sampling
    double timestep         = 0.01;      // Time step used in importance sampling movement
    // Other parameters
    int seed                = 42;        // Seed for the random number generator. -1 means random seed
    int num_threads         = -1;        // Number of threads for calculation. -1 means max (automatic)
    string systemName = "temp/temp";

    // ---------------SYSTEM SETUP-----------------------
    System* system = new System(num_threads, seed);
    system->setHamiltonian    (new HarmonicOscillator(system, omega, beta));
    system->setWaveFunction   (new SimpleGaussian(system, alpha, a, beta));
    system->setInitialState   (new RandomUniform(system, numberOfDimensions, numberOfParticles, a));
    system->setStepLength     (stepLength);
    system->setTimeStep       (timestep);
    system->setChoice         (ImpSampling);
    system->setNumberOfSteps(metroSteps, equiSteps);

    // Testing alpha less than 0.5
    Sampler* sampler = new Sampler(system, system->getNumerOfThreads());
    system->setSampler(sampler);

    system->runMetropolisSteps(system->getChoice());
    REQUIRE( sampler->getParamDer() < 0 );

    // Testing alpha greater than 0.5
    alpha = 0.6;
    system->setWaveFunction   (new SimpleGaussian(system, alpha, a, beta));

    sampler = new Sampler(system, system->getNumerOfThreads());
    system->setSampler(sampler);

    system->runMetropolisSteps(system->getChoice());
    REQUIRE( sampler->getParamDer() > 0 );
}

TEST_CASE("energy increases with number of dimensitons", "[sampler->getEnergy()]") {
    // Physical system parameters
    int numberOfDimensions  = 1;
    int numberOfParticles   = 2;
    double omega            = 1.0;       // Oscillator frequency
    double alpha            = 0.5;       // Variational parameter, initial value
    double a                = 0;     // Interaction radius
    double beta             = 1;   // Shift of trap in the z-direction
    // Metropolis parameters
    int metroSteps = 10000;
    int equiSteps = 10000;           // Steps are specified down with the type of calculation
    bool ImpSampling        = 1;         // cout << "Perform importance sampling? [0,1] : "; cin >> ImpSampling;
    double stepLength       = 0.1;       // Metropolis step length used without importance sampling
    double timestep         = 0.1;      // Time step used in importance sampling movement
    // Other parameters
    int seed                = 42;        // Seed for the random number generator. -1 means random seed
    int num_threads         = -1;        // Number of threads for calculation. -1 means max (automatic)
    string systemName = "temp/temp";

    // ---------------SYSTEM SETUP-----------------------
    System* system = new System(num_threads, seed);
    system->setHamiltonian    (new HarmonicOscillator(system, omega, beta));
    system->setWaveFunction   (new SimpleGaussian(system, alpha, a, beta));
    system->setInitialState   (new RandomUniform(system, numberOfDimensions, numberOfParticles, a));
    system->setStepLength     (stepLength);
    system->setTimeStep       (timestep);
    system->setChoice         (ImpSampling);
    system->setNumberOfSteps(metroSteps, equiSteps);

    // 1 dim
    Sampler* sampler = new Sampler(system, system->getNumerOfThreads());
    system->setSampler(sampler);
    system->runMetropolisSteps(system->getChoice());
    double dim1 = sampler->getEnergy();

    // 2 dim
    system->setInitialState   (new RandomUniform(system, 2, numberOfParticles, a));
    sampler = new Sampler(system, system->getNumerOfThreads());
    system->setSampler(sampler);
    system->runMetropolisSteps(system->getChoice());
    double dim2 = sampler->getEnergy();

    // 3 dim
    system->setInitialState   (new RandomUniform(system, 3, numberOfParticles, a));
    sampler = new Sampler(system, system->getNumerOfThreads());
    system->setSampler(sampler);
    system->runMetropolisSteps(system->getChoice());
    double dim3 = sampler->getEnergy();

    REQUIRE( dim1 < dim2 );
    REQUIRE( dim2 < dim3 );
}

TEST_CASE("particles are placed correctly", "[RandomUniform()]") {
    // Physical system parameters
    int numberOfDimensions  = 3;
    int numberOfParticles   = 500;
    double omega            = 1.0;       // Oscillator frequency
    double alpha            = 0.5;       // Variational parameter, initial value
    double a                = 0.01;     // Interaction radius
    double beta             = 1;   // Shift of trap in the z-direction
    // Metropolis parameters
    int metroSteps = 10000;
    int equiSteps = 10000;           // Steps are specified down with the type of calculation
    bool ImpSampling        = 1;         // cout << "Perform importance sampling? [0,1] : "; cin >> ImpSampling;
    double stepLength       = 0.1;       // Metropolis step length used without importance sampling
    double timestep         = 0.1;      // Time step used in importance sampling movement
    // Other parameters
    int seed                = 42;        // Seed for the random number generator. -1 means random seed
    int num_threads         = -1;        // Number of threads for calculation. -1 means max (automatic)
    string systemName = "temp/temp";

    // ---------------SYSTEM SETUP-----------------------
    System* system = new System(num_threads, seed);
    system->setHamiltonian    (new HarmonicOscillator(system, omega, beta));
    system->setWaveFunction   (new SimpleGaussian(system, alpha, a, beta));
    system->setInitialState   (new RandomUniform(system, numberOfDimensions, numberOfParticles, a));
    system->setStepLength     (stepLength);
    system->setTimeStep       (timestep);
    system->setChoice         (ImpSampling);
    system->setNumberOfSteps(metroSteps, equiSteps);

    double mistakes = 0;
    vector<Particle*> particles = system->getInitialState()->newParticles();
    for (int p1 = 0; p1 < (int) particles.size(); p1++) {
        for (int d = 0; d < numberOfDimensions; d++) {
            if (abs(particles[p1]->getPosition()[d]) > 1) { // All particles within [-1, 1]
                mistakes++;
            }
        }

        for (int p2 = 0; p2 < (int) particles.size(); p2++) {
            double distSq = 0;
            for (int d = 0; d < numberOfDimensions; d++) {
                distSq += pow(particles[p1]->getPosition()[d] - particles[p2]->getPosition()[d], 2);
            }
            if (distSq < a * a && p1 != p2) { // All particles at least a apart
                mistakes++;
            }
        }
    }
    REQUIRE( mistakes == 0 );
}

TEST_CASE("particles have their length squared updated when moved", "[Particle->setPosition()]") {
    // Physical system parameters
    int numberOfDimensions  = 3;
    int numberOfParticles   = 500;
    double omega            = 1.0;       // Oscillator frequency
    double alpha            = 0.5;       // Variational parameter, initial value
    double a                = 0.01;     // Interaction radius
    double beta             = 1;   // Shift of trap in the z-direction
    // Metropolis parameters
    int metroSteps = 10000;
    int equiSteps = 10000;           // Steps are specified down with the type of calculation
    bool ImpSampling        = 1;         // cout << "Perform importance sampling? [0,1] : "; cin >> ImpSampling;
    double stepLength       = 0.1;       // Metropolis step length used without importance sampling
    double timestep         = 0.1;      // Time step used in importance sampling movement
    // Other parameters
    int seed                = 42;        // Seed for the random number generator. -1 means random seed
    int num_threads         = -1;        // Number of threads for calculation. -1 means max (automatic)
    string systemName = "temp/temp";

    // ---------------SYSTEM SETUP-----------------------
    System* system = new System(num_threads, seed);
    system->setHamiltonian    (new HarmonicOscillator(system, omega, beta));
    system->setWaveFunction   (new SimpleGaussian(system, alpha, a, beta));
    system->setInitialState   (new RandomUniform(system, numberOfDimensions, numberOfParticles, a));
    system->setStepLength     (stepLength);
    system->setTimeStep       (timestep);
    system->setChoice         (ImpSampling);
    system->setNumberOfSteps(metroSteps, equiSteps);

    vector<Particle*> particles = system->getInitialState()->newParticles();
    Particle* p = particles[0];
    vector<double> oldPos = p->getPosition();
    oldPos[0] += 1;

    double l = p->getLengthSq();
    p->setPosition(oldPos);
    REQUIRE( l != p->getLengthSq() );

    l = p->getLengthSq();
    p->adjustPosition(0.3, 0);
    REQUIRE( l != p->getLengthSq() );

    l = p->getLengthSq();
    p->adjustLangevin(oldPos);
    REQUIRE( l != p->getLengthSq() );
}

TEST_CASE("Correlated WaveFunction does not have minumum at 0.5", "[Particle->setPosition()]") {
    // Physical system parameters
    int numberOfDimensions  = 1;
    int numberOfParticles   = 2;
    double omega            = 1.0;       // Oscillator frequency
    double alpha            = 0.5;       // Variational parameter, initial value
    double a                = 0.01;     // Interaction radius
    double beta             = 1;   // Shift of trap in the z-direction
    // Metropolis parameters
    int metroSteps = 10000;
    int equiSteps = 10000;           // Steps are specified down with the type of calculation
    bool ImpSampling        = 1;         // cout << "Perform importance sampling? [0,1] : "; cin >> ImpSampling;
    double stepLength       = 0.1;       // Metropolis step length used without importance sampling
    double timestep         = 0.1;      // Time step used in importance sampling movement
    // Other parameters
    int seed                = 42;        // Seed for the random number generator. -1 means random seed
    int num_threads         = -1;        // Number of threads for calculation. -1 means max (automatic)
    string systemName = "temp/temp";

    // ---------------SYSTEM SETUP-----------------------
    System* system = new System(num_threads, seed);
    system->setHamiltonian    (new HarmonicOscillator(system, omega, beta));
    system->setWaveFunction   (new Correlated(system, alpha, a, beta));
    system->setInitialState   (new RandomUniform(system, numberOfDimensions, numberOfParticles, a));
    system->setStepLength     (stepLength);
    system->setTimeStep       (timestep);
    system->setChoice         (ImpSampling);

    // Testing sampler
    system->setNumberOfSteps(metroSteps, equiSteps);
    Sampler* sampler = new Sampler(system, system->getNumerOfThreads());
    system->setSampler(sampler);

    system->runMetropolisSteps(system->getChoice());

    REQUIRE( sampler->getParamDer() != 0 );
}
