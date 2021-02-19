#include <iostream>
#include "system.h"
#include "particle.h"
#include "WaveFunctions/wavefunction.h"
#include "WaveFunctions/simplegaussian.h"
#include "Hamiltonians/hamiltonian.h"
#include "Hamiltonians/harmonicoscillator.h"
#include "InitialStates/initialstate.h"
#include "InitialStates/randomuniform.h"
#include "Math/random.h"
#include "paramTester.h"
#include <string>

using namespace std;


int main() {

    { // Simplest possible system simulation
        // Seed for the random number generator
        int seed = 2020;
        // Parameters for system
        int numberOfDimensions  = 1;
        int numberOfParticles   = 1;
        int numberOfSteps       = (int) 1e6;
        double omega            = 1.0;          // Oscillator frequency.
        double alpha            = 0.4;        // Variational parameter.
        double stepLength       = 0.1;          // Metropolis step length.
        double equilibration    = 0.375;          // Amount of the total steps used for equilibration.

        // System setup
        System* system = new System(6); // no seed usage
        system->setHamiltonian              (new HarmonicOscillator(system, omega));
        system->setWaveFunction             (new SimpleGaussian(system, alpha));
        system->setInitialState             (new RandomUniform(system, numberOfDimensions, numberOfParticles));
        system->setNumberOfSteps            (numberOfSteps);
        system->setEquilibrationFraction    (equilibration);
        system->setStepLength               (stepLength);

        // Alpha testing, parameters and setup
        double alpha_end = 0.6;
        double alpha_step = 0.01;

        ParamTester* paramTester = new ParamTester(system, "HO_Gauss_RU.txt");
        paramTester->alphaGrid(alpha, alpha_end, alpha_step);
    }

    return 0;
}
