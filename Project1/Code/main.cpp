#include "system.h"
#include "WaveFunctions/wavefunction.h"
#include "WaveFunctions/simplegaussian.h"
#include "WaveFunctions/correlated.h"
#include "Hamiltonians/hamiltonian.h"
#include "Hamiltonians/harmonicoscillator.h"
#include "InitialStates/initialstate.h"
#include "InitialStates/randomuniform.h"
#include "paramTester.h"

#include <string>

using namespace std;

void testing();
void energiPerAlpha();
void corrResultat();
void benchmark(); // DONE!
void benchmarkDim(); // DONE!
void medUtenImp();
void medUtenParallell();

int main() {
    //testing();
    return 0;
}

void testing() {
    { // Simplest possible system simulation
        // ----------------SYSTEM PARAMETERS---------------------
        // Physical system parameters
        int numberOfDimensions  = 2;
        int numberOfParticles   = 1;
        double omega            = 1.0;       // Oscillator frequency
        double alpha            = 0.45;       // Variational parameter, initial value
        double a                = 0.0042;     // Interaction radius
        double beta             = 1;   // Shift of trap in the z-direction
        // Metropolis parameters
        int metroSteps, equiSteps;           // Steps are specified down with the type of calculation
        bool ImpSampling        = 1;         // cout << "Perform importance sampling? [0,1] : "; cin >> ImpSampling;
        double stepLength       = 0.05;       // Metropolis step length used without importance sampling
        double timestep         = 0.01;      // Time step used in importance sampling movement
        // Other parameters
        int seed                = -1;        // Seed for the random number generator. -1 means random seed
        int num_threads         = -1;        // Number of threads for calculation. -1 means max (automatic)
        string systemName = "HO_Gauss_RU";

        // ---------------SYSTEM SETUP-----------------------
        System* system = new System(num_threads, seed);
        system->setHamiltonian              (new HarmonicOscillator(system, omega));
        system->setWaveFunction             (new Correlated(system, alpha, a, beta));
        system->setInitialState             (new RandomUniform(system, numberOfDimensions, numberOfParticles, a));
        system->setStepLength               (stepLength);
        system->setTimeStep                 (timestep);
        system->setChoice                   (ImpSampling);

        // -------------Calculations---------------------
        ParamTester* paramTester = new ParamTester(system, systemName);

        // Doing grid search to produce plot of energy for each alpha
        metroSteps = (int) 1e6; // Number of metropolis steps
        equiSteps = (int) 1e5;  // Amount of the total steps used for equilibration
        system->setNumberOfSteps(metroSteps, equiSteps);
        double alpha_end = 0.55;    // final alpha parameter to test, the first is the alpha defined earlier
        double alpha_step = 0.01;  // step length in alpha search
        paramTester->alphaGrid(alpha, alpha_end, alpha_step);

        // Using gradient descent to find optimal alpha. Also prints results along the way.
        metroSteps = (int) 1e5; // Number of metropolis steps
        equiSteps = (int) 1e4;  // Amount of the total steps used for equilibration
        system->setNumberOfSteps(metroSteps, equiSteps);
        double lr = 0.1;
        double tol = 0.0000001;
        double max_iter = 20;
        double alpha_opt = paramTester->alphaGD(alpha, lr, tol, max_iter);

        // Doing a large scale calculation for optimal alpha
        metroSteps = (int) 1e7; // Number of metropolis steps
        equiSteps = (int) 1e6;    // Amount of the total steps used for equilibration
        system->setNumberOfSteps(metroSteps, equiSteps);
        paramTester->bigCalc(alpha_opt);
    }
}

void energiPerAlpha() {
    return;
}

void corrResultat() {
    return;
}

void benchmark() {
    double numPart[4] = {1, 10, 100, 500};
    for (int i = 0; i < 4; i++) {
        int numberOfParticles = numPart[i];
        // ----------------SYSTEM PARAMETERS---------------------
        // Physical system parameters
        int numberOfDimensions  = 3;
        //int numberOfParticles   = 1;
        double omega            = 1.0;       // Oscillator frequency
        double alpha            = 0.1;       // Variational parameter, initial value
        double a                = 0.0042;     // Interaction radius
        double beta             = 1;   // Shift of trap in the z-direction
        // Metropolis parameters
        int metroSteps, equiSteps;           // Steps are specified down with the type of calculation
        bool ImpSampling        = 1;         // cout << "Perform importance sampling? [0,1] : "; cin >> ImpSampling;
        double stepLength       = 0.1;       // Metropolis step length used without importance sampling
        double timestep         = 0.1;      // Time step used in importance sampling movement
        // Other parameters
        int seed                = 42;        // Seed for the random number generator. -1 means random seed
        int num_threads         = -1;        // Number of threads for calculation. -1 means max (automatic)
        string systemName = "Benchmarks/Bench" + to_string(numberOfParticles);

        // ---------------SYSTEM SETUP-----------------------
        System* system = new System(num_threads, seed);
        system->setHamiltonian    (new HarmonicOscillator(system, omega));
        system->setWaveFunction   (new SimpleGaussian(system, alpha, a, beta));
        system->setInitialState   (new RandomUniform(system, numberOfDimensions, numberOfParticles, a));
        system->setStepLength     (stepLength);
        system->setTimeStep       (timestep);
        system->setChoice         (ImpSampling);

        // -------------Calculations---------------------
        ParamTester* paramTester = new ParamTester(system, systemName);

        // Using gradient descent to find optimal alpha. Also prints results along the way.
        metroSteps = (int) 2e4; // Number of metropolis steps
        equiSteps = (int) 1e4;  // Amount of the total steps used for equilibration
        system->setNumberOfSteps(metroSteps, equiSteps);
        double lr = 0.02;
        double tol = 0.00000001;
        double max_iter = 200;
        double alpha_opt = paramTester->alphaGD(alpha, lr, tol, max_iter);

        // Doing a large scale calculation for optimal alpha
        metroSteps = (int) 1e7; // Number of metropolis steps
        equiSteps = (int) 1e5;    // Amount of the total steps used for equilibration
        system->setNumberOfSteps(metroSteps, equiSteps);
        paramTester->bigCalc(alpha_opt);

    }
}

void benchmarkDim() {
    double numDims[3] = {1, 2, 3};
    for (int i = 0; i < 3; i++) {
        int numberOfDimensions = numDims[i];
        // ----------------SYSTEM PARAMETERS---------------------
        // Physical system parameters
        //int numberOfDimensions  = 3;
        int numberOfParticles   = 10;
        double omega            = 1.0;       // Oscillator frequency
        double alpha            = 0.1;       // Variational parameter, initial value
        double a                = 0.0042;     // Interaction radius
        double beta             = 1;   // Shift of trap in the z-direction
        // Metropolis parameters
        int metroSteps, equiSteps;           // Steps are specified down with the type of calculation
        bool ImpSampling        = 1;         // cout << "Perform importance sampling? [0,1] : "; cin >> ImpSampling;
        double stepLength       = 0.1;       // Metropolis step length used without importance sampling
        double timestep         = 0.1;      // Time step used in importance sampling movement
        // Other parameters
        int seed                = 42;        // Seed for the random number generator. -1 means random seed
        int num_threads         = -1;        // Number of threads for calculation. -1 means max (automatic)
        string systemName = "Benchmarks/BenchDim" + to_string(numberOfDimensions);

        // ---------------SYSTEM SETUP-----------------------
        System* system = new System(num_threads, seed);
        system->setHamiltonian    (new HarmonicOscillator(system, omega));
        system->setWaveFunction   (new SimpleGaussian(system, alpha, a, beta));
        system->setInitialState   (new RandomUniform(system, numberOfDimensions, numberOfParticles, a));
        system->setStepLength     (stepLength);
        system->setTimeStep       (timestep);
        system->setChoice         (ImpSampling);

        // -------------Calculations---------------------
        ParamTester* paramTester = new ParamTester(system, systemName);

        // Using gradient descent to find optimal alpha. Also prints results along the way.
        metroSteps = (int) 2e4; // Number of metropolis steps
        equiSteps = (int) 1e4;  // Amount of the total steps used for equilibration
        system->setNumberOfSteps(metroSteps, equiSteps);
        double lr = 0.1;
        double tol = 0.00000001;
        double max_iter = 200;
        double alpha_opt = paramTester->alphaGD(alpha, lr, tol, max_iter);

        // Doing a large scale calculation for optimal alpha
        metroSteps = (int) 1e7; // Number of metropolis steps
        equiSteps = (int) 1e5;    // Amount of the total steps used for equilibration
        system->setNumberOfSteps(metroSteps, equiSteps);
        paramTester->bigCalc(alpha_opt);

    }
}

void medUtenImp() {
    return;
}

void medUtenParallell() {
    return;
}
