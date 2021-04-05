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

#include <iostream>
#include <cmath> // log10
#include <string> // output
#include <iomanip>
#include <fstream>
#include <sstream>

using namespace std;

void testing();
void energyPerAlpha(); // DONE!
void corrResultat(); // DOING!
void benchmark(); // DONE!
void benchmarkDim(); // DONE!
void medUtenImp();
void medUtenParallell();
void loopinteraction();

int main() {
    medUtenParallell();
    return 0;
}

void testing() {
    { // Simplest possible system simulation
        // ----------------SYSTEM PARAMETERS---------------------
        // Physical system parameters
        int numberOfDimensions  = 3;
        int numberOfParticles   = 10;
        double omega            = 1.0;       // Oscillator frequency
        double alpha            = 0.5;       // Variational parameter, initial value
        double a                = 0.00433;     // Interaction radius
        double beta             = 1;   // Shift of trap in the z-direction
        // Metropolis parameters
        int metroSteps, equiSteps;           // Steps are specified down with the type of calculation
        bool ImpSampling        = 1;         // cout << "Perform importance sampling? [0,1] : "; cin >> ImpSampling;
        double stepLength       = 0.05;       // Metropolis step length used without importance sampling
        double timestep         = 0.01;      // Time step used in importance sampling movement
        // Other parameters
        int seed                = 42;        // Seed for the random number generator. -1 means random seed
        int num_threads         = -1;        // Number of threads for calculation. -1 means max (automatic)
        string systemName = "CorrOB_test";

        // ---------------SYSTEM SETUP-----------------------
        System* system = new System(num_threads, seed);
        system->setHamiltonian              (new HarmonicOscillator(system, omega, beta));
        system->setWaveFunction             (new Correlated(system, alpha, a, beta));
        system->setInitialState             (new RandomUniform(system, numberOfDimensions, numberOfParticles, a));
        system->setStepLength               (stepLength);
        system->setTimeStep                 (timestep);
        system->setChoice                   (ImpSampling);

        // -------------Calculations---------------------
        ParamTester* paramTester = new ParamTester(system, systemName);
        /*
        // Doing grid search to produce plot of energy for each alpha
        metroSteps = (int) 1e6; // Number of metropolis steps
        equiSteps = (int) 1e5;  // Amount of the total steps used for equilibration
        system->setNumberOfSteps(metroSteps, equiSteps);
        double alpha_end = 0.4989;    // final alpha parameter to test, the first is the alpha defined earlier
        double alpha_step = 0.0001;  // step length in alpha search
        paramTester->alphaGrid(alpha, alpha_end, alpha_step);
        */
        // Using gradient descent to find optimal alpha. Also prints results along the way.
        /*
        metroSteps = (int) 1e5; // Number of metropolis steps
        equiSteps = (int) 1e4;  // Amount of the total steps used for equilibration
        system->setNumberOfSteps(metroSteps, equiSteps);
        double lr = 0.01;
        double tol = 0.0000001;
        double max_iter = 200;
        double alpha_opt = paramTester->alphaGD(alpha, lr, tol, max_iter);
        */

        //double alpha_opt = 0.468;
        double alpha_opt = 0.497692;
        // Doing a large scale calculation for optimal alpha
        metroSteps = (int) 1e7; // Number of metropolis steps
        equiSteps = (int) 1e5;    // Amount of the total steps used for equilibration
        system->setNumberOfSteps(metroSteps, equiSteps);
        paramTester->bigCalc(alpha_opt);

    }
}

void energyPerAlpha() {
    cout << "Running simulations to find energy vs. alpha relationship" << endl;
    // ----------------SYSTEM PARAMETERS---------------------
    // Physical system parameters
    int numberOfDimensions  = 3;
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
    string systemName = "EnergyPerAlpha/EnergyPerAlpha";

    // ---------------SYSTEM SETUP-----------------------
    System* system = new System(num_threads, seed);
    system->setHamiltonian    (new HarmonicOscillator(system, omega, beta));
    system->setWaveFunction   (new SimpleGaussian(system, alpha, a, beta));
    system->setInitialState   (new RandomUniform(system, numberOfDimensions, numberOfParticles, a));
    system->setStepLength     (stepLength);
    system->setTimeStep       (timestep);
    system->setChoice         (ImpSampling);

    // -------------Calculations---------------------
    ParamTester* paramTester = new ParamTester(system, systemName);

    // Doing grid search to produce plot of energy for each alpha
    metroSteps = (int) 2e6; // Number of metropolis steps
    equiSteps = (int) 1e5;  // Amount of the total steps used for equilibration
    system->setNumberOfSteps(metroSteps, equiSteps);
    double alpha_end = 0.9;    // final alpha parameter to test, the first is the alpha defined earlier
    double alpha_step = 0.005;  // step length in alpha search
    paramTester->alphaGrid(alpha, alpha_end, alpha_step);

    double lrs[5] = {0.001, 0.01, 0.1, 0.2, 0.5};
    for (int i = 0; i < 5; i++) {
        paramTester = new ParamTester(system, "EnergyPerAlpha/EnergyPerLr" + to_string(i));
        // Using gradient descent to find optimal alpha. Also prints results along the way.
        metroSteps = (int) 2e4; // Number of metropolis steps
        equiSteps = (int) 1e4;  // Amount of the total steps used for equilibration
        system->setNumberOfSteps(metroSteps, equiSteps);
        double lr = lrs[i];
        double tol = 0.00000001;
        double max_iter = 200;
        double alpha_opt = paramTester->alphaGD(alpha, lr, tol, max_iter);
    }
    return;
}

void corrResultat() {
    cout << "Running iterative simulations for different numbers of bosons.";
    cout << endl;


    double numPart[5] = {2,3,5,10,20};
    for (int i = 0; i < 5; i++) {
        int numberOfParticles = numPart[i];
        // ----------------SYSTEM PARAMETERS---------------------
        // Physical system parameters
        int numberOfDimensions  = 3;
        //int numberOfParticles   = 1;
        double omega            = 1.0;       // Oscillator frequency
        double alpha            = 0.5;       // Initial guess taken from non-int case
        double a                = 0.00433;   // Interaction radius
        double beta             = 1.;       // Shift of trap in the z-direction
        // Metropolis parameters
        int metroSteps, equiSteps;           // Steps are specified down with the type of calculation
        bool ImpSampling        = 1;         // cout << "Perform importance sampling? [0,1] : "; cin >> ImpSampling;
        double stepLength       = 0.1;       // Metropolis step length used without importance sampling
        double timestep         = 0.01;     // Time step used in importance sampling movement
        // Other parameters
        vector<double> alphas(10, 0);
        for (int i = 0; i < alphas.size(); i++ ) {
            int seed                = i;        // Seed for the random number generator. -1 means random seed
            int num_threads         = -1;        // Number of threads for calculation. -1 means max (automatic)
            string systemName = "BigSims" + to_string(numberOfParticles);

            // ---------------SYSTEM SETUP-----------------------
            System* system = new System(num_threads, seed);
            system->setHamiltonian    (new HarmonicOscillator(system, omega, beta));
            system->setWaveFunction   (new Correlated(system, alpha, a, beta));
            system->setInitialState   (new RandomUniform(system, numberOfDimensions, numberOfParticles, a));
            system->setStepLength     (stepLength);
            system->setTimeStep       (timestep);
            system->setChoice         (ImpSampling);

            // -------------Calculations---------------------
            ParamTester* paramTester = new ParamTester(system, systemName);

            // Using gradient descent to find optimal alpha. Also prints results along the way.
            metroSteps = (int) 1e5; // Number of metropolis steps
            equiSteps = (int) 1e4;  // Amount of the total steps used for equilibration
            system->setNumberOfSteps(metroSteps, equiSteps);
            double lr = 0.01;
            double tol = 1e-6;
            double max_iter = 200;
            double alpha_opt = paramTester->alphaGD(alpha, lr, tol, max_iter);
            alphas[i] = alpha_opt;
        }
        double sum = 0.;
        double ave;
        for (int i = 0; i < alphas.size(); i++) {
            sum += alphas[i];
        }
        ave = sum / alphas.size();
        cout << "Average optimal alpha for " + to_string(numberOfParticles) + " bosons: " << ave << endl;

        int seed                = i;        // Seed for the random number generator. -1 means random seed
        int num_threads         = -1;        // Number of threads for calculation. -1 means max (automatic)
        string systemName = "BigSims" + to_string(numberOfParticles);

        // ---------------SYSTEM SETUP-----------------------
        System* system = new System(num_threads, seed);
        system->setHamiltonian    (new HarmonicOscillator(system, omega, beta));
        system->setWaveFunction   (new Correlated(system, alpha, a, beta));
        system->setInitialState   (new RandomUniform(system, numberOfDimensions, numberOfParticles, a));
        system->setStepLength     (stepLength);
        system->setTimeStep       (timestep);
        system->setChoice         (ImpSampling);

        // -------------Calculations---------------------
        ParamTester* paramTester = new ParamTester(system, systemName);
        // Doing a large scale calculation for optimal alpha
        metroSteps = (int) 1e7; // Number of metropolis steps
        equiSteps = (int) 1e5;    // Amount of the total steps used for equilibration
        system->setNumberOfSteps(metroSteps, equiSteps);
        paramTester->bigCalc(ave);

    }
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
        system->setHamiltonian    (new HarmonicOscillator(system, omega, beta));
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
        system->setHamiltonian    (new HarmonicOscillator(system, omega, beta));
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

void loopinteraction() {
      cout << "Running iterative simulation with different hard shell diameters 'a'." << endl;
      // ----------------SYSTEM PARAMETERS---------------------
      // Physical system parameters
      int numberOfDimensions  = 3;
      int numberOfParticles   = 5;
      double omega            = 1.0;       // Oscillator frequency
      double alpha            = 0.5;       // Initial guess taken from non-int case
      double beta             = 1.;       // Shift of trap in the z-direction
      // Metropolis parameters
      int metroSteps, equiSteps;           // Steps are specified down with the type of calculation
      bool ImpSampling        = 1;         // cout << "Perform importance sampling? [0,1] : "; cin >> ImpSampling;
      double stepLength       = 0.1;       // Metropolis step length used without importance sampling
      double timestep         = 0.01;     // Time step used in importance sampling movement
      // Other parameters
      double a_vals[3] = {0.0 , 0.00443, 0.443};
      for (int i = 0; i < 3; i++ ) {
          int seed                = i;        // Seed for the random number generator. -1 means random seed
          int num_threads         = -1;        // Number of threads for calculation. -1 means max (automatic)
          string systemName = "Different_a_vals/a_loop_finaltry" + to_string(a_vals[i]);

          // ---------------SYSTEM SETUP-----------------------
          System* system = new System(num_threads, seed);
          system->setHamiltonian    (new HarmonicOscillator(system, omega, beta));
          system->setWaveFunction   (new Correlated(system, alpha, a_vals[i], beta));
          system->setInitialState   (new RandomUniform(system, numberOfDimensions, numberOfParticles, a_vals[i]));
          system->setStepLength     (stepLength);
          system->setTimeStep       (timestep);
          system->setChoice         (ImpSampling);

          // -------------Calculations---------------------
          ParamTester* paramTester = new ParamTester(system, systemName);

          // Doing a large scale calculation for optimal alpha
          metroSteps = (int) 1e7; // Number of metropolis steps
          equiSteps = (int) 1e6;    // Amount of the total steps used for equilibration
          system->setNumberOfSteps(metroSteps, equiSteps);
          paramTester->bigCalc(0.498944);

  }
}

void medUtenImp() {

  cout << "Running with/without comparison of importance sampling performance.";
  cout << endl;

  // ----------------SYSTEM PARAMETERS---------------------
  // Physical system parameters
  int numberOfDimensions  = 3;
  int numberOfParticles   = 5;
  double omega            = 1.0;       // Oscillator frequency
  double alpha            = 0.5;       // Initial guess taken from non-int case
  double a                = 0.;
  double beta             = 1.;       // Shift of trap in the z-direction
  // Metropolis parameters
  int metroSteps, equiSteps;           // Steps are specified down with the type of calculation
  bool ImpSampling        = 1;         // cout << "Perform importance sampling? [0,1] : "; cin >> ImpSampling;
  double stepLength       = 1.;       // Metropolis step length used without importance sampling
  double timestep         = 0.01;     // Time step used in importance sampling movement
  // Other parameters

      int seed                = 2021;        // Seed for the random number generator. -1 means random seed
      int num_threads         = -1;        // Number of threads for calculation. -1 means max (automatic)
      string systemName = "WandWOImp/a_loop" + to_string(1);

      // ---------------SYSTEM SETUP-----------------------
      System* system = new System(num_threads, seed);
      system->setHamiltonian    (new HarmonicOscillator(system, omega, beta));
      system->setWaveFunction   (new SimpleGaussian(system, alpha, a, beta));
      system->setInitialState   (new RandomUniform(system, numberOfDimensions, numberOfParticles, a));
      system->setStepLength     (stepLength);
      system->setTimeStep       (timestep);
      system->setChoice         (ImpSampling);

      // -------------Calculations---------------------
      ParamTester* paramTester = new ParamTester(system, systemName);

      // Doing a large scale calculation for optimal alpha
      metroSteps = (int) 1e5; // Number of metropolis steps
      equiSteps = (int) 1e4;    // Amount of the total steps used for equilibration
      system->setNumberOfSteps(metroSteps, equiSteps);
      paramTester->bigCalc(0.498944);


    return;
}

void medUtenParallell() {

  cout << "Running performance comparison for different numbers of threads";
  cout << endl;

  // ----------------SYSTEM PARAMETERS---------------------
  // Physical system parameters
  int numberOfDimensions  = 3;
  int numberOfParticles   = 5;
  double omega            = 1.0;       // Oscillator frequency
  double alpha            = 0.5;       // Initial guess taken from non-int case
  double a                = 0.;
  double beta             = 1.;       // Shift of trap in the z-direction
  // Metropolis parameters
  int metroSteps, equiSteps;           // Steps are specified down with the type of calculation
  bool ImpSampling        = 0;         // cout << "Perform importance sampling? [0,1] : "; cin >> ImpSampling;
  double stepLength       = 0.1;       // Metropolis step length used without importance sampling
  double timestep         = 0.01;     // Time step used in importance sampling movement
  // Other parameters

      int seed                = 2021;        // Seed for the random number generator. -1 means random seed
      int num_threads         = -1;        // Number of threads for calculation. -1 means max (automatic)
      string systemName = "WandWOPara/thread_max";

      // ---------------SYSTEM SETUP-----------------------
      System* system = new System(num_threads, seed);
      system->setHamiltonian    (new HarmonicOscillator(system, omega, beta));
      system->setWaveFunction   (new SimpleGaussian(system, alpha, a, beta));
      system->setInitialState   (new RandomUniform(system, numberOfDimensions, numberOfParticles, a));
      system->setStepLength     (stepLength);
      system->setTimeStep       (timestep);
      system->setChoice         (ImpSampling);

      // -------------Calculations---------------------
      ParamTester* paramTester = new ParamTester(system, systemName);

      // Doing a large scale calculation for optimal alpha
      metroSteps = (int) 1e7; // Number of metropolis steps
      equiSteps = (int) 1e5;    // Amount of the total steps used for equilibration
      system->setNumberOfSteps(metroSteps, equiSteps);
      paramTester->bigCalc(0.498944);


    return;
    return;
}
