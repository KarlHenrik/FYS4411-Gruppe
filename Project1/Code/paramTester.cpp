#include "paramTester.h"
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

#include <string>
#include <fstream>
#include <iomanip>
#include <cstdlib>

using namespace std;

ParamTester::ParamTester(System* system, string fileName) {
    m_ofile.open("Output/" + fileName);
    m_ofile << setiosflags(ios::showpoint | ios::uppercase);

    m_system = system;
}

void ParamTester::alphaGrid(double alpha, double alpha_end, double alpha_step) {
    printOutputToTerminal();

    vector<double> parameters {0};
    for (;alpha <= alpha_end + 1e-12; alpha += alpha_step) {
        parameters.at(0) = alpha;

        m_system->getWaveFunction()->setParameters(parameters);
        m_system->runMetropolisSteps();
    }
}

void ParamTester::printOutputToTerminal() {
    int     np = m_system->getNumberOfParticles();
    int     nd = m_system->getNumberOfDimensions();
    int     ms = m_system->getNumberOfMetropolisSteps();
    int     p  = m_system->getWaveFunction()->getNumberOfParameters();
    double  ef = m_system->getEquilibrationFraction();

    cout << endl;
    cout << "  -- System info -- " << endl;
    cout << " Number of particles  : " << np << endl;
    cout << " Number of dimensions : " << nd << endl;
    cout << " Number of Metropolis steps run : 10^" << std::log10(ms) << endl;
    cout << " Number of equilibration steps  : 10^" << std::log10(std::round(ms*ef)) << endl;
    cout << " Number of parameters : " << p << endl;
    cout << endl;
    cout << "  -- Results -- " << endl;
    cout << setiosflags(ios::showpoint | ios::uppercase);
    for (int i = 0; i < p; i++) {
        cout << left << " " << setw(15) << "Param " + to_string(i + 1);
    }
    cout << setw(15) << "E";
    cout << setw(15) << "E^2";
    cout << endl;

}
