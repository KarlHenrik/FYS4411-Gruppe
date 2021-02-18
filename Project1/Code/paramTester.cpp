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
#include <chrono> // timing
#include <ctime>
#include <ratio>

#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cstdlib>

using namespace std;

ParamTester::ParamTester(System* system, string fileName) {
    m_system = system;
    m_fileName = fileName;
}

void ParamTester::alphaGrid(double alpha, double alpha_end, double alpha_step) {
    m_system->addOutput(systemInfo());

    vector<double> parameters {0};
    double n = (alpha_end - alpha) / alpha_step;
    double i = 0;

    chrono::high_resolution_clock::time_point start = chrono::high_resolution_clock::now();

    for (;alpha <= alpha_end + 1e-12; alpha += alpha_step) {
        parameters.at(0) = alpha;

        m_system->getWaveFunction()->setParameters(parameters);
        m_system->runMetropolisSteps();

        cout << "\r" + to_string(alpha) + " - " + to_string((int) (i++ / n * 100) ) + "%" << flush;
    }

    chrono::high_resolution_clock::time_point stop = chrono::high_resolution_clock::now();
    chrono::duration<double> time_span = chrono::duration_cast<chrono::duration<double>>(start - stop);

    m_system->addOutput("Execution time: " + to_string(time_span.count()) + "\n");
    cout << endl;
    printOutputToTerminal();
    writeOutputToFile();
}

void ParamTester::printOutputToTerminal() {
    cout << m_system->getOutput();
}

void ParamTester::writeOutputToFile() {
    ofstream ofile;
    ofile.open("Output/" + m_fileName);
    ofile << setiosflags(ios::showpoint | ios::uppercase);

    ofile << m_system->getOutput();

    ofile.close();
}

string ParamTester::systemInfo() {
    int     np = m_system->getNumberOfParticles();
    int     nd = m_system->getNumberOfDimensions();
    int     ms = m_system->getNumberOfMetropolisSteps();
    int     p  = m_system->getWaveFunction()->getNumberOfParameters();
    double  ef = m_system->getEquilibrationFraction();

    stringstream buffer;

    buffer << "  -- System info -- " << endl;
    buffer << " Number of particles  : " << np << endl;
    buffer << " Number of dimensions : " << nd << endl;
    buffer << " Number of Metropolis steps run : 10^" << std::log10(ms) << endl;
    buffer << " Number of equilibration steps  : 10^" << std::log10(std::round(ms*ef)) << endl;
    buffer << " Number of parameters : " << p << endl;
    buffer << endl;
    buffer << "  -- Results -- " << endl;
    buffer << setiosflags(ios::showpoint | ios::uppercase);
    for (int i = 0; i < p; i++) {
        buffer << left << " " << setw(15) << "Param " + to_string(i + 1);
    }
    buffer << setw(15) << "E";
    buffer << setw(15) << "E^2";
    buffer << setw(15) << "VAR";
    buffer << endl;

    return buffer.str();
}
