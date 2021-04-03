#include "paramTester.h"
#include "system.h"
#include "sampler.h"
#include "WaveFunctions/wavefunction.h"
#include "WaveFunctions/simplegaussian.h"

#include <chrono> // timing
#include <ctime>
#include <ratio>

#include <cmath> // log10
#include <string> // output
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

using namespace std;

ParamTester::ParamTester(System* system, string fileName) {
    m_system = system;
    m_fileName = fileName;
}

// alphaGD is used for finding the alpha which gives the lowest energy with gradient descent, and then doing a large calculation for that alpha
double ParamTester::alphaGD(double alpha, double lr, double tol, int max_iter) {
    m_system->addOutput(systemInfo());

    vector<double> parameters {0};
    double delta = tol + 1;
    int iter = 0;

    chrono::high_resolution_clock::time_point start = chrono::high_resolution_clock::now();
    // looping through the alphas of interest
    while (delta > tol && iter < max_iter) {
        Sampler* sampler = new Sampler(m_system, m_system->getNumerOfThreads());
        m_system->setSampler(sampler);
        parameters.at(0) = alpha;

        m_system->getWaveFunction()->setParameters(parameters);
        m_system->runMetropolisSteps(m_system->getChoice());

        cout << "\r" + to_string(alpha) + " - " + to_string((int) (1.0 * iter / max_iter * 100.0) ) + "%" << flush;

        delta = abs(sampler->getParamDer() * lr);
        alpha = alpha - sampler->getParamDer() * lr;
        iter++;
    }
    cout << endl;
    chrono::high_resolution_clock::time_point stop = chrono::high_resolution_clock::now();
    chrono::duration<double> time_span = chrono::duration_cast<chrono::duration<double>>(stop - start);
    m_system->addOutput("Execution time: " + to_string(time_span.count()) + "\n\n");

    printOutputToTerminal();
    writeOutputToFile("_GD");
    m_system->clearOutput();
    return alpha;
}


void ParamTester::bigCalc(double alpha) {
    m_system->addOutput("         Energy\n");
    vector<double> parameters {0};

    cout << "Running large scale calculation for alpha = " + to_string(alpha) << endl;

    chrono::high_resolution_clock::time_point start = chrono::high_resolution_clock::now();
    // looping through the alphas of interest
    m_system->setSampler(new SavingSampler(m_system, m_system->getNumerOfThreads(), m_fileName));
    parameters.at(0) = alpha;

    m_system->getWaveFunction()->setParameters(parameters);
    m_system->runMetropolisSteps(m_system->getChoice());

    chrono::high_resolution_clock::time_point stop = chrono::high_resolution_clock::now();
    chrono::duration<double> time_span = chrono::duration_cast<chrono::duration<double>>(stop - start);
    cout << "Execution time: " + to_string(time_span.count()) + "\n\n" << endl;

    m_system->addOutput(to_string(time_span.count()));
    writeOutputToFile("_Big");
    m_system->clearOutput();
}


// alphaGrid is used for finding finding decent expected values for a range of alphas, not for a proper calculation or search
void ParamTester::alphaGrid(double alpha, double alpha_end, double alpha_step) {
    m_system->addOutput(systemInfo());

    vector<double> parameters {0};
    double n = (alpha_end - alpha) / alpha_step;
    double i = 0;

    chrono::high_resolution_clock::time_point start = chrono::high_resolution_clock::now();
    // looping through the alphas of interest
    for (;alpha <= alpha_end + 1e-12; alpha += alpha_step) {
        m_system->setSampler(new Sampler(m_system, m_system->getNumerOfThreads()));
        parameters.at(0) = alpha;

        m_system->getWaveFunction()->setParameters(parameters);
        m_system->runMetropolisSteps(m_system->getChoice());


        cout << "\r" + to_string(alpha) + " - " + to_string((int) (i++ / n * 100) ) + "%" << flush;
    }
    cout << endl;
    chrono::high_resolution_clock::time_point stop = chrono::high_resolution_clock::now();
    chrono::duration<double> time_span = chrono::duration_cast<chrono::duration<double>>(stop - start);
    m_system->addOutput("Execution time: " + to_string(time_span.count()) + "\n\n");

    printOutputToTerminal();
    writeOutputToFile("_Grid");
    m_system->clearOutput();
}

void ParamTester::printOutputToTerminal() {
    cout << m_system->getOutput();
}

void ParamTester::writeOutputToFile(string ending) {
    ofstream ofile;
    ofile.open("Output/" + m_fileName + ending + ".txt");
    ofile << setiosflags(ios::showpoint | ios::uppercase);

    ofile << m_system->getOutput();

    ofile.close();
}

string ParamTester::systemInfo() {
    int     np = m_system->getNumberOfParticles();
    int     nd = m_system->getNumberOfDimensions();
    int     ms = m_system->getMetroSteps();
    int     p  = m_system->getWaveFunction()->getNumberOfParameters();
    int     eq = m_system->getEquiSteps();

    stringstream buffer;

    buffer << "  -- System info -- " << endl;
    buffer << " Number of particles  : " << np << endl;
    buffer << " Number of dimensions : " << nd << endl;
    buffer << " Number of Metropolis steps run : 10^" << log10(ms) << endl;
    buffer << " Number of equilibration steps  : 10^" << log10(eq) << endl;
    buffer << " Number of parameters : " << p << endl;
    buffer << endl;
    buffer << "  -- Results -- " << endl;
    buffer << setiosflags(ios::showpoint | ios::uppercase);
    for (int i = 0; i < p; i++) {
        buffer << left << " " << setw(15) << "Param_" + to_string(i + 1);
    }
    buffer << setw(15) << "E";
    buffer << setw(15) << "E^2";
    buffer << setw(15) << "VAR";
    buffer << setw(15) << "% Accepted Transitions";
    buffer << endl;

    return buffer.str();
}
