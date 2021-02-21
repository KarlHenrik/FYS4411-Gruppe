#include "simplegaussian.h"
#include <cmath>
#include <cassert>
#include "wavefunction.h"
#include "../system.h"
#include "../particle.h"

using namespace std;

SimpleGaussian::SimpleGaussian(System* system, double alpha) :
        WaveFunction(system) {
    assert(alpha >= 0);
    m_numberOfParameters = 1;
    m_parameters.reserve(1);
    m_parameters.push_back(alpha);
}

double SimpleGaussian::evaluate(std::vector<class Particle*> particles) {
    double wf = 1;
    double alpha = m_parameters[0];
    for (int i = 0; i < particles.size(); i++) {
        wf *= exp(-alpha * particles[i]->getLengthSq());
    }

    return wf;
}

double SimpleGaussian::evaluateChange(Particle* randParticle, double waveFunctionValue, double oldLengthSq) {
    double alpha = m_parameters[0];
    waveFunctionValue /= exp(-alpha * oldLengthSq);
    waveFunctionValue *= exp(-alpha * randParticle->getLengthSq());

    return waveFunctionValue;
}

double SimpleGaussian::computeDoubleDerivative(std::vector<class Particle*> particles) {
    return computeLocalDoubleDerivative(particles) * evaluate(particles);
}

double SimpleGaussian::computeLocalDoubleDerivative(std::vector<class Particle*> particles) {
    double doubleDerivative = 0;
    double alpha = m_parameters[0];

    for (int i = 0; i < particles.size(); i++) {
        doubleDerivative += particles[i]->getDims() - 2.0 * alpha * particles[i]->getLengthSq();
    }
    doubleDerivative *= -2.0 * alpha;

    return doubleDerivative;
}

// Computation of the Quantum Force with the special case of a
// Gaussian trial WF
std::vector <double> SimpleGaussian::ComputeQF(Particle* randParticle, std::vector<double> oldPos) {
  std::vector <double> QuantumForce;
  double alpha = m_parameters[0];

  QuantumForce = oldPos;
  for (int i = 0; i < randParticle->getDims(); i++) {
    QuantumForce[i] *= -4*alpha;
  }

  return QuantumForce;
}
