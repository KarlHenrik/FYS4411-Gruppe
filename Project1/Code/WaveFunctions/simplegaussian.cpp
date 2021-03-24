#include "simplegaussian.h"
#include <cmath>
#include <cassert>
#include "wavefunction.h"
#include "../system.h"
#include "../particle.h"

#include <iostream>

using namespace std;

SimpleGaussian::SimpleGaussian(System* system, double alpha) :
        WaveFunction(system) {
    assert(alpha >= 0);
    m_numberOfParameters = 1;
    m_parameters.reserve(1);
    m_parameters.push_back(alpha);
}

double SimpleGaussian::computeRatio(vector<Particle*> particles, int particle_idx, Particle* randParticle, vector<double> oldPos, double oldLengthSq, int thread_num) {
    (void)particle_idx, (void)oldPos, (void)thread_num, (void)particles;
    double alpha = m_parameters[0];
    double waveFuncRatio = exp(-alpha * randParticle->getLengthSq() + alpha * oldLengthSq);

    //Old wavefunc value term: exp(-alpha * oldLengthSq)
    //New wavefunc value term: exp(-alpha * randParticle->getLengthSq())
    //All other terms are the same and cancel. Since they are both exponentials, we can subtract the exponents

    return waveFuncRatio;
}

double SimpleGaussian::computeDoubleDerivative(vector<class Particle*> particles) {
    double doubleDerivative = 0;
    double alpha = m_parameters[0];

    for (int i = 0; i < (int) particles.size(); i++) {
        doubleDerivative += particles[i]->getDims() - 2.0 * alpha * particles[i]->getLengthSq();
    }
    doubleDerivative *= -2.0 * alpha;

    return doubleDerivative;
}

// Computation of the Quantum Force with the special case of a Gaussian trial WF
vector<double> SimpleGaussian::computeQF(vector<Particle*> particles, int particle_idx, Particle* randParticle, vector<double> oldPos, int thread) {
    (void)particle_idx, (void)particles, (void)thread;
    vector <double> QuantumForce(randParticle->getDims(), 0);
    double alpha = m_parameters[0];

    for (int i = 0; i < randParticle->getDims(); i++) {
        QuantumForce[i] = -4 * alpha * oldPos[i];
    }
    return QuantumForce;
}


double SimpleGaussian::computeParamDer(vector<Particle*> particles) {
    // Only works for this very specific problem!
    double der = 0;
    for (int i = 0; i < (int) particles.size(); i++) {
        der += particles[i]->getLengthSq();
    }
    return der;
}

void SimpleGaussian::setup(vector<Particle*> particles, int thread) {
    (void)particles, (void)thread;
    return;
}
