#include "simplegaussian.h"
#include <cmath>
#include <cassert>
#include "wavefunction.h"
#include "../system.h"
#include "../particle.h"

#include <iostream>

using namespace std;

SimpleGaussian::SimpleGaussian(System* system, double alpha, double /*a*/, double beta) : WaveFunction(system) {
    assert(alpha >= 0);
    m_numberOfParameters = 1;
    m_parameters.reserve(1);
    m_parameters.push_back(alpha);
    ell = {1, 1, beta};
}

double SimpleGaussian::computeRatio(vector<Particle*> /*particles*/, int /*particle_idx*/, Particle* randParticle, vector<double> oldPos, int /*thread_num*/) {
    double alpha = m_parameters[0];
    double ell_r2 = 0;
    double old_ell_r2 = 0;
    for (int d = 0; d < randParticle->getDims(); d++) {
        ell_r2 += randParticle->getPosition().at(d) * randParticle->getPosition().at(d) * ell.at(d);
        old_ell_r2 += oldPos.at(d) * oldPos.at(d) * ell.at(d);
    }
    double waveFuncRatio = exp(-alpha * ell_r2 + alpha * old_ell_r2);

    //Old wavefunc value term: exp(-alpha * old_ell_r)
    //New wavefunc value term: exp(-alpha * ell_r)
    //All other terms are the same and cancel. Since they are both exponentials, we can subtract the exponents

    return waveFuncRatio;
}

double SimpleGaussian::computeDoubleDerivative(vector<class Particle*> particles, int /*thread*/) {
    double doubleDerivative = 0;
    double alpha = m_parameters[0];
    double d_pos;
    for (int i = 0; i < (int) particles.size(); i++) {
        for (int d = 0; d < particles[i]->getDims(); d++) {
            d_pos = particles[i]->getPosition().at(d);
            doubleDerivative += ell[d] - 2 * alpha * d_pos * d_pos * ell[d] * ell[d];
        }
    }
    doubleDerivative *= -2.0 * alpha;
    return doubleDerivative;
}

// Computation of the Quantum Force with the special case of a Gaussian trial WF
vector<double> SimpleGaussian::computeQF(vector<Particle*> /*particles*/, int /*particle_idx*/, Particle* randParticle, vector<double> oldPos, int /*thread*/) {
    vector <double> QuantumForce(randParticle->getDims(), 0);
    double alpha = m_parameters[0];

    for (int i = 0; i < randParticle->getDims(); i++) {
        QuantumForce[i] = -4 * alpha * oldPos[i] * ell[i];
    }
    return QuantumForce;
}


double SimpleGaussian::computeParamDer(vector<Particle*> particles) {
    // Only works for this very specific problem!
    double der = 0;
    for (int i = 0; i < (int) particles.size(); i++) {
        vector<double> ipos = particles[i]->getPosition();
        for (int d = 0; d < particles[i]->getDims(); d++) {
            der += ipos[d] * ipos[d] * ell[d];
        }
    }
    return -der / particles.size();
}
