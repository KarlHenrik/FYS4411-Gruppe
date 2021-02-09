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
    /* You need to implement a Gaussian wave function here. The positions of
     * the particles are accessible through the particle[i].getPosition()
     * function.
     *
     * For the actual expression, use exp(-alpha * r^2), with alpha being the
     * (only) variational parameter.
     */

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
    /* All wave functions need to implement this function, so you need to
     * find the double derivative analytically. Note that by double derivative,
     * we actually mean the sum of the Laplacians with respect to the
     * coordinates of each particle.
     *
     * This quantity is needed to compute the (local) energy (consider the
     * Schrödinger equation to see how the two are related).
     */

    double doubleDerivative = 0;
    double alpha = m_parameters[0];

    for (int i = 0; i < particles.size(); i++) {
        doubleDerivative += particles[i]->getDims() - 2.0 * alpha * particles[i]->getLengthSq();
    }
    doubleDerivative *= -2.0 * alpha;

    return doubleDerivative;
}
