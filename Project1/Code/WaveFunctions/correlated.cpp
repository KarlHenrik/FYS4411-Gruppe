#include "correlated.h"
#include <cmath>
#include <cassert>
#include "wavefunction.h"
#include "../system.h"
#include "../particle.h"

#include <iostream>

using namespace std;

Correlated::Correlated(System* system, double alpha) :
        WaveFunction(system) {
    assert(alpha >= 0);
    m_numberOfParameters = 1;
    m_parameters.reserve(1);
    m_parameters.push_back(alpha);

    // initializing matrix for storing relative distances between bosons
    old_Dist = new double **[m_system->getNumerOfThreads()];
    new_Dist = new double **[m_system->getNumerOfThreads()];
    unit_Vec = new double ***[m_system->getNumerOfThreads()];
}

double Correlated::computeRatio(vector<Particle*> particles, int particle_idx, Particle* randParticle, vector<double> oldPos, double oldLengthSq, int thread_num) {
    double alpha = m_parameters[0];

    double** old_R = old_Dist[thread_num];
    double** new_R = new_Dist[thread_num];
    // Starting with the uncorrelated part
    double oldTerm = exp(-alpha * oldLengthSq);
    double newTerm = exp(-alpha * particles[particle_idx]->getLengthSq());
    for (int i = 0; i < particle_idx; i++) { // particles before randParticle
        oldTerm *= 1 - m_a / old_R[particle_idx][i];
        newTerm *= 1 - m_a / new_R[particle_idx][i];
    }
    for (int i = particle_idx + 1; i < (int) particles.size(); i++) { // particles after randParticle
        oldTerm *= 1 - m_a / old_R[particle_idx][i];
        newTerm *= 1 - m_a / new_R[particle_idx][i];
    }
    double waveFuncRatio = newTerm / oldTerm;
    //All other terms are the same and cancel

    return waveFuncRatio;
}

double Correlated::computeDoubleDerivative(vector<class Particle*> particles) {
    double doubleDerivative = 0;
    double alpha = m_parameters[0];

    for (int i = 0; i < (int) particles.size(); i++) {
        doubleDerivative += particles[i]->getDims() - 2.0 * alpha * particles[i]->getLengthSq();
    }
    doubleDerivative *= -2.0 * alpha;

    return doubleDerivative;
}

// Computation of the Quantum Force with the special case of a Gaussian trial WF
vector <double> Correlated::computeQF(vector<Particle*> particles, int particle_idx, Particle* randParticle, vector<double> oldPos, int thread) {
    vector<double> QuantumForce(randParticle->getDims(), 0);
    double alpha = m_parameters[0];

    double temp_fac = 0;
    vector<double> vecSum(randParticle->getDims(), 0);
    for (int i = 0; i < particle_idx; i++) {
        temp_fac = m_a / (pow(old_Dist[particle_idx][i], 2) - m_a * old_Dist[particle_idx][i])
        for (int d = 0; d < (int) particles[0]->getDims(); d++) {
            vecSum[d] += unit_Vec[thread][particle_idx][i][d] * temp_fac;
        }
    } for (int i = particle_idx + 1; i < (int) particles.size(); i++) {
        temp_fac = m_a / (pow(old_Dist[particle_idx][i], 2) - m_a * old_Dist[particle_idx][i])
        for (int d = 0; d < (int) particles[0]->getDims(); d++) {
            vecSum[d] += unit_Vec[thread][particle_idx][i][d] * temp_fac;
        }
    }

    for (int d = 0; d < (int) randParticle->getDims(); d++) {
        QuantumForce[d] = -4 * alpha * oldPos[d] + 2 * vecSum[d];
    }
    return QuantumForce;
}


double Correlated::computeParamDer(vector<Particle*> particles) {
    // Only works for this very specific problem!
    double der = 0;
    for (int i = 0; i < (int) particles.size(); i++) {
        der += particles[i]->getLengthSq();
    }
    return der;
}

void Correlated::setup(vector<Particle*> particles, int thread) {
    old_Dist[thread] = new double *[particles.size()];
    new_Dist[thread] = new double *[particles.size()];
    unit_Vec[thread] = new double **[particles.size()];
    for (int p1 = 0; p1 < (int) particles.size(); p1++) {
        old_Dist[thread][p1] = new double [particles.size()];
        new_Dist[thread][p1] = new double [particles.size()];
        unit_Vec[thread][p1] = new double *[particles.size()];
        for (int p2 = 0; p2 < (int) particles.size(); p2++) {
            double interdistance = 0;
            for (int d = 0; d < (int) particles[0]->getDims(); d++) {
                interdistance += pow(particles[p1]->getPosition().at(d) - particles[p2]->getPosition().at(d), 2);
            }
            old_Dist[thread][p1][p2] = interdistance;
            new_Dist[thread][p1][p2] = interdistance;

            unit_Vec[thread][p1][p2] = new double [particles[0]->getDims()];
            for (int d = 0; d < (int) particles[0]->getDims(); d++) {
                unit_Vec[thread][p1][p2][d] = (particles[p2]->getPosition().at(d) - particles[p1]->getPosition().at(d)) / interdistance;
            }
        }
    }
}
/*
double Correlated::ComputeLocalFullDer(vector<class Particle*> particles) {
    double alpha = m_parameters[0];
    double beta = 2;
    double a = 1;
    double DoubleDer = 0;
    double TotalDoubleDer = 0;

    //int np = particles.size();
    //int dims = particles[0]->getDims();

    double **R = MakeDistanceMatrix(particles);

    vector<double> SUM1(particles.size(),0);
    vector<double> SUM2(particles.size(),0);
    double SUM3 = 0;

    for (int p = 0; p < particles.size(); p++) {
        for (int i = 0; i < particles.size(); i++) {
            for (int j = 0; i < particles.size(); j++) {
                for (int d = 0; d < particles[p]->getDims(); d++) {
                    double rval = R[p][i];
                    double rval_j = R[p][j];
                    SUM1[d] += a / (rval * (rval*rval - a*rval))*(particles[p]->
                    getPosition()[d] - particles[i]->getPosition()[d]);
                    SUM2[d] += a / (rval_j*(rval_j*rval_j - a*rval_j))*(particles[p]->
                    getPosition()[d] - particles[j]->getPosition()[d]);
                }
            }
            SUM3 += (a*a - 2*a*R[p][i])/pow(R[p][i]*R[p][i] - a*R[p][i],2) +
            a/(R[p][i]*R[p][i]*R[p][i] - a*R[p][i]*R[p][i]);

        }
        for (int d = 0; d < particles[p]->getDims(); d++) {
            DoubleDer += -4*alpha*particles[p]->getPosition()[d]*SUM1[d] + SUM1[d]*SUM2[d];
        }
        DoubleDer += -4*alpha - 2*alpha*beta + 4*alpha*alpha*particles[p]->getLengthSq() + SUM3;
    }
    TotalDoubleDer += DoubleDer;

    return TotalDoubleDer;
}
*/
