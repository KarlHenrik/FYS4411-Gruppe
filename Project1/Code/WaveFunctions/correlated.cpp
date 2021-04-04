#include "correlated.h"
#include <cmath>
#include <cassert>
#include "wavefunction.h"
#include "../system.h"
#include "../particle.h"

#include <iostream>

using namespace std;

Correlated::Correlated(System* system, double alpha, double a, double beta) : WaveFunction(system) {
    assert(alpha >= 0);
    m_numberOfParameters = 1;
    m_parameters.reserve(1);
    m_parameters.push_back(alpha);
    m_a = a;
    ell = {1, 1, beta};

    // initializing matrix for storing relative distances between bosons
    old_dists = new double **[m_system->getNumerOfThreads()];
    dists = new double **[m_system->getNumerOfThreads()];
    unit_vecs = new double ***[m_system->getNumerOfThreads()];
    old_unit_vecs = new double ***[m_system->getNumerOfThreads()];
}

double Correlated::computeRatio(vector<Particle*> particles, int particle_idx, Particle* randParticle, vector<double> oldPos, int thread_num) {
    computeDists(particles, particle_idx, randParticle, thread_num);

    double alpha = m_parameters[0];
    // Starting with the uncorrelated part
    double ell_r2 = 0;
    double old_ell_r2 = 0;
    for (int d = 0; d < randParticle->getDims(); d++) {
        ell_r2 += randParticle->getPosition().at(d) * randParticle->getPosition().at(d) * ell.at(d);
        old_ell_r2 += oldPos.at(d) * oldPos.at(d) * ell.at(d);
    }
    double newTerm = exp(-alpha * ell_r2);
    double oldTerm = exp(-alpha * old_ell_r2);
    // Then the uncorrelated part which comes from the distances to the moved particle
    for (int p2 = 0; p2 < (int) particles.size(); p2++) { // particles before randParticle
        if (p2 == particle_idx) {
            continue;
        }
        if (dists[thread_num][particle_idx][p2] < m_a) {
            return 0; // Deny all moves which put particles closer than a
        }
        oldTerm *= 1 - m_a / old_dists[thread_num][particle_idx][p2];
        newTerm *= 1 - m_a / dists[thread_num][particle_idx][p2];
    }

    return newTerm / oldTerm; //All other terms are the same and cancel
}

void Correlated::computeDists(vector<Particle*> particles, int particle_idx, Particle* randParticle, int thread) {
    int p1 = particle_idx;
    vector<double> p1Pos = randParticle->getPosition();
    double interdistance;
    for (int p2 = 0; p2 < (int) particles.size(); p2++) {
        if (p2 == particle_idx) {
            continue;
        }
        interdistance = 0;
        for (int d = 0; d < (int) particles[0]->getDims(); d++) {
            interdistance += pow(p1Pos.at(d) - particles[p2]->getPosition().at(d), 2);
        }
        interdistance = sqrt(interdistance);

        old_dists[thread][p1][p2] = dists[thread][p1][p2];
        old_dists[thread][p2][p1] = dists[thread][p1][p2];
        dists[thread][p1][p2] = interdistance;
        dists[thread][p2][p1] = interdistance;

        for (int d = 0; d < (int) particles[0]->getDims(); d++) {
            old_unit_vecs[thread][p1][p2][d] = unit_vecs[thread][p1][p2][d];
            old_unit_vecs[thread][p2][p1][d] = unit_vecs[thread][p2][p1][d];
            unit_vecs[thread][p1][p2][d] = (particles[p2]->getPosition().at(d) - p1Pos.at(d)) / interdistance;
            unit_vecs[thread][p2][p1][d] = -unit_vecs[thread][p1][p2][d];
        }
    }
    dists[thread][p1][p1] = 1; // used to avoid zerodivision elsewhere
}

void Correlated::revertDists(vector<Particle*> particles, int particle_idx, int thread) {
    int p1 = particle_idx;
    for (int p2 = 0; p2 < (int) particles.size(); p2++) {
        if (p2 == particle_idx) {
            continue;
        }
        dists[thread][p1][p2] = old_dists[thread][p1][p2];
        dists[thread][p2][p1] = old_dists[thread][p2][p1];
        for (int d = 0; d < (int) particles[0]->getDims(); d++) {
            unit_vecs[thread][p1][p2][d] = old_unit_vecs[thread][p1][p2][d];
            unit_vecs[thread][p2][p1][d] = old_unit_vecs[thread][p2][p1][d];
        }
    }
    return;
}

// Computation of the Quantum Force with the special case of a Gaussian trial WF
vector <double> Correlated::computeQF(vector<Particle*> particles, int particle_idx, Particle* randParticle, vector<double> oldPos, int thread) {
    vector<double> QuantumForce(randParticle->getDims(), 0);
    double alpha = m_parameters[0];

    double temp_fac = 0;
    vector<double> vecSum(randParticle->getDims(), 0);
    for (int i = 0; i < (int) particles.size(); i++) {
        if (i == particle_idx) {
            continue;
        }
        temp_fac = m_a / (pow(old_dists[thread][particle_idx][i], 2) - m_a * old_dists[thread][particle_idx][i]);
        for (int d = 0; d < (int) particles[0]->getDims(); d++) {
            vecSum[d] += unit_vecs[thread][particle_idx][i][d] * temp_fac;
        }
    }

    for (int d = 0; d < (int) randParticle->getDims(); d++) {
        QuantumForce[d] = -4 * alpha * oldPos[d] * ell[d] + 2 * vecSum[d];
    }
    return QuantumForce;
}


double Correlated::computeParamDer(vector<Particle*> particles) {
    double der = 0;
    for (int i = 0; i < (int) particles.size(); i++) {
        vector<double> ipos = particles[i]->getPosition();
        for (int d = 0; d < particles[i]->getDims(); d++) {
            der += ipos[d] * ipos[d] * ell[d];
        }
    }
    return -der / particles.size();
}

void Correlated::clear(vector<Particle*> particles, int thread) {
    for (int p1 = 0; p1 < (int) particles.size(); p1++) {
        for (int p2 = 0; p2 < (int) particles.size(); p2++) {
            delete[] unit_vecs[thread][p1][p2];
            delete[] old_unit_vecs[thread][p1][p2];
        }
    }
    for (int p1 = 0; p1 < (int) particles.size(); p1++) {
        delete[] old_dists[thread][p1];
        delete[] dists[thread][p1];
        delete[] unit_vecs[thread][p1];
        delete[] old_unit_vecs[thread][p1];
    }
    delete[] old_dists[thread];
    delete[] dists[thread];
    delete[] unit_vecs[thread];
    delete[] old_unit_vecs[thread];
}

void Correlated::setup(vector<Particle*> particles, int thread) {
    old_dists[thread] = new double *[particles.size()];
    dists[thread] = new double *[particles.size()];
    unit_vecs[thread] = new double **[particles.size()];
    old_unit_vecs[thread] = new double **[particles.size()];
    for (int p1 = 0; p1 < (int) particles.size(); p1++) {
        old_dists[thread][p1] = new double [particles.size()];
        dists[thread][p1] = new double [particles.size()];
        unit_vecs[thread][p1] = new double *[particles.size()];
        old_unit_vecs[thread][p1] = new double *[particles.size()];
        for (int p2 = 0; p2 < (int) particles.size(); p2++) {
            double interdistance = 0;
            for (int d = 0; d < (int) particles[0]->getDims(); d++) {
                interdistance += pow(particles[p1]->getPosition().at(d) - particles[p2]->getPosition().at(d), 2);
            }
            interdistance = sqrt(interdistance);
            old_dists[thread][p1][p2] = interdistance;
            dists[thread][p1][p2] = interdistance;

            unit_vecs[thread][p1][p2] = new double [particles[0]->getDims()];
            old_unit_vecs[thread][p1][p2] = new double [particles[0]->getDims()];
            if (p1 == p2) {
                for (int d = 0; d < (int) particles[0]->getDims(); d++) {
                    unit_vecs[thread][p1][p2][d] = 0;
                    old_unit_vecs[thread][p1][p2][d] = 0;
                }
            } else {
                for (int d = 0; d < (int) particles[0]->getDims(); d++) {
                    unit_vecs[thread][p1][p2][d] = (particles[p2]->getPosition().at(d) - particles[p1]->getPosition().at(d)) / interdistance;
                    old_unit_vecs[thread][p1][p2][d] = unit_vecs[thread][p1][p2][d];
                }
            }
        }
        old_dists[thread][p1][p1] = 1; // these should never appear in a non-zero expression!
        dists[thread][p1][p1] = 1; // the unit vectors should zero these out in all cases
    }
}

double Correlated::computeDoubleDerivative(vector<class Particle*> particles, int thread) {
    double dblDer = 0;
    double alpha = m_parameters[0];
    double a = m_a;
    int dims = particles[0]->getDims();
    double temp_fac;
    double r2;
    double d_pos;

    for (int k = 0; k < (int) particles.size(); k++) {
        // Term from uncorrelated part
        for (int d = 0; d < particles[k]->getDims(); d++) {
            d_pos = particles[k]->getPosition().at(d);
            dblDer += 4 * alpha * alpha * d_pos * d_pos * ell[d] * ell[d];
            dblDer += -2 * alpha * ell[d];
        }
        // Calculating the vector sum
        vector<double> vecSum(dims, 0.0);
        for (int p2 = 0; p2 < (int) particles.size(); p2++) {
            temp_fac = a / (pow(dists[thread][k][p2], 2) - a * dists[thread][k][p2]);
            for (int d = 0; d < dims; d++) {
                vecSum.at(d) += temp_fac * unit_vecs[thread][k][p2][d];
            }
        }
        // Adding the vector products
        for (int d = 0; d < dims; d++) {
            dblDer += 4 * alpha * particles[k]->getPosition().at(d) * ell.at(d) * vecSum.at(d);
            dblDer += vecSum.at(d) * vecSum.at(d);
        }
        // Adding the final sum
        for (int p2 = 0; p2 < k; p2++) {
            r2 = dists[thread][k][p2] * dists[thread][k][p2];
            dblDer += a * (a - 2 * dists[thread][k][p2]) / pow(a * dists[thread][k][p2] - r2, 2);
            dblDer += 2 / dists[thread][k][p2] * (-a) / (a * dists[thread][k][p2] - r2);
        } for (int p2 = k + 1; p2 < (int) particles.size(); p2++) {
            r2 = dists[thread][k][p2] * dists[thread][k][p2];
            dblDer += a * (a - 2 * dists[thread][k][p2]) / pow(a * dists[thread][k][p2] - r2, 2);
            dblDer += 2 / dists[thread][k][p2] * (-a) / (a * dists[thread][k][p2] - r2);
        }
    }
    return dblDer;
}
