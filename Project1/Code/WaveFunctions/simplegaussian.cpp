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

// Computation of the Quantum Force with the special case of a Gaussian trial WF
std::vector <double> SimpleGaussian::ComputeQF(Particle* randParticle, std::vector<double> oldPos) {
    std::vector <double> QuantumForce(randParticle->getDims(), 0);
    double alpha = m_parameters[0];

    for (int i = 0; i < randParticle->getDims(); i++) {
        QuantumForce[i] = -4 * alpha * oldPos[i];
    }

    return QuantumForce;
}


double SimpleGaussian::computeParamDer(std::vector<Particle*> particles) {
    // Only works for this very specific problem!
    double der = 0;
    for (int i = 0; i < (int) particles.size(); i++) {
        der += particles[i]->getLengthSq();
    }
    return der;
}

// Code taken from Morten H. Jensen's GITHUB
/*
double ** SimpleGaussian::AllocateMatrix(int m, int n) {
double ** Matrix;
Matrix = new double*[m];
for(int i=0;i<m;i++){
  Matrix[i] = new double[n];
  for(int j=0;j<m;j++)
    Matrix[i][j] = 0.0;
}
return Matrix;
}
*/
double ** SimpleGaussian::MakeDistanceMatrix(std::vector<Particle*> particles) {

  //int np = particles.size();

  // initializing matrix for storing relative distances between bosons
  double **R;
  R = new double*[particles.size()+2];
  for (int i = 0; i < particles.size(); i++) {
    R[i] = new double [particles.size()+2];
    for (int j = 0; j < particles.size(); j++) {
      R[j][j] = 0.;
      for (int k = j + 1; j < particles.size(); j++) {
        double interdistance = 0.;
        for (int d = 0; d < particles[i]->getDims(); d++) {
          interdistance += particles[i]->getPosition()[k] - particles[j]->
          getPosition()[k];
          R[j][k] = sqrt(interdistance*interdistance);
          R[k][j] = R[j][k];
        }
      }
    }
  }
  return R;
}

double SimpleGaussian::ComputeLocalFullDer(std::vector<class Particle*> particles) {
  double alpha = m_parameters[0];
  double beta = 2.;
  double a = 1.;
  double DoubleDer = 0.;
  double TotalDoubleDer = 0.;

  //int np = particles.size();
  //int dims = particles[0]->getDims();

  double **R = MakeDistanceMatrix(particles);

  std::vector<double> SUM1(particles.size(),0);
  std::vector<double> SUM2(particles.size(),0);
  double SUM3 = 0;

  for (int p = 0; p < particles.size(); p++) {
    for (int i = 0; i < particles.size(); i++) {
      for (int j = 0; i < particles.size(); j++) {
        for (int d = 0; d < particles[p]->getDims(); d++) {
          double rval = R[p][i];
          double rval_j = R[p][j];
          SUM1[d] += a/(rval * (rval*rval - a*rval))*(particles[p]->
            getPosition()[d] - particles[i]->getPosition()[d]);
          SUM2[d] += a/(rval_j*(rval_j*rval_j - a*rval_j))*(particles[p]->
            getPosition()[d] - particles[j]->getPosition()[d]);
        }
      }
      SUM3 += (a*a - 2*a*R[p][i])/pow(R[p][i]*R[p][i] - a*R[p][i],2) +
      a/(R[p][i]*R[p][i]*R[p][i] - a*R[p][i]*R[p][i]);

    }
    for (int d = 0; d < particles[p]->getDims(); d++) {
      DoubleDer += -4*alpha*particles[p]->getPosition()[d]*SUM1[d] + SUM1[d]*SUM2[d];
    }
    DoubleDer += -4*alpha - 2*alpha*beta + 4*alpha*alpha*particles[p]->getLengthSq()
    + SUM3;
  }
  TotalDoubleDer += DoubleDer;

  return TotalDoubleDer;
}
