#pragma once
#include "wavefunction.h"

class SimpleGaussian : public WaveFunction {
public:
    SimpleGaussian(class System* system, double alpha);
    double evaluate(std::vector<class Particle*> particles);
    double computeDoubleDerivative(std::vector<class Particle*> particles);
    double ComputeLocalFullDer(std::vector<class Particle*> particles);
    std::vector<double> ComputeQF(Particle*, std::vector<double>);
    double evaluateChange(Particle*, double, double);
    double computeParamDer(std::vector<Particle*> particles);
    double **AllocateMatrix(int,int);
    double **MakeDistanceMatrix(std::vector<Particle*> particles);
};
