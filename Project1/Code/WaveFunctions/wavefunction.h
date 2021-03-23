#pragma once
#include <vector>


class WaveFunction {
public:
    WaveFunction(class System* system);
    int     getNumberOfParameters() { return m_numberOfParameters; }
    std::vector<double> getParameters() { return m_parameters; }
    virtual void setParameters(std::vector<double>);
    virtual double evaluate(std::vector<class Particle*> particles) = 0;
    virtual double computeDoubleDerivative(std::vector<class Particle*> particles) = 0;
    virtual double computeLocalDoubleDerivative(std::vector<class Particle*> particles) = 0;
    virtual double ComputeLocalFullDer(std::vector<class Particle*> particles) = 0;
    virtual double evaluateChange(Particle*, double, double) = 0;
    virtual std::vector<double> ComputeQF(Particle*, std::vector<double>) = 0;
    virtual double computeParamDer(std::vector<Particle*> particles) = 0;

protected:
    int     m_numberOfParameters = 0;
    std::vector<double> m_parameters = std::vector<double>();
    class System* m_system = nullptr;
};
