#pragma once
#include <vector>

class Hamiltonian {
public:
    Hamiltonian(class System* system);
    virtual double computeEnergy(std::vector<class Particle*> particles) = 0;
    double computeKinetic(std::vector<Particle*> particles);

protected:
    class System* m_system = nullptr;
};
