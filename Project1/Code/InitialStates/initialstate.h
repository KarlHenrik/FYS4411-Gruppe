#pragma once
#include <vector>
#include "../particle.h"

class InitialState {
public:
    InitialState(class System* system);
    virtual std::vector<class Particle*> newParticles() = 0;

protected:
    class System* m_system = nullptr;
    int m_numberOfDimensions = 0;
    int m_numberOfParticles = 0;
};
