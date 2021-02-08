#pragma once
#include "initialstate.h"
#include "../particle.h"

class RandomUniform : public InitialState {
public:
    RandomUniform(System* system, int numberOfDimensions, int numberOfParticles);
    std::vector<class Particle*> newParticles();
};
