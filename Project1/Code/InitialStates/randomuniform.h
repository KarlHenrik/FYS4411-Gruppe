#pragma once
#include "initialstate.h"
#include "../particle.h"
#include "Math/random.h"

using namespace std;

class RandomUniform : public InitialState {
public:
    RandomUniform(System* system, int numberOfDimensions, int numberOfParticles, double a);
    vector<class Particle*> newParticles();
    vector<double> findPos(vector<class Particle*> particles, Random* rng);
private:
    double m_a;
};
