#include "randomuniform.h"
#include <cassert>
#include "Math/random.h"
#include "../particle.h"
#include "../system.h"
#include <cmath>

#include <iostream>

using namespace std;

RandomUniform::RandomUniform(System* system, int numberOfDimensions, int numberOfParticles, double a) : InitialState(system) {
    assert(numberOfDimensions > 0 && numberOfParticles > 0);
    m_numberOfDimensions = numberOfDimensions;
    m_numberOfParticles  = numberOfParticles;
    m_a = a;
    if (m_a > 0.5) {
        cout << "a is too large! You will not be able to place the particles with the current scheme" << endl;
    }

    m_system->setNumberOfDimensions(numberOfDimensions);
    m_system->setNumberOfParticles(numberOfParticles);
}

vector<class Particle*> RandomUniform::newParticles() {
    Random* rng = m_system->getRandomEngine();
    vector<class Particle*> particles;
    for (int i=0; i < m_numberOfParticles; i++) {
        vector<double> position = findPos(particles, rng);

        particles.push_back(new Particle());
        particles.at(i)->setNumberOfDimensions(m_numberOfDimensions);
        particles.at(i)->setPosition(position);
    }
    return particles;
}

vector<double> RandomUniform::findPos(vector<class Particle*> particles, Random* rng) {
    vector<double> position = vector<double>();
    for (int j=0; j < m_numberOfDimensions; j++) {
        position.push_back((rng->nextDouble() - 0.5) * 2.0); // uniform from -1 to 1
    }

    double limit = m_a * m_a;
    double min_dist = m_a + 1;
    for (int k = 0; k < (int) particles.size(); k++) {
        double dist = 0;
        for (int d=0; d < m_numberOfDimensions; d++) {
            dist += pow(position.at(d) - particles.at(k)->getPosition().at(d), 2);
        }
        if (dist < min_dist) {
            min_dist = dist;
        }
    }
    if (min_dist > limit) {
        return position;
    } else {
        return findPos(particles, rng);
    }
}
