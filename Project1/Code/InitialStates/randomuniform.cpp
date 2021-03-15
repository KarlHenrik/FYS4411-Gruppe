#include "randomuniform.h"
#include <cassert>
#include "Math/random.h"
#include "../particle.h"
#include "../system.h"

RandomUniform::RandomUniform(System* system, int numberOfDimensions, int numberOfParticles) : InitialState(system) {
    assert(numberOfDimensions > 0 && numberOfParticles > 0);
    m_numberOfDimensions = numberOfDimensions;
    m_numberOfParticles  = numberOfParticles;

    m_system->setNumberOfDimensions(numberOfDimensions);
    m_system->setNumberOfParticles(numberOfParticles);
}

std::vector<class Particle*> RandomUniform::newParticles() {
    Random* rng = m_system->getRandomEngine();
    std::vector<class Particle*> particles;
    for (int i=0; i < m_numberOfParticles; i++) {
        std::vector<double> position = std::vector<double>();

        for (int j=0; j < m_numberOfDimensions; j++) {
            position.push_back(rng->nextDouble());
        }
        particles.push_back(new Particle());
        particles.at(i)->setNumberOfDimensions(m_numberOfDimensions);
        particles.at(i)->setPosition(position);
    }
    return particles;
}
