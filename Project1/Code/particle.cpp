#include "particle.h"
#include <cassert>

Particle::Particle() {
}

void Particle::setPosition(const std::vector<double> &position) {
    assert(position.size() == (unsigned int) m_numberOfDimensions);
    m_position = position;
}

void Particle::adjustPosition(double change, int dimension) {
    m_position.at(dimension) += change;
}

void Particle::setNumberOfDimensions(int numberOfDimensions) {
    m_numberOfDimensions = numberOfDimensions;
}

double Particle::getLengthSq() {
    double lengthSq = 0;
    for (int i = 0; i < m_numberOfDimensions; i++) {
        lengthSq += m_position[i] * m_position[i];
    }
    return lengthSq;
}

int Particle::getDims() {
    return m_numberOfDimensions;
}
