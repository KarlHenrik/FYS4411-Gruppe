#include "particle.h"
#include <cassert>
#include <iostream>

using namespace std;

Particle::Particle() {
}

void Particle::setPosition(vector<double> position) {
    m_position = position;
    updateLengthSq();
}

void Particle::adjustPosition(double change, int dimension) { // adjust position in one dimension
    m_position.at(dimension) += change;
    updateLengthSq();
}

void Particle::adjustLangevin(vector<double> move) { // adjust position by a vector
    for (int i = 0; i < m_numberOfDimensions; i++) {
        m_position[i] += move[i];
    }
    updateLengthSq();
}

void Particle::setNumberOfDimensions(int numberOfDimensions) {
    m_numberOfDimensions = numberOfDimensions;
}

void Particle::updateLengthSq() {
    m_lengthSq = 0;
    for (int i = 0; i < m_numberOfDimensions; i++) {
        m_lengthSq += m_position[i] * m_position[i];
    }
}
