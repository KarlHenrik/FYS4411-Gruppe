#pragma once
#include <vector>

class Particle {
public:
    Particle();
    void setPosition(const std::vector<double> &position);
    void adjustPosition(double change, int dimension);
    void adjustLangevin(std::vector<double> move);
    void setNumberOfDimensions(int numberOfDimensions);
    std::vector<double> getPosition() { return m_position; }
    double getLengthSq() { return m_lengthSq; }
    int getDims() { return m_numberOfDimensions; }

private:
    int m_numberOfDimensions = 0;
    std::vector<double> m_position = std::vector<double>();
    double m_lengthSq;
    void updateLengthSq();
};
