#pragma once
#include <vector>

using namespace std;

class Particle {
public:
    Particle();
    void setPosition(vector<double> position);
    void adjustPosition(double change, int dimension);
    void adjustLangevin(vector<double> move);
    void setNumberOfDimensions(int numberOfDimensions);
    vector<double> getPosition() { return m_position; }
    double getLengthSq() { return m_lengthSq; }
    int getDims() { return m_numberOfDimensions; }

private:
    int m_numberOfDimensions = 0;
    vector<double> m_position = std::vector<double>();
    double m_lengthSq;
    void updateLengthSq();
};
