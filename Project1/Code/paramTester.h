#pragma once
#include <string>
#include <fstream>

using namespace std;

class ParamTester {
public:
    ParamTester(class System* system, string fileName);
    void alphaGrid(double alpha, double alpha_end, double alpha_step);
    void printOutputToTerminal();
private:
    class System* m_system = nullptr;
    ofstream m_ofile;
};
