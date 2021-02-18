#pragma once

#include <string>

using namespace std;

class ParamTester {
public:
    ParamTester(class System* system, string fileName);
    void alphaGrid(double alpha, double alpha_end, double alpha_step);
    void printOutputToTerminal();
    void writeOutputToFile();
    string systemInfo();
private:
    class System* m_system = nullptr;
    string m_fileName;
};
