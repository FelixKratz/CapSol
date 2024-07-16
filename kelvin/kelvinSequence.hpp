#pragma once

#include "kelvin.hpp"
#include "pressureShooting.hpp"

class KelvinSequence {
  public:
    double stressAmplitude = 0.;
    std::vector<Kelvin> shapes;
    std::vector<double> dt;
    std::vector<double> stress;
    uint32_t count = 0;

    KelvinSequence();
    KelvinSequence(uint32_t count);

    void SetCount(uint32_t _count);
    void SetLaplace(Laplace* laplace);
    void SetApexStresses(std::vector<double>* stresses);
    void SetParameters(double p, double nu, double k, double eta, double dt);
    void SetParameters(double p, double nu, double k, double eta, std::vector<double>* _dt);
    bool Solve(bool parallel = false);
    void SaveToDisk(std::string path);
    void Interpolate();
    void GetPointSets(std::vector<PointSet>* _shapes);
};
