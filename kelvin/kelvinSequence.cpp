#include "kelvinSequence.hpp"
#include "interface.hpp"

KelvinSequence::KelvinSequence() : count(0) {
  shapes.resize(count);
  dt.resize(count);
  stress.resize(count);
}

KelvinSequence::KelvinSequence(uint32_t count) : count(count) {
  SetCount(count);
}

void KelvinSequence::SetCount(uint32_t _count) {
  count = _count;

  shapes.resize(count);
  dt.resize(count);
  stress.resize(count);

  for (uint32_t i = 0; i < count; i++) {
    shapes[i].SetPreviousShape(i > 0 ? &shapes[i - 1] : NULL);
    shapes[i].SetIntegrationStep(INTEGRATION_STEP);
  }
}

void KelvinSequence::SetLaplace(Laplace* laplace) {
  for (uint32_t i = 0; i < count; i++) {
    shapes[i].SetLaplace(laplace);
  }
}

void KelvinSequence::SetParameters(double p, double nu, double k, double eta, double _dt) {
  for (uint32_t i = 0; i < count; i++) {
    dt[i] = _dt;
    shapes[i].SetParameters(p, nu, k, eta, _dt);
  }
}

void KelvinSequence::SetParameters(double p, double nu, double k, double eta, std::vector<double>* _dt) {
  for (uint32_t i = 0; i < count; i++) {
    dt[i] = (*_dt)[i];
    shapes[i].SetParameters(p, nu, k, eta, dt[i]);
  }
}

void KelvinSequence::SetApexStresses(std::vector<double>* stresses) {
  double lowest = 1e10;
  double highest = 0;
  for (uint32_t i = 0; i < count; i++) {
    lowest = ((*stresses)[i] < lowest) ? (*stresses)[i] : lowest;
    highest = ((*stresses)[i] > highest) ? (*stresses)[i] : highest;

    shapes[i].SetInitialConditions(0.0, 0.0, 0.0, (*stresses)[i]);
  }

  stressAmplitude = (highest - lowest) / 2.;
}

bool KelvinSequence::Solve(bool parallel) {
  bool success = true;
  for (uint32_t i = 0; i < count; i++) {
    success &= shapes[i].PressureShooting(parallel);
    if (!success || shapes[i].IsInvalid()) {
      printf("Failed to generate the %dth shape\n", i);
      return false; 
    }
  }

  return success;
}

void KelvinSequence::Interpolate() {
  for (uint32_t i = 0; i < count; i++) {
    shapes[i].Interpolate(true, true);
  }
}

void KelvinSequence::GetPointSets(std::vector<PointSet>* _shapes) {
  _shapes->resize(count);
  for (uint32_t i = 0; i < count; i++) {
    generatePointSet(&shapes[i], &(*_shapes)[i], true);
  }
}
