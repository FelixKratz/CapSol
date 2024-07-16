#include "kelvin.hpp"
#include <math.h>

void Kelvin::Clone(Kelvin* _clone) {
  SetLaplace(_clone->undeformed);
  SetParameters(_clone->parameters.p, _clone->parameters.nu,
                _clone->parameters.k, _clone->parameters_kelvin.eta,
                _clone->parameters_kelvin.dt);

  SetPreviousShape(_clone->previous);
  SetIntegrationStep(_clone->integration_step);
  SetInitialConditions(0.0, 0.0, 0.0, _clone->GetApexStress());
}

void Kelvin::SetParameters(double p, double nu, double k, double eta, double dt) {
  Clear();
  parameters.p = p;
  parameters.nu = nu;
  parameters.k = k;
  parameters_kelvin.eta = eta;
  parameters_kelvin.dt = dt;
}

double Kelvin::GetViscosity() {
  return parameters_kelvin.eta;
}

double Kelvin::GetTimeStep() {
  return parameters_kelvin.dt;
}

void Kelvin::SetPreviousShape(Kelvin* _previous) {
  Clear();
  previous = _previous;
  previousCount = 0;
  while (_previous) {
    previousCount++;
    _previous = _previous->previous;
  }
}

double Kelvin::J_dot_s_s(double t) {
  return 1.0 / (1.0 - parameters.nu * parameters.nu) * 1 / parameters_kelvin.eta * exp(- GetEH0() / parameters_kelvin.eta * t);
}

double Kelvin::J_dot_phi_phi(double t) {
  return 1.0 / (1.0 - parameters.nu * parameters.nu) * 1 / parameters_kelvin.eta * exp(- GetEH0() / parameters_kelvin.eta * t);
}

double Kelvin::J_dot_s_phi(double t) {
  return - parameters.nu * J_dot_s_s(t);
}

double Kelvin::J_dot_phi_s(double t) {
  return - parameters.nu * J_dot_phi_phi(t);
}
