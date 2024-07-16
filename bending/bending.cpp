#include "bending.hpp"

void Bending::Clone(Bending* clone) {
  SetIntegrationStep(clone->integration_step);
  SetLaplace(clone->undeformed);
  SetParameters(clone->parameters);
  SetInitialConditions(0, 0, 0, 0, clone->parameters.tau_s_0);
}

void Bending::SetInitialConditions(double s, double r, double z, double m_s_0, double tau_s_0) {
  initial_conditions.s = s;
  initial_conditions.r = r;
  initial_conditions.z = z;
  initial_conditions.m_s = m_s_0;
  initial_conditions.tau_s = tau_s_0;

  parameters.tau_s_0 = tau_s_0;
}

void Bending::SetParameters(BendingParameters params) {
  parameters = params;
  initial_conditions.tau_s = params.tau_s_0;
};

void Bending::SetApexStress(double tau_s_0) {
  initial_conditions.tau_s = tau_s_0;
  parameters.tau_s_0 = tau_s_0;
}

void Bending::Integrate() {
  Clear();
  current = initial_conditions;
  double s_0 = 0;
  L0 = undeformed->L0;
  contact_region = true;

  if (parameters.contact_length < integration_step)
    contact_region = false;

  BendingVariables initial_derivatives;
  ShapeEquations(s_0, &current, &initial_derivatives);
  Push(s_0, current.s, current.r, current.z, 0);

  // Integrate the contact region
  for(;;) {
    if (!contact_region) break;

    if (!RungeKutta(s_0)) {
      MarkInvalid();
      return;
    }

    s_0 += integration_step;

    if (std::isnan(current.r)
        || std::isnan(current.z)
        || std::isnan(current.m_s)
        || std::isnan(current.tau_s)) {
      MarkInvalid();
      return;
    }

    if (constitutive_parameters.lambda_s <= 0) {
      MarkInvalid();
      return;
    }

    if (current.s >= parameters.contact_length) break;
    Push(s_0, current.s, current.r, current.z, current.psi);
  }

  // Handle the splitting point
  if (contact_region) SplittingPoint(s_0, &current);

  contact_region = false;
  // Integrate the free region
  int step_count = (L0 - s_0) / integration_step + 1;
  for(int i = 0; i < step_count; i++) {
    if (i == step_count - 1) SetIntegrationStep(L0 - s_0);

    if (!RungeKutta(s_0)) {
      MarkInvalid();
      return;
    }

    s_0 += integration_step;

    if (std::isnan(current.r)
        || std::isnan(current.z)
        || std::isnan(current.m_s)
        || std::isnan(current.tau_s)) {
      MarkInvalid();
      return;
    }

    if (constitutive_parameters.lambda_s <= 0) {
      MarkInvalid();
      return;
    }

    if (s_0 >= L0) break;
    Push(s_0, current.s, current.r, current.z, current.psi);
  }

  MarkValid();
}
