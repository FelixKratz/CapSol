#include "shapeEquations.hpp"

#define P parameters.p
#define NU parameters.nu
// Fit with K_2D and G_2D
// #define NU ((parameters.k - parameters.nu) / (parameters.k + parameters.nu))
#define K parameters.k
#define ETA parameters_kelvin.eta
#define DT parameters_kelvin.dt

#define LAMBDA_S constitutive_parameters.lambda_s
#define LAMBDA_PHI constitutive_parameters.lambda_phi
#define KAPPA_PHI constitutive_parameters.kappa_phi
#define TAU_PHI constitutive_parameters.tau_phi
#define TAU_PHI_V (constitutive_parameters.tau_phi - GAMMA_SCALE)

#define R current->r
#define Z current->z
#define PSI current->psi
#define TAU_S current->tau_s
#define TAU_S_V (current->tau_s - GAMMA_SCALE)

#define SD_PRIME derivative->sd
#define R_PRIME derivative->r
#define Z_PRIME derivative->z
#define PSI_PRIME derivative->psi
#define TAU_S_PRIME derivative->tau_s

bool Kelvin::ShapeEquations(double s, HookeVariables* current, HookeVariables* derivative) {
  if (use_generalized_equation)
    return ShapeEquationsGeneralized(s, current, derivative);

  double Y_2D = K * (2.0 * (1.0 - NU));
  double rho = undeformed->GetDensity();
  double r0 = undeformed->r_integration.values[integration_count];

  double lambda_s_old = 1.0, lambda_phi_old = 1.0;
  if (previous) {
    lambda_s_old = previous->lambda_s_detailed.values[integration_count];
    lambda_phi_old = previous->lambda_phi_detailed.values[integration_count];

    // The spline implementation introduces numerical instabilities, thus the
    // above method is used, where all intermediate steps are saved to an
    // additional array.
    // double lambda_s_old_ = gsl_spline_eval(previous->lambda_s.spline,
    //                                        s,
    //                                        acc);
    // double lambda_phi_old_ = gsl_spline_eval(previous->lambda_phi.spline,
    //                                          s,
    //                                          acc);
  }

  if (s < 1.0e-8) {
    LAMBDA_S = (ETA / Y_2D * lambda_s_old/DT
        + 1.0)/(1.0 + (GAMMA_SCALE - TAU_S) * (1.0 - NU) / Y_2D
          + ETA / Y_2D / DT);

    LAMBDA_PHI = LAMBDA_S;
    R_PRIME = LAMBDA_S * cos(PSI);
    Z_PRIME = LAMBDA_S * sin(PSI);
    KAPPA_PHI = (P - rho * Z) / (2. * TAU_S);
    PSI_PRIME = LAMBDA_S * KAPPA_PHI;
    TAU_PHI = TAU_S;
    TAU_S_PRIME = 0.0;
  }
  else {
    // If the radius is too small, we might encounter numerical problems in the
    // below equations.
    if (R < 1e-10) return false;

    double sin_psi = sin(PSI);
    double cos_psi = cos(PSI);
    LAMBDA_PHI = R / r0;
    KAPPA_PHI = sin_psi / R;
    LAMBDA_S = ((TAU_S - GAMMA_SCALE)*(1.0 - NU*NU)*LAMBDA_PHI / Y_2D + 1.0 + NU
        + ETA/Y_2D /DT*(lambda_s_old + NU*lambda_phi_old))/(1.0 + ETA/Y_2D/DT)
        - NU*LAMBDA_PHI;
    
    TAU_PHI = 1.0/LAMBDA_S * (Y_2D / (1.0 - NU * NU) * (LAMBDA_PHI - 1.0
          + NU*(LAMBDA_S - 1.0)) + ETA/DT/(1.0 - NU * NU)*(LAMBDA_PHI
          - lambda_phi_old + NU*(LAMBDA_S - lambda_s_old))) + GAMMA_SCALE;

    R_PRIME = LAMBDA_S * cos_psi;
    Z_PRIME = LAMBDA_S * sin_psi;
    PSI_PRIME = (LAMBDA_S / TAU_S) * (P - rho * Z - TAU_PHI * KAPPA_PHI);
    TAU_S_PRIME = - LAMBDA_S * cos_psi * (TAU_S - TAU_PHI) / R;
  }

  SD_PRIME = LAMBDA_S;
  lambda_s_detailed.Push(LAMBDA_S);
  lambda_phi_detailed.Push(LAMBDA_PHI);
  integration_count++;
  return true;
}

bool Kelvin::ShapeEquationsGeneralized(double s, HookeVariables* current, HookeVariables* derivative) {
  double rho = undeformed->GetDensity();
  double r0 = undeformed->r_integration.values[integration_count];

  double R_s_s = 0., R_phi_phi = 0., R_s_phi = 0., R_phi_s = 0.;
  Kelvin* curr_prev = previous;
  int n = 1;
  while(curr_prev && curr_prev->previous) {
    R_s_s += J_dot_s_s(n*DT) * (curr_prev->tau_s_detailed.values[integration_count] - GAMMA_SCALE) * curr_prev->lambda_phi_detailed.values[integration_count];
    R_phi_phi += J_dot_phi_phi(n*DT) * (curr_prev->tau_phi_detailed.values[integration_count] - GAMMA_SCALE) * curr_prev->lambda_s_detailed.values[integration_count];
    R_phi_s += J_dot_phi_s(n*DT) * (curr_prev->tau_s_detailed.values[integration_count] - GAMMA_SCALE) * curr_prev->lambda_phi_detailed.values[integration_count];
    R_s_phi += J_dot_s_phi(n*DT) * (curr_prev->tau_phi_detailed.values[integration_count] - GAMMA_SCALE) * curr_prev->lambda_s_detailed.values[integration_count];
    curr_prev = curr_prev->previous;
    n++;
  }

  if (s < 1.0e-8) {
    if (n > 1) {
      LAMBDA_S = (1.0 - J_dot_s_phi(0) / J_dot_phi_phi(0) + DT * (R_s_s - J_dot_s_phi(0) / J_dot_phi_phi(0) * R_phi_phi + R_s_phi - J_dot_s_phi(0) / J_dot_phi_phi(0) * R_phi_s)) / (1.0 - J_dot_s_phi(0) / J_dot_phi_phi(0)
      - DT / 2.0 * TAU_S_V * (J_dot_s_s(0) - J_dot_s_phi(0) / J_dot_phi_phi(0) * J_dot_phi_s(0)));

      TAU_PHI = TAU_S;
    } else {
      LAMBDA_S = 1.0;
      TAU_PHI = GAMMA_SCALE;
    }
    LAMBDA_PHI = LAMBDA_S;
    R_PRIME = LAMBDA_S * cos(PSI);
    Z_PRIME = LAMBDA_S * sin(PSI);
    KAPPA_PHI = (P - rho * Z) / (2. * TAU_S);
    PSI_PRIME = LAMBDA_S * KAPPA_PHI;
    TAU_S_PRIME = 0.0;
  }
  else {
    // If the radius is too small, we might encounter numerical problems in the
    // below equations.
    if (R < 1e-10) return false;

    double sin_psi = sin(PSI);
    double cos_psi = cos(PSI);
    LAMBDA_PHI = R / r0;
    KAPPA_PHI = sin_psi / R;
    if (n > 1) {
      LAMBDA_S = 1 + J_dot_s_phi(0) / J_dot_phi_phi(0) * (LAMBDA_PHI - 1)
       + DT/2.0 * (
       LAMBDA_PHI * TAU_S_V * (J_dot_s_s(0) - J_dot_s_phi(0) / J_dot_phi_phi(0) * J_dot_phi_s(0))
       + 2.0 * (R_s_s + R_s_phi - (J_dot_s_phi(0) / J_dot_phi_phi(0)) * (R_phi_phi + R_phi_s))); 
  
      TAU_PHI = 1.0 / (LAMBDA_S * DT * J_dot_phi_phi(0)) * ( 2.0 * (LAMBDA_PHI - 1) - DT * (
        J_dot_phi_s(0) * LAMBDA_PHI * TAU_S_V + 2.0*R_phi_s + 2.0*R_phi_phi)) + GAMMA_SCALE;
    } else {
      LAMBDA_S = 1.0;
      TAU_PHI = GAMMA_SCALE;
    }

    R_PRIME = LAMBDA_S * cos_psi;
    Z_PRIME = LAMBDA_S * sin_psi;
    PSI_PRIME = (LAMBDA_S / TAU_S) * (P - rho * Z - TAU_PHI * KAPPA_PHI);
    TAU_S_PRIME = - LAMBDA_S * cos_psi * (TAU_S - TAU_PHI) / R;
  }

  SD_PRIME = LAMBDA_S;
  lambda_s_detailed.Push(LAMBDA_S);
  lambda_phi_detailed.Push(LAMBDA_PHI);
  tau_s_detailed.Push(TAU_S);
  tau_phi_detailed.Push(TAU_PHI);
  integration_count++;
  return true;
}
