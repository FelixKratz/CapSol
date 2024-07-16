#include "shapeEquations.hpp"

#define P parameters.p
#define NU parameters.nu
// Fit with K_2D and G_2D
// #define NU ((parameters.k - parameters.nu) / (parameters.k + parameters.nu))
#define K parameters.k

#define LAMBDA_S constitutive_parameters.lambda_s
#define LAMBDA_PHI constitutive_parameters.lambda_phi
#define KAPPA_PHI constitutive_parameters.kappa_phi
#define TAU_PHI constitutive_parameters.tau_phi

#define SD current->sd
#define R current->r
#define Z current->z
#define PSI current->psi
#define TAU_S current->tau_s

#define SD_PRIME derivative->sd
#define R_PRIME derivative->r
#define Z_PRIME derivative->z
#define PSI_PRIME derivative->psi
#define TAU_S_PRIME derivative->tau_s

bool Hooke::ShapeEquations(double s_0, HookeVariables* current,
                                     HookeVariables* derivative) {

  double Y_2D = K * (2.0 * (1.0 - NU));
  if (s_0 < 1e-8) {
    LAMBDA_S = 1.0 / (1.0 - (TAU_S - gamma_ratio) * (1.0 - NU) / Y_2D);

    LAMBDA_PHI = LAMBDA_S;
    R_PRIME = LAMBDA_S * cos(PSI);
    Z_PRIME = LAMBDA_S * sin(PSI);
    KAPPA_PHI = P / (2.0 * TAU_S);
    PSI_PRIME = LAMBDA_S * KAPPA_PHI;
    TAU_S_PRIME = 0.0;
    TAU_PHI = TAU_S;
  }
  else {
    // If the radius is too small, we might encounter numerical problems in the
    // below equations.
    if (R < 1e-10) return false;

    double r0;

    if (use_splines)
      r0 = gsl_spline_eval(undeformed->r.spline, s_0, undeformed->acc);
    else
      r0 = undeformed->r_integration.values[integration_count];

    double rho = (flip_gravity ? -1. : 1.)*undeformed->GetDensity()*gamma_ratio;

    double sin_psi = sin(PSI);
    double cos_psi = cos(PSI);
    double one_minus_nu_squared = 1.0 - NU * NU;

    LAMBDA_PHI = R / r0;
    KAPPA_PHI = sin_psi / R;

    LAMBDA_S = (one_minus_nu_squared * LAMBDA_PHI * (TAU_S - gamma_ratio))
                / Y_2D + 1.0 - NU * (LAMBDA_PHI - 1.0);
    TAU_PHI = (Y_2D / one_minus_nu_squared) * (1.0 / LAMBDA_S)
              * ((LAMBDA_PHI - 1.0) + NU * (LAMBDA_S - 1.0)) + gamma_ratio;

    if (enable_wrinkling && TAU_PHI < 0.0) {
      if (!is_in_wrinkling_domain) {
        is_in_wrinkling_domain = true;
        wrinkling_start = s_0;
      }
      LAMBDA_S = (TAU_S * LAMBDA_PHI + Y_2D - gamma_ratio * (1.0 + NU))
                  /(Y_2D - gamma_ratio * 2.0 * NU
                    - one_minus_nu_squared / Y_2D);
      TAU_PHI = 0.0;
    } else if (is_in_wrinkling_domain) {
      is_in_wrinkling_domain = false;
      wrinkling_stop = s_0;
    }

    R_PRIME = LAMBDA_S * cos_psi;
    Z_PRIME = LAMBDA_S * sin_psi;
    PSI_PRIME = (LAMBDA_S / TAU_S) * (P - rho * Z - TAU_PHI * KAPPA_PHI);
    TAU_S_PRIME = - LAMBDA_S * cos_psi * (TAU_S - TAU_PHI) / R;
  }

  SD_PRIME = LAMBDA_S;
  integration_count++;
  return true;
}
