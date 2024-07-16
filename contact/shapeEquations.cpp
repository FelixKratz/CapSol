#include "shapeEquations.hpp"

#define P_U hooke_u.parameters.p
#define NU_U hooke_u.parameters.nu
#define K_U hooke_u.parameters.k

#define P_D hooke_d.parameters.p
#define NU_D hooke_d.parameters.nu
#define K_D hooke_d.parameters.k

#define LAMBDA_S_U hooke_u.constitutive_parameters.lambda_s
#define LAMBDA_PHI_U hooke_u.constitutive_parameters.lambda_phi
#define KAPPA_PHI_U hooke_u.constitutive_parameters.kappa_phi
#define TAU_PHI_U hooke_u.constitutive_parameters.tau_phi

#define LAMBDA_S_D hooke_d.constitutive_parameters.lambda_s
#define LAMBDA_PHI_D hooke_d.constitutive_parameters.lambda_phi
#define KAPPA_PHI_D hooke_d.constitutive_parameters.kappa_phi
#define TAU_PHI_D hooke_d.constitutive_parameters.tau_phi

#define R current->r
#define Z current->z
#define PSI current->psi
#define TAU_S_U current->tau_s_u
#define TAU_S_D current->tau_s_d
#define S_0_D current->s_0_d
#define S current->s

#define R_PRIME derivative->r
#define Z_PRIME derivative->z
#define PSI_PRIME derivative->psi
#define TAU_S_PRIME_U derivative->tau_s_u
#define TAU_S_PRIME_D derivative->tau_s_d
#define S_0_PRIME_D derivative->s_0_d
#define S_PRIME derivative->s

#define GAMMA_RATIO parameters.tension_ratio
#define GAMMA_I_U parameters.tension_inner_u
#define GAMMA_I_D parameters.tension_inner_d
#define GAMMA_UD parameters.tension_ud

bool SlipContact::SplittingPoint(double s_0, SlipContactVariables* current) {
  double tau_s_u_l_plus = TAU_S_U - GAMMA_UD
                          + GAMMA_SCALE - GAMMA_I_U;

  double tau_s_d_l_plus = TAU_S_D - GAMMA_UD
                          + GAMMA_RATIO*(GAMMA_SCALE - GAMMA_I_D);

  double delta_psi_u = acos(1. + (pow(TAU_S_U + TAU_S_D
                                 - tau_s_u_l_plus, 2.0)
                                - tau_s_d_l_plus * tau_s_u_l_plus)
                                / (2.*tau_s_u_l_plus * (TAU_S_U + TAU_S_D)));

  double delta_psi_d = asin(sin(delta_psi_u) * tau_s_u_l_plus / tau_s_d_l_plus);

  double psi_u_l_plus = PSI + abs(delta_psi_u); 
  double psi_d_l_plus = - PSI + abs(delta_psi_d); 

  parameters.force = M_PI * P_U * R * R
                     - 2. * M_PI * R * tau_s_u_l_plus * sin(psi_u_l_plus);
  hooke_u.current.tau_s = tau_s_u_l_plus;

  hooke_u.current.psi = psi_u_l_plus;
  hooke_u.current.r = R;
  hooke_u.current.z = Z;

  hooke_d.current.tau_s = tau_s_d_l_plus;

  hooke_d.current.psi = psi_d_l_plus;
  hooke_d.current.r = R;
  hooke_d.current.z = -Z;

  return true;
}

bool SlipContact::ShapeEquations(double s_0, SlipContactVariables* current,
                                             SlipContactVariables* derivative) {

  // Before accessing the spline we need to check if we already beyond the
  // iterpolation range and abort accordingly.
  // We need to do this explicitly here, because we can not predict S_0_D
  // during an RK4 step as we are integrating in s_0_u.
  if (s_0 > hooke_u.undeformed->L0 || S_0_D > hooke_d.undeformed->L0)
    return false;

  double Y_2D_U = K_U * (2.0 * (1.0 - NU_U));
  double Y_2D_D = K_D * (2.0 * (1.0 - NU_D));
  double gamma_u = GAMMA_I_U + GAMMA_UD;
  double gamma_d = GAMMA_I_D + GAMMA_UD;
  double p_u = P_U;
  double p_d = symmetric ? P_U : P_D;

  if (s_0 < 1e-8) {
    LAMBDA_S_U = 1.0 / (1.0 - (TAU_S_U - gamma_u) * (1.0 - NU_U) / Y_2D_U);
    LAMBDA_S_D = 1.0 / (1.0 - (TAU_S_D - gamma_d) * (1.0 - NU_D) / Y_2D_D);
    LAMBDA_PHI_U = LAMBDA_S_U;
    LAMBDA_PHI_D = LAMBDA_S_D;

    R_PRIME = LAMBDA_S_U * cos(PSI);
    Z_PRIME = LAMBDA_S_U * sin(PSI);

    KAPPA_PHI_U = (p_u - p_d) / (2.0 * (TAU_S_U + TAU_S_D));

    PSI_PRIME = LAMBDA_S_U * KAPPA_PHI_U;
    TAU_S_PRIME_U = 0.0;
    TAU_S_PRIME_D = 0.0;
    TAU_PHI_U = TAU_S_U;
    TAU_PHI_D = TAU_S_D;
  } else {
    double sin_psi = sin(PSI);
    double cos_psi = cos(PSI);

    // Calculate all available properties of the upper shape
    double rho_u = hooke_u.undeformed->GetDensity();
    double one_minus_nu_squared_u = 1.0 - NU_U * NU_U;
    double r0_u = gsl_spline_eval(hooke_u.undeformed->r.spline,
                                  s_0,
                                  hooke_u.undeformed->acc);

    LAMBDA_PHI_U = R / r0_u;
    KAPPA_PHI_U = sin_psi / R;
    LAMBDA_S_U = (one_minus_nu_squared_u * LAMBDA_PHI_U
                  * (TAU_S_U - gamma_u))
                 / Y_2D_U + 1.0 - NU_U * (LAMBDA_PHI_U - 1.0);

    TAU_PHI_U = (Y_2D_U / one_minus_nu_squared_u) * (1.0 / LAMBDA_S_U)
                * ((LAMBDA_PHI_U - 1.0) + NU_U * (LAMBDA_S_U - 1.0))
                + gamma_u;

    R_PRIME = LAMBDA_S_U * cos_psi;
    Z_PRIME = LAMBDA_S_U * sin_psi;
    TAU_S_PRIME_U = - LAMBDA_S_U * cos_psi * (TAU_S_U - TAU_PHI_U) / R;

    double r0_d = gsl_spline_eval(hooke_d.undeformed->r.spline,
                                  S_0_D,
                                  hooke_d.undeformed->acc);

    // Calculate all available properties of the lower shape
    // We need to scale the density by the tension ratio, because the
    // density in the undeformed state is defined with \gamma_d as a scale.
    // Additionally, we flip the axis of gravity, since the sign is already
    // contained in the shape equations.
    double rho_d = - hooke_d.undeformed->GetDensity() * GAMMA_RATIO;
    double one_minus_nu_squared_d = 1.0 - NU_D * NU_D;

    LAMBDA_PHI_D = R / r0_d;
    LAMBDA_S_D = (one_minus_nu_squared_d * LAMBDA_PHI_D
                 * (TAU_S_D - gamma_d))
                / Y_2D_D + 1.0 - NU_D * (LAMBDA_PHI_D - 1.0);

    TAU_PHI_D = (Y_2D_D / one_minus_nu_squared_d) * (1.0 / LAMBDA_S_D)
                * ((LAMBDA_PHI_D - 1.0) + NU_D * (LAMBDA_S_D - 1.0))
                + gamma_d;

    TAU_S_PRIME_D = - LAMBDA_S_U * cos_psi * (TAU_S_D - TAU_PHI_D) / R;

    // Finally calculate the arc angle derivative
    PSI_PRIME = LAMBDA_S_U / (TAU_S_U + TAU_S_D)
                * ((p_u - rho_u * Z - p_d + rho_d * Z)
                - KAPPA_PHI_U * (TAU_PHI_U + TAU_PHI_D));
  }

  // Lower undeformed coordinate evolution equation
  S_0_PRIME_D = LAMBDA_S_U / LAMBDA_S_D;

  // Deformed arc length coordinate evolution equation
  S_PRIME = LAMBDA_S_U;

  // Trivial coupling condition
  KAPPA_PHI_D = - KAPPA_PHI_U;
  return true;
}
