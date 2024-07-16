#include "shapeEquations.hpp"

#define P parameters.p
#define NU parameters.nu
#define K parameters.k
#define E_B parameters.e_b
#define F parameters.force
#define GAMMA_CONTACT parameters.gamma_contact
#define GAMMA_WALL parameters.gamma_wall
#define Q_PLUS parameters.q_plus;

#define LAMBDA_S constitutive_parameters.lambda_s
#define LAMBDA_PHI constitutive_parameters.lambda_phi
#define KAPPA_PHI constitutive_parameters.kappa_phi
#define TAU_PHI constitutive_parameters.tau_phi
#define K_S constitutive_parameters.K_s
#define K_PHI constitutive_parameters.K_phi
#define KAPPA_S constitutive_parameters.kappa_s
#define KAPPA_PHI constitutive_parameters.kappa_phi
#define M_PHI constitutive_parameters.m_phi

#define R current->r
#define Z current->z
#define PSI current->psi
#define TAU_S current->tau_s
#define M_S current->m_s
#define S current->s
#define Q current->q

#define R_PRIME derivative->r
#define Z_PRIME derivative->z
#define PSI_PRIME derivative->psi
#define M_S_PRIME derivative->m_s
#define TAU_S_PRIME derivative->tau_s
#define S_PRIME derivative->s
#define Q_PRIME derivative->q

bool Bending::SplittingPoint(double s_0, BendingVariables* current) {
  double kappa_s_0 = undeformed->kappa_s_integration.values[integration_count];

  double Y_2D = K * (2.0 * (1.0 - NU));
  double E = Y_2D / (1.0 - NU*NU);
  
  double delta = (LAMBDA_S / LAMBDA_PHI * E + GAMMA_WALL - (GAMMA_SCALE - GAMMA_CONTACT));
  double alpha = K_S - (LAMBDA_S*LAMBDA_S) / M_S / LAMBDA_PHI * E;
  double xi = NU * K_PHI + kappa_s_0;
  double eta = NU * K_PHI * kappa_s_0;

  double epsilon = E_B * (alpha*alpha + alpha*xi + eta);
  double omega = E_B * delta*delta / (M_S * M_S) + E;
  
  double chi = delta * (E_B / M_S * (2.0 * alpha + xi) - LAMBDA_PHI);

  double p = chi / omega;
  double q = epsilon / omega;
  double lambda_s_plus_1 = - 0.5* p + sqrt(0.25*p*p - q);
  double lambda_s_plus_2 = - 0.5* p - sqrt(0.25*p*p - q);

  double K_s_plus_1 = alpha + lambda_s_plus_1 * delta / M_S;
  double K_s_plus_2 = alpha + lambda_s_plus_2 * delta / M_S;

  double kappa_s_1 = (K_s_plus_1 + kappa_s_0) / lambda_s_plus_1;
  double kappa_s_2 = (K_s_plus_2 + kappa_s_0) / lambda_s_plus_2;

  if (kappa_s_1 < 0 && kappa_s_2 > 0) {
    KAPPA_S = kappa_s_2;
    LAMBDA_S = lambda_s_plus_2;
    K_S = K_s_plus_2;
    M_S = E_B / LAMBDA_PHI * (K_S + NU*K_PHI);
    TAU_S = E / LAMBDA_PHI * (LAMBDA_S - 1 + NU* (LAMBDA_PHI - 1)) + GAMMA_SCALE;
  } else if (kappa_s_1 > 0 && kappa_s_2 < 0) {
    KAPPA_S = kappa_s_1;
    KAPPA_S = kappa_s_1;
    LAMBDA_S = lambda_s_plus_1;
    K_S = K_s_plus_1;
    M_S = E_B / LAMBDA_PHI * (K_S + NU*K_PHI);
    TAU_S = E / LAMBDA_PHI * (LAMBDA_S - 1 + NU* (LAMBDA_PHI - 1)) + GAMMA_SCALE;
  } else {
    printf("We have a problem at the contact point... %f %f\n", kappa_s_1,
                                                                kappa_s_2 );
  }

  Q = Q_PLUS;
  F = M_PI * R * R * P - 2.0 * M_PI * R * Q_PLUS;

  return true;
}

bool Bending::ShapeEquations(double s_0, BendingVariables* current, BendingVariables* derivative) {

  double Y_2D = K * (2.0 * (1.0 - NU));

  double kappa_s_0 = undeformed->kappa_s_integration.values[integration_count];
  double psi_0 = undeformed->psi_integration.values[integration_count];
  double r_0 = undeformed->r_integration.values[integration_count];

  if (s_0 < 1e-8) {

    if (contact_region) {
      // This is only the contact region
      LAMBDA_S = 1.0 / (1.0 - (TAU_S - GAMMA_CONTACT) * (1.0 - NU) / Y_2D);
      K_S = -kappa_s_0;
      K_PHI = K_S;
      M_S = E_B / LAMBDA_S * (K_S + NU * K_PHI);
      M_PHI = E_B / LAMBDA_S * (K_PHI + NU * K_S);
      TAU_S_PRIME = 0;
      M_S_PRIME = 0;
      KAPPA_S = 0;
      KAPPA_PHI = 0;
      PSI = 0;
      Q = 0;
    } else {
      // M_S is a free parameter here.
      // This is just to do sanity checking of the free shape solution
      LAMBDA_S = 1.0 / (1.0 - (TAU_S - GAMMA_SCALE) * (1.0 - NU) / Y_2D);
      K_S = LAMBDA_S * M_S / (E_B * (1.0 + NU));
      K_PHI = K_S;
      M_PHI = E_B / LAMBDA_S * (K_PHI + NU*K_S);
      TAU_S_PRIME = 0;
      M_S_PRIME = 0;
      Q = 0;
    }
  } else {
    double sin_psi = sin(PSI);
    double cos_psi = cos(PSI);
    double kappa_phi_0 = sin(psi_0) / r_0;
    double one_minus_nu_squared = 1.0 - NU * NU;

    LAMBDA_PHI = R / r_0;
    KAPPA_PHI = sin_psi / R;

    // Switch between the shape equations in the contact region and those in the
    // free region (only q changes)
    if (contact_region) {
      LAMBDA_S = (one_minus_nu_squared * LAMBDA_PHI * (TAU_S - GAMMA_CONTACT))
                  / Y_2D + 1.0 - NU * (LAMBDA_PHI - 1.0);

      TAU_PHI = Y_2D / one_minus_nu_squared / LAMBDA_S * (LAMBDA_PHI - 1 
                + NU* (LAMBDA_S - 1)) + GAMMA_CONTACT;

      K_S = -kappa_s_0;
      K_PHI = -kappa_phi_0;
      M_S = E_B / LAMBDA_PHI * (K_S + NU * K_PHI);
      M_PHI = E_B / LAMBDA_S * (K_PHI + NU * K_S);
      KAPPA_S = 0;
      KAPPA_PHI = 0;
      TAU_S_PRIME = LAMBDA_S * (cos_psi / R * (TAU_PHI - TAU_S));
    } else {
      LAMBDA_S = (one_minus_nu_squared * LAMBDA_PHI * (TAU_S - GAMMA_SCALE))
                  / Y_2D + 1.0 - NU * (LAMBDA_PHI - 1.0);

      TAU_PHI = Y_2D / one_minus_nu_squared / LAMBDA_S * (LAMBDA_PHI - 1.0
                + NU* (LAMBDA_S - 1.0)) + GAMMA_SCALE;

      K_PHI = sin_psi / r_0 - kappa_phi_0;
      K_S = LAMBDA_PHI * M_S / E_B - NU * K_PHI;
      KAPPA_S = (K_S + kappa_s_0) / LAMBDA_S;
      M_PHI = E_B / LAMBDA_S * (K_PHI + NU*K_S);

      Q_PRIME = LAMBDA_S * (P - TAU_S * KAPPA_S - TAU_PHI * KAPPA_PHI - Q / R * cos_psi);
      // double q = (- TAU_S * sin_psi
      //            + 0.5*(P * R - F / (M_PI * R))) / cos_psi;

      TAU_S_PRIME = LAMBDA_S * (cos_psi / R * (TAU_PHI - TAU_S) + KAPPA_S * Q);
      M_S_PRIME = LAMBDA_S * (cos_psi / R * (M_PHI - M_S) - Q);
    }
  }

  R_PRIME = LAMBDA_S * cos(PSI);
  Z_PRIME = LAMBDA_S * sin(PSI);
  PSI_PRIME = K_S + kappa_s_0;
  S_PRIME = LAMBDA_S;

  integration_count++;
  return true;
}
