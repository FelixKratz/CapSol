#include "laplace.hpp"

bool Laplace::ShapeEquations(double s, LaplaceVariables* current,
                                       LaplaceVariables* derivative) {
  double p = parameters.p_L;
  double rho = parameters.rho;

  if (current->r < 1e-8 && s < 1e-8) {
    derivative->r = cos(current->psi);
    derivative->z = sin(current->psi);
    derivative->psi = 0.5*(p-rho*current->z);
  }
  else {
    double sin_y2 = sin(current->psi);
    derivative->r = cos(current->psi);
    derivative->z = sin_y2;
    derivative->psi = (p-rho*current->z)-sin_y2/(current->r);
  }

  if (!fitting_dummy) {
    r_integration.Push(current->r);
    psi_integration.Push(current->psi);
    kappa_s_integration.Push(derivative->psi);
    integration_count++;
  }
  return true;
}

void Laplace::Clone(Laplace* _clone) {
  SetParameters(_clone->GetPressure(),
                _clone->GetDensity(),
                _clone->GetOmega());

  SetIntegrationStep(_clone->integration_step);
}

LaplaceFitParameters Laplace::SetParameters(LaplaceFitParameters* fitParameters,
                                            bool clampToFit) {
  SetParameters(fitParameters->p_L, fitParameters->rho);
  return *fitParameters;
}

LaplaceFitParameters Laplace::GetFitParameters() {
  return LaplaceFitParameters(parameters);
}

void Laplace::SetParameters(double p, double rho, int omega) {
  Clear();
  parameters.p_L = p;
  parameters.rho = rho;
  parameters.omega = omega;

  initial_conditions.r = 0;
  initial_conditions.z = 0;
  initial_conditions.psi = 0;
}

double Laplace::CalculateArea() {
  return 2.0 * M_PI * gsl_spline_eval_integ(r.spline, 0., L0, acc);
}

#define LAPLACE_MAX_SHAPE_LENGTH 20.0
bool Laplace::Solve(bool parallel) {
  current = initial_conditions;

  double s = initial_conditions.s;
  int count_bc_crossing_from_left = 0;
  int count_bc_crossing_from_right = 0;

  double z_before_step = 0.0;
  double z_after_step = 0.0;

  double r0_before_step = 0.0;
  double r0_after_step = 0.0;

  bool inNeck = false;
  int neckCount = 0;
  bool inBulge = false;
  int bulgeCount = 0;

  Clear();
  MarkValid();

  // Initial conditions
  LaplaceVariables initial_derivatives;
  ShapeEquations(s, &current, &initial_derivatives);
  Push(s, s, current.r, current.z, current.psi);

  while (s < LAPLACE_MAX_SHAPE_LENGTH) {
    z_before_step = current.z;
    r0_before_step = current.r;

    RungeKutta(s);

    s += integration_step;

    r0_after_step = current.r;
    z_after_step = current.z;

    if (current.r < 0
        // || current.psi > M_PI
        // || current.psi < - M_PI
        ) {
        // || z_after_step < z_before_step) {
      MarkInvalid();
      return false;
    }

    Push(s, s, current.r, current.z, current.psi);

    // Check boundary conditions
    if ((r0_before_step < parameters.capillary_radius)
        && (r0_after_step > parameters.capillary_radius)) {
      count_bc_crossing_from_left++;
      if (parameters.omega == bulgeCount + neckCount + 1) break;
    }
    if ((r0_before_step > parameters.capillary_radius)
        && (r0_after_step < parameters.capillary_radius)) {
      count_bc_crossing_from_right++;
      if (parameters.omega == bulgeCount + neckCount + 1) break;
    }

    // Count necks and bulges
    if (r0_before_step > r0_after_step) {
      inNeck = true;
      if (inBulge) {
        inBulge = false;
        bulgeCount++;
      }
    }
    if (r0_before_step < r0_after_step) {
      inBulge = true;
      if (inNeck) {
        inNeck = false;
        neckCount++;
      }
    }

    if (bulgeCount + neckCount > parameters.omega
        || s >= LAPLACE_MAX_SHAPE_LENGTH - integration_step) {
      MarkInvalid();
      return false;
    };
  }

  // Check if the boundary condition is satisfied and if target omega is found
  if (count_bc_crossing_from_right + count_bc_crossing_from_left > 0
      && (parameters.omega == bulgeCount + neckCount + 1)) {
    // Only perform sparse inverpolation here, since we shift the solution once
    // more and perform a full interpolation afterwards
    Interpolate(true);

    double h0 = integration_step / 2.0;
    // All even bulge + neck counts go through the capillary from left to right
    double s_search = s - h0;
    bool sign = current.r > parameters.capillary_radius;
    bool left_sign = gsl_spline_eval(r.spline, s - integration_step, acc)
                     > parameters.capillary_radius;
    bool right_sign = sign;
    bool center_sign;

    // Bisection algorithm to find the actual bc zero crossing.
    // Checking the truth value of h0 is a machine precision agnostic check. It
    // will return true as long as the machine can resove it and false when h0
    // is smaller than machine precision.
    do {
      center_sign = gsl_spline_eval(r.spline, s_search, acc)
                    > parameters.capillary_radius;

      if (left_sign == center_sign && center_sign != right_sign) {
        left_sign = center_sign;
        s_search += 0.5*h0;
      } else if (left_sign != center_sign && center_sign == right_sign) {
        right_sign = center_sign;
        s_search -= 0.5*h0;
      } else {
        // Should never happen, except if the spline can no longer resolve h0.
        break;
      }
      h0 *= 0.5;
    } while (h0);

    L0 = s_search;

    if (L0 > s - integration_step) {
      // Repeat the last integration step with the appropriate step width
      // to get the actual values of all quantities at the capillary
      Pop();
      s -= integration_step;
      double h = integration_step;
      SetIntegrationStep(L0 - s);
      current.r = r.values.back();
      current.z = z.values.back();
      current.psi = psi.values.back();
      RungeKutta(s);
      s += integration_step;
      Push(s, s, current.r, current.z, current.psi);
      SetIntegrationStep(h);
    }

    MarkValid();
    MoveDown();

    if (IsValid()) Interpolate();
    return true;
  }

  // The shape is too long and thus considered invalid.
  MarkInvalid();
  return false;
}

void Laplace::MoveDown() {
  double dz = z.values.back();
  if (std::isnan(dz)) {
    std::cout << "Error in spline eval" << std::endl;
    std::cout << "p: " << GetPressure() << " "
              << " rho: " << GetDensity() << " "
              << " omega: " << GetOmega() << std::endl;
    MarkInvalid();
    return;
  }

  if (dz < -1e-10) {
    MarkInvalid();
    return;
  }

  for (int i = 0; i < size; i++) z.values[i] -= dz;
}

void Laplace::Interpolate(bool sparse) {
  if (acc) gsl_interp_accel_free(acc);
  acc = gsl_interp_accel_alloc();

  r.Interpolate(&s);
  z.Interpolate(&s);

  // These are only evaluated if the shape is not a dummy and the
  // interpolation is not explicitly requested as being sparse
  if (!sparse && !fitting_dummy) {
    // Allocate the splines
    psi.Interpolate(&s);
    z.Interpolate(&s);
  }
}

double p_L_init = INITIAL_P_L;
double rho_init = INITIAL_RHO;
void SetLaplaceInitialGuess(double p_L, double rho) {
  p_L_init = p_L;
  rho_init = rho;
}
