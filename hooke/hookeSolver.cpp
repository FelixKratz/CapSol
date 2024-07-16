#include "hookeSolver.hpp"

void Hooke::Integrate(double s_0, bool reset) {
  if (reset) {
    // Clear the current solution and set the initial conditions.
    // This should only be done when starting the integration
    // from the apex, otherwise we would delete the already existing solution.
    Clear();
    current = initial_conditions;
  }

  // We integrate in the undeformed coordinates
  L0 = undeformed->L0;
  double s = s_0;

  // Store the integration step, as we may change it in the loop below.
  double h = integration_step;

  // Populate all "secondary" variables for the initial condition.
  HookeVariables initial_derivatives;
  ShapeEquations(s, &current, &initial_derivatives);

  if (reset) {
    // Push initial conditions
    Push(s);
  }

  int step_count = (L0 - s_0) / integration_step + 1;
  for (int i = 0; i < step_count; i++) {
    // Adapt the integration step when close to the upper integration bound,
    // such that we integrate exactly up to the capillary in the last
    // integration step.
    if (i == step_count - 1) SetIntegrationStep(L0 - s);

    // Perform the integration step.
    if (!RungeKutta(s)) {
      // A problem has occured in the evaluation of the shape equations.
      MarkInvalid();
      return;
    }

    s += integration_step;

    if (std::isnan(current.r)
        || std::isnan(current.z)
        || std::isnan(current.psi)
        || std::isnan(current.tau_s)
        || std::isnan(constitutive_parameters.lambda_phi)
        || std::isnan(constitutive_parameters.lambda_s)  ) {
      // If anything goes nan we must concede the integration.
      MarkInvalid();
      return;
    }

    if (current.psi < 0
       || current.psi > 1.5*M_PI
       || constitutive_parameters.lambda_s < 0) {
      MarkUnphysical();
    }

    // Save the integration step to the DataSeries objects.
    Push(s);
  }

  // Reset the integration step to the original integration step
  SetIntegrationStep(h);
  MarkValid();

  if (reset) {
    // Only move down the solution if we integrate the full solution from the
    // apex.
    MoveDown();
  }
}

void Hooke::MoveDown() {
  // Shooting dummys are not moved down, nor are they being interpolated,
  // because only the radius at the capillary is important for those shapes
  if (!shooting_dummy) {
    if (z.values.size() > 0) {
      double dz = z.values[z.values.size() - 1];
      for (int i = 0; i < size; i++) z.values[i] -= dz;
      // Only perform sparse interpolation
      Interpolate(true, false);
    }
    else MarkInvalid();
  }
}

void Hooke::Interpolate(bool sparse, bool deformed) {
  // Interpolation is not performed for shooting dummys, as it is not needed.
  if (shooting_dummy) return;

  DataSeries* x = &s;

  // A fitting dummy can never be interpolated in the deformed arc coordinates.
  if (deformed && !fitting_dummy) x = &sd;
  if (deformed && unphysical) {
    printf("Trying to interpolate an unphysical solution...abort\n");
  }

  if (acc) gsl_interp_accel_free(acc);
  acc = gsl_interp_accel_alloc();

  z.Interpolate(x);
  r.Interpolate(x);

  if (!sparse && !fitting_dummy) {
    // Interpolate additional quantities for analysis, this is rarely ever
    // needed, as the data is available at the resolution of the integration in
    // the data series arrays anyways.
    lambda_s.Interpolate(x);
    lambda_phi.Interpolate(x);
    psi.Interpolate(x);
    tau_s.Interpolate(x);
    tau_phi.Interpolate(x);
    lambda_s.Interpolate(x);
    lambda_phi.Interpolate(x);
    kappa_s.Interpolate(x);
    kappa_phi.Interpolate(x);
    ls_times_rs.Interpolate(&s);
  }
}
