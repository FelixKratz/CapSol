#include "contact.hpp"
#include "config.hpp"

void SlipContact::Clone(SlipContact* clone) {
  if (clone->symmetric) MarkSymmetric();
  SetIntegrationStep(clone->integration_step);
  SetLaplace(clone->hooke_u.undeformed, clone->hooke_d.undeformed);
  SetParameters(clone->parameters);
  SetInitialConditions(0, 0, 0, 0,
                       clone->parameters.hooke_u_parameters.tau_s_0,
                       clone->parameters.hooke_d_parameters.tau_s_0);
}

void SlipContact::SetParameters(SlipContactParameters params) {
  hooke_u.SetParameters(params.hooke_u_parameters);
  hooke_d.SetParameters(params.hooke_d_parameters);
  parameters = params;
  SetApexStresses(params.hooke_u_parameters.tau_s_0,
                  params.hooke_d_parameters.tau_s_0);
}

void SlipContact::SetInitialConditions(double s, double r, double z, double psi, double tau_s_0_u, double tau_s_0_d) {
  initial_conditions.s = s;
  initial_conditions.r = r;
  initial_conditions.z = z;
  initial_conditions.psi = psi;
  initial_conditions.tau_s_u = tau_s_0_u;
  initial_conditions.tau_s_d = tau_s_0_d;

  parameters.hooke_u_parameters.tau_s_0 = tau_s_0_u;
  parameters.hooke_d_parameters.tau_s_0 = tau_s_0_d;

  hooke_u.SetInitialConditions(r, z, psi, tau_s_0_u);
  hooke_d.SetInitialConditions(r, -z, -psi, tau_s_0_d);
}

void SlipContact::SetApexStresses(double tau_s_0_u, double tau_s_0_d) {
  initial_conditions.tau_s_u = tau_s_0_u;
  initial_conditions.tau_s_d = tau_s_0_d;

  parameters.hooke_u_parameters.tau_s_0 = tau_s_0_u;
  parameters.hooke_d_parameters.tau_s_0 = tau_s_0_d;

  hooke_u.SetApexStress(tau_s_0_u);
  hooke_d.SetApexStress(tau_s_0_d);
}

void SlipContact::SetLaplace(Laplace* laplace_u, Laplace* laplace_d) {
  hooke_u.SetLaplace(laplace_u);
  hooke_d.SetLaplace(laplace_d);
}

bool SlipContact::PressureShooting(bool parallel) {
  return pressure_cannon.PressureShooting(parallel);
}

bool SlipContact::VolumeShooting(bool parallel) {
  return volume_cannon.VolumeShooting(parallel);
}

bool SlipContact::AngleShooting(bool parallel) {
  return angle_cannon.AngleShooting(parallel);
}

void SlipContact::SyncAndPush(double s) {
  hooke_u.current.r = current.r;
  hooke_u.current.z = current.z;
  hooke_u.current.psi = current.psi;
  hooke_u.current.tau_s = current.tau_s_u;
  hooke_u.current.sd = current.s;
  hooke_u.Push(s);

  hooke_d.current.r = current.r;
  hooke_d.current.z = - current.z;
  hooke_d.current.psi = - current.psi;
  hooke_d.current.tau_s = current.tau_s_d;
  hooke_d.current.sd = current.s;
  hooke_d.Push(current.s_0_d);
}

void SlipContact::IntegrateLiquid() {
  // Here we only handle the symmetric liquid shape because it is especially
  // easy. The contact region is a straight line of length l, the rest of the
  // shape is simply that of a pendant liquid droplet with appropriate initial
  // conditions given by the contact equations.

  current = initial_conditions;
  SyncAndPush(0.0);

  // Setup all parameters such that we can reuse the elastic splitting point
  // equations for the liquid case
  current.r = parameters.contact_length;
  current.z = 0;
  current.psi = 0;
  current.tau_s_u = parameters.tension_inner_u + parameters.tension_ud;
  current.tau_s_d = parameters.tension_inner_d + parameters.tension_ud;
  hooke_u.parameters.p = hooke_u.undeformed->GetPressure();
  SplittingPoint(parameters.contact_length, &current);

  Laplace dummy(hooke_u.undeformed->GetPressure(),
                hooke_u.undeformed->GetDensity());

  dummy.SetIntegrationStep(integration_step);
  dummy.initial_conditions.psi = hooke_u.current.psi;
  dummy.initial_conditions.r = parameters.contact_length;
  dummy.initial_conditions.z = 0;
  dummy.initial_conditions.s = integration_step;
  if (dummy.Solve()) {
    for (int i = 0; i < dummy.psi.values.size(); i++) {
      hooke_u.current.psi = dummy.psi.values[i];
      hooke_u.current.r = dummy.r.values[i];
      hooke_u.current.z = dummy.z.values[i] - dummy.z.values[0];
      hooke_u.current.tau_s = 1.0;
      hooke_u.current.sd = dummy.s.values[i];
      hooke_u.Push(dummy.s.values[i]);

      hooke_d.current = hooke_u.current;
      hooke_d.Push(dummy.s.values[i]);
    }
    hooke_u.Interpolate(true);
    hooke_d.Interpolate(true);
    hooke_u.MarkValid();
    hooke_d.MarkValid();
  } else {
    hooke_u.MarkInvalid();
    hooke_d.MarkInvalid();
  }
}

void SlipContact::Integrate() {
  if (parameters.hooke_d_parameters.k == 0.0
      && parameters.hooke_u_parameters.k == 0.0) {
    if (!symmetric) { printf("Unimplemented!\n"); exit(1); }
    // This is the purely liquid case. It must be handled explicitly, because
    // there no longer exists a notion of lambda_s, lambda_phi, s_0_u and s_0_d
    IntegrateLiquid();
    return;
  }

  hooke_u.Clear();
  hooke_d.Clear();

  current = initial_conditions;

  double s_0_u = 0;
  SlipContactVariables initial_derivatives;
  ShapeEquations(s_0_u, &current, &initial_derivatives);
  SyncAndPush(s_0_u);

  // Integrate the contact region
  for(;;) {
    if (!RungeKutta(s_0_u)) {
      hooke_u.MarkInvalid();
      hooke_d.MarkInvalid();
      return;
    }

    s_0_u += integration_step;

    if (std::isnan(current.r)
        || std::isnan(current.z)
        || std::isnan(current.psi)
        || std::isnan(current.tau_s_u)
        || std::isnan(current.tau_s_d)) {
      hooke_u.MarkInvalid();
      hooke_d.MarkInvalid();
      return;
    }

    if (hooke_u.constitutive_parameters.lambda_s <= 0
        || hooke_d.constitutive_parameters.lambda_s <= 0) {
      hooke_u.MarkInvalid();
      hooke_d.MarkInvalid();
      return;
    }

    if (current.s >= parameters.contact_length) break;
    SyncAndPush(s_0_u);
  }

  // Handle the splitting point
  SplittingPoint(s_0_u, &current);

  // Force hooke shapes to use the "slow" spline approach to generate
  // the reference radius instead of the fast lookup approach, since we have
  // bricked the lookup table during the contact integration
  // (the upper shape could use the lookup table if the counter is
  // incremented properly)
  hooke_u.use_splines = true;
  hooke_d.use_splines = true;

  // Wrinkling makes the volume shooting algorithm unstable
  hooke_u.enable_wrinkling = true;
  hooke_d.enable_wrinkling = true;

  hooke_d.gamma_ratio = parameters.tension_ratio;

  hooke_u.SetIntegrationStep(integration_step);
  hooke_d.SetIntegrationStep(integration_step);

  // Integrate the remaining solution for the upper and lower shapes,
  // starting from their respective *undeformed* coordinates at the
  // splitting point

  hooke_u.Integrate(s_0_u, false);
  hooke_d.Integrate(current.s_0_d, false);
}
