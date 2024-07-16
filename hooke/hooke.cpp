#include "hooke.hpp"
#include <math.h>

#define min(a, b) ((a) < (b) ? (a) : (b))
#define max(a, b) ((a) < (b) ? (b) : (a))
#define clamp(x, l, u) (min(max((x), (l)), (u)))

Hooke::~Hooke() {
  Clear();
  undeformed = nullptr;
}

void Hooke::Clone(Hooke* _clone) {
  gamma_ratio = _clone->gamma_ratio;
  SetLaplace(_clone->undeformed);
  SetParameters(_clone->GetPressure(),
                _clone->GetPoisson(),
                _clone->GetCompression());

  SetIntegrationStep(_clone->integration_step);
  SetInitialConditions(0.0, 0.0, 0.0, _clone->GetApexStress());

  if (_clone->fitting_dummy) MarkAsFittingDummy();
  if (_clone->shooting_dummy) MarkAsShootingDummy();
}

void Hooke::SetParameters(double p, double nu, double k) {
  Clear();
  parameters.p = p;
  parameters.nu = nu;
  parameters.g = k * (1. - nu) / (1. + nu);
  parameters.k = k;
}

HookeFitParameters Hooke::GetFitParameters() {
  return HookeFitParameters(parameters);
}

HookeFitParameters Hooke::SetParameters(HookeFitParameters* fitParameters, bool clampToFit) {
  Clear();
  HookeFitParameters validParameters;
  validParameters.g = max(fitParameters->g, 1e-3);
  validParameters.tau_s_0 = max(fitParameters->tau_s_0, 1e-3);
  validParameters.k = max(fitParameters->k, 1e-3);

  if (clampToFit) {
    SetParameters(0., (validParameters.k - validParameters.g)/(validParameters.k + validParameters.g), validParameters.k);
    SetApexStress(validParameters.tau_s_0);
  } else {
    SetParameters(0., (fitParameters->k - fitParameters->g)/(fitParameters->k + validParameters.g), fitParameters->k);
    SetApexStress(fitParameters->tau_s_0);
  }
  return validParameters;
};

void Hooke::ResetAdditionalShapeVectors() {
  tau_s.Clear();
  lambda_s.Clear();
  lambda_phi.Clear();
  kappa_phi.Clear();
  tau_phi.Clear();
  kappa_s.Clear();
  ls_times_rs.Clear();
}

void Hooke::Clear() {
  Shape::Clear();
  ResetAdditionalShapeVectors();
  integration_count = 0;
}

void Hooke::Push(double s) {
  // If this is only a dummy used for shooting purposes, we only care about
  // the boundary value and do not need to store these additional data points
  if (shooting_dummy) return;

  Shape::Push(s, current.sd, current.r, current.z, current.psi);
  tau_s.Push(current.tau_s);
  lambda_s.Push(constitutive_parameters.lambda_s);
  lambda_phi.Push(constitutive_parameters.lambda_phi);
  kappa_phi.Push(constitutive_parameters.kappa_phi);
  tau_phi.Push(constitutive_parameters.tau_phi);
  kappa_s.Push(constitutive_parameters.kappa_s);
  ls_times_rs.Push(constitutive_parameters.lambda_s * current.r);
}

void Hooke::SetInitialConditions(double r, double z, double psi, double tau_s) {
  Clear();
  initial_conditions.r = r;
  initial_conditions.z = z;
  initial_conditions.psi = psi;
  initial_conditions.tau_s = tau_s;
  parameters.tau_s_0 = tau_s;
}

std::pair<double, double> Hooke::GetXYPairFromS(double s) {
  if (!r.spline || !z.spline || !acc || s <= 0) return {0, 0};
  return { gsl_spline_eval(r.spline, s, acc),
           gsl_spline_eval(z.spline, s, acc) };
}

double K_init = INITIAL_K;
double nu_init = INITIAL_NU;
double tau_s_0_init = INITIAL_STRESS;

void SetHookeInitialGuess(double _K_init, double _nu_init, double _tau_s_0_init) {
  K_init = _K_init;
  nu_init = _nu_init;
  tau_s_0_init = _tau_s_0_init;
}
