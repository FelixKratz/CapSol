#include "interface.hpp"
#include <string>

extern double p_L_init;
extern double rho_init;
bool fitLaplaceShapes(std::vector<PointSet*> laplacePointSets, Laplace* laplace, struct laplace_parameters* fit_parameters, double* err) {
  laplace->SetParameters(p_L_init, rho_init);
  laplace->SetIntegrationStep(INTEGRATION_STEP);
  laplace->Solve();

  if (laplace->IsInvalid() && (p_L_init != INITIAL_P_L
                               || rho_init != INITIAL_RHO)) {
    laplace->SetParameters(INITIAL_P_L, INITIAL_RHO);
    laplace->SetIntegrationStep(INTEGRATION_STEP);
    laplace->Solve();
  }

  if (!laplace->Fit(laplacePointSets, err)) {
    *fit_parameters = {0, 0};
    return false;
  }
  else {
    *fit_parameters = { laplace->GetPressure(), laplace->GetDensity() };
    return true;
  }

  *fit_parameters = {0, 0};
  return false;
}
bool fitLaplaceShapes(std::vector<PointSet*> laplacePointSets, struct laplace_parameters* fit_parameters, double* err) {
  Laplace laplace;
  return fitLaplaceShapes(laplacePointSets, &laplace, fit_parameters, err);
}

bool fitLaplaceShape(PointSet* laplacePointSet, Laplace* laplace, struct laplace_parameters* fit_parameters, double* err) {
  std::vector<PointSet*> pointSets;
  pointSets.push_back(laplacePointSet);
  return fitLaplaceShapes(pointSets, laplace, fit_parameters, err);
}

bool fitLaplaceShape(PointSet* laplacePointSet, struct laplace_parameters* fit_parameters, double* err) {
  Laplace laplace;
  return fitLaplaceShape(laplacePointSet, &laplace, fit_parameters, err);
}

extern double K_init;
extern double nu_init;
extern double tau_s_0_init;
bool fitHookeShapes(Laplace* laplace, std::vector<PointSet*> hookePointSet, struct hooke_parameters* fit_parameters, double* err) {
  bool success = false;

  Hooke hooke;
  hooke.SetLaplace(laplace);
  hooke.SetParameters(laplace->GetPressure(), nu_init, K_init);
  hooke.SetIntegrationStep(INTEGRATION_STEP);
  hooke.SetInitialConditions(0.0, 0.0, 0.0, tau_s_0_init);
  bool hookeShooting = hooke.PressureShooting(true);
  if (!hookeShooting && (K_init != INITIAL_K
                         || nu_init != INITIAL_NU
                         || tau_s_0_init != INITIAL_STRESS)) {
    hooke.SetParameters(laplace->GetPressure(), INITIAL_NU, INITIAL_K);
    hooke.SetIntegrationStep(INTEGRATION_STEP);
    hooke.SetInitialConditions(0.0, 0.0, 0.0, INITIAL_STRESS);
    hookeShooting = hooke.PressureShooting(true);
  }

  if (!hookeShooting) {
    printf("[!] Invalid Initial Hooke Shape!\n");
    *fit_parameters = {0, 0, 0, 0, 0};
    return false;
  }

  success = hooke.Fit(hookePointSet, err);

  if (success) {
    *fit_parameters = { laplace->GetPressure(), laplace->GetDensity(), hooke.GetCompression(), hooke.GetPoisson(), hooke.GetApexStress() };
    return true;
  }
  *fit_parameters = {0, 0, 0, 0, 0};
  return false;
}

bool fitHookeShape(Laplace* laplace, PointSet* hookePointSet, struct hooke_parameters* fit_parameters, double* err) {
  std::vector<PointSet*> pointSets;
  pointSets.push_back(hookePointSet);
  return fitHookeShapes(laplace, pointSets, fit_parameters, err);
}

extern double eta_init;
extern double stress_amplitude_init;
bool fitKelvinSequence(std::vector<PointSet>* kelvinPointSet, Laplace* laplace, struct kelvin_parameters* fit_parameters, int count, int periods, double* err) {
  KelvinSequence kelvinSequence(count);

  std::vector<double> apexStresses;

  double dt = ((double)periods) * 2.*M_PI / ((double)(count) - 1.0);
  for (int i = 0; i < count; i++) {
    apexStresses.push_back(1. + stress_amplitude_init*sin(((double)i) * dt));
  }
  kelvinSequence.SetLaplace(laplace);
  kelvinSequence.SetParameters(1., nu_init, K_init, eta_init, dt);
  kelvinSequence.SetApexStresses(&apexStresses);
  if (!kelvinSequence.Solve(true)) {
    printf("\n[!] Intial Kelvin Invalid!\n");
    return false;
  }

  bool success = KelvinFit(&kelvinSequence, kelvinPointSet, count, err);

  if (success) {
    *fit_parameters = { laplace->GetPressure(),
                        laplace->GetDensity(),
                        kelvinSequence.shapes[0].GetCompression(),
                        kelvinSequence.shapes[0].GetPoisson(),
                        kelvinSequence.stressAmplitude,
                        kelvinSequence.shapes[0].GetViscosity() };
    return true;
  }

  return false;
}
