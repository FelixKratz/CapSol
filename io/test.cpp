#ifdef TEST
#include "test.hpp"
#include "phenomenological.h"

void kelvin_fit() {
  printf("\n--- Testing Basic Non Linear Kelvin Simulation ---\n");
  fflush(stdout);

  Laplace laplace;

  laplace.SetParameters(4., 1.);
  laplace.SetIntegrationStep(INTEGRATION_STEP);
  laplace.Solve();

  if (!laplace.IsValid()) {
    printf("[x] Laplace Invalid\n");
    exit(1);
  } 
  printf("[+] Laplace Valid\n");

  printf("[?] Testing Kelvin Sequence                       ");
  fflush(stdout);

  int periods = 1;
  int count = 30;
  double nu = 0.5;
  double K = 2.0;
  double eta = 1.;
  double dtau_s_0 = 0.4;

  KelvinSequence sequence(count);
  sequence.SetLaplace(&laplace);
  double dt = ((double)periods) * 2.* M_PI / ((double)count - 1.);
  sequence.SetParameters(1., nu, K, eta, dt);
  std::vector<double> apexStresses;
  for (int i = 0; i < sequence.count; i++) {
    apexStresses.push_back(1.0 + dtau_s_0*sin(((double)i) * dt));
  }
  sequence.SetApexStresses(&apexStresses);

  if (!sequence.Solve(true)) {
    printf("\n[x] Error generating Kelvin Shape\n");
    exit(1);
  }

  sequence.Interpolate();
  
  std::vector<PointSet> shapes;
  sequence.GetPointSets(&shapes);

  printf("\r[+] Kelvin Sequence Valid                \n");
  fflush(stdout);

  struct kelvin_parameters fit_parameters;
  printf("--- Trying the Kelvin Fit ---\n");

  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
  double err = 0;
  if (fitKelvinSequence(&shapes, &laplace, &fit_parameters, count, periods, &err)) {
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << " took: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;
    fflush(stdout);
  } else {
    printf("\r[x] Kelvin Fit did not converge      \n");
    fflush(stdout);
    exit(4);
  }
}

void hooke_fit() {
  struct laplace_parameters laplace_parameters = { 2.5, 0.3 };
  struct hooke_parameters hooke_parameters = { laplace_parameters.p_L,
                                               laplace_parameters.rho,
                                               1., 0.33, 0.8           };

  printf("\n--- Testing Basic Non Linear Hooke Simulation ---\n");
  fflush(stdout);

  Laplace laplace;
  Hooke hooke;

  laplace.SetParameters(hooke_parameters.p_L, hooke_parameters.rho);
  laplace.SetIntegrationStep(1e-5);
  laplace.Solve();

  if (!laplace.IsValid()) {
    printf("[x] Laplace Invalid\n");
    fflush(stdout);
    exit(1);
  }
  printf("[+] Laplace Valid\n");
  fflush(stdout);

  hooke.SetLaplace(&laplace);
  hooke.SetParameters(1, hooke_parameters.Psi, hooke_parameters.K);
  hooke.SetIntegrationStep(1e-5);
  hooke.SetInitialConditions(0.0, 0.0, 0.0, 1.0);
  bool hookeShooting = hooke.PressureShooting(true);

  if (!hookeShooting) {
    printf("[x] Stress Free Hooke Invalid\n");
    fflush(stdout);
    exit(1);
  }
  printf("[+] Stress Free Hooke Valid\n");
  fflush(stdout);

  hooke.SetLaplace(&laplace);
  hooke.SetParameters(1, hooke_parameters.Psi, hooke_parameters.K);
  hooke.SetIntegrationStep(1e-5);
  hooke.SetInitialConditions(0.0, 0.0, 0.0, hooke_parameters.tau_s_0);
  hookeShooting = hooke.PressureShooting(true);

  if (!hookeShooting) {
    printf("[x] Stressed Hooke Invalid\n");
    fflush(stdout);
    exit(1);
  }
  printf("[+] Stressed Hooke Valid\n");
  fflush(stdout);

  printf("[!] Interpolating Shape");
  fflush(stdout);
  bool deformed = true;
  hooke.Interpolate(true, deformed);

  printf("\r[-] Saving Shapes to disk");
  fflush(stdout);

  PointSet hookePoints;
  generatePointSet(&hooke, &hookePoints, true);

  printf("\r[+] Non Linear Hooke Simulation Working.\n");
  fflush(stdout);
  struct hooke_parameters fit_parameters;
  printf("--- Trying the Hooke Fit ---\n");

  double error = 0;
  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
  laplace.SetParameters(hooke_parameters.p_L, hooke_parameters.rho);
  laplace.SetIntegrationStep(INTEGRATION_STEP);
  laplace.Solve();
  if (fitHookeShape(&laplace, &hookePoints, &fit_parameters, &error)) {
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    printf("\r[+] Hooke Fit Working. SE: %.3e; Benchmark: %.3e;", pow(hooke_parameters.Psi - fit_parameters.Psi, 2) + pow(hooke_parameters.tau_s_0 - fit_parameters.tau_s_0, 2) + pow(hooke_parameters.K - fit_parameters.K, 2), HOOKE_TARGET_PRECISION);
    std::cout << " took: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;
    fflush(stdout);
  } else {
    printf("\r[x] Hooke Fit did not converge      \n");
    fflush(stdout);
    exit(4);
  }
}

void laplace_fit() {
  struct laplace_parameters laplace_parameters = {2.5, 0.1};

  printf("\n--- Testing Basic Laplace Simulation ---\n");

  Laplace laplace(laplace_parameters.p_L, laplace_parameters.rho);
  laplace.SetIntegrationStep(INTEGRATION_STEP);
  laplace.Solve();

  if (!laplace.IsValid()) {
    printf("[x] Laplace Invalid\n");
    fflush(stdout);
    exit(1);
  }

  printf("[+] Laplace Valid\n");
  fflush(stdout);

  printf("[!] Saving Shape to disk");
  fflush(stdout);

  PointSet laplacePoints;
  generatePointSet(&laplace, &laplacePoints);

  printf("\r[+] Laplace Simulation Working.\n");
  fflush(stdout);
  struct laplace_parameters fit_parameters;
  printf("--- Trying the Laplace Fit ---\n");
  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

  double error;
  Laplace laplaceFit;
  if (fitLaplaceShape(&laplacePoints, &laplaceFit, &fit_parameters, &error)) {
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    printf("[+] Laplace Fit Working. SE %.3e; Benchmark %.3e;", pow(laplace_parameters.p_L - fit_parameters.p_L, 2) + pow(laplace_parameters.rho - fit_parameters.rho, 2), LAPLACE_TARGET_PRECISION);
    std::cout << " took: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;
    fflush(stdout);
  }
  else {
    printf("[x] Laplace Fit Failed.");
    fflush(stdout);
    exit(1);
  }
  if (pow(laplace_parameters.p_L - fit_parameters.p_L, 2) + pow(laplace_parameters.rho - fit_parameters.rho, 2) > LAPLACE_TARGET_PRECISION) {
    printf("\r[!] Laplace Fit not precise enough: %f\n", pow(laplace_parameters.p_L - fit_parameters.p_L, 2) + pow(laplace_parameters.rho - fit_parameters.rho, 2));
    fflush(stdout);
    // exit(2);
  }
}

void area_and_volume() {
  struct laplace_parameters laplace_parameters = { 2.5, 0.3 };
  struct hooke_parameters hooke_parameters = { laplace_parameters.p_L,
                                               laplace_parameters.rho,
                                               1., 0.33, 1.0          };

  Laplace laplace;
  Hooke hooke;

  laplace.SetParameters(hooke_parameters.p_L, hooke_parameters.rho);
  laplace.SetIntegrationStep(INTEGRATION_STEP);
  laplace.Solve();

  if (!laplace.IsValid()) {
    printf("[x] Laplace Invalid\n");
    fflush(stdout);
    exit(1);
  }
  printf("[+] Laplace Valid\n");
  fflush(stdout);

  hooke.SetLaplace(&laplace);
  hooke.SetParameters(1, hooke_parameters.Psi, hooke_parameters.K);
  hooke.SetIntegrationStep(INTEGRATION_STEP);
  hooke.SetInitialConditions(0.0, 0.0, 0.0, 1.0);
  bool hookeShooting = hooke.PressureShooting(true);
  hooke.Interpolate(false, true);

  if (!hookeShooting) {
    printf("[x] Stress Free Hooke Invalid\n");
    fflush(stdout);
    exit(1);
  }
  printf("[+] Stress Free Hooke Valid\n");
  fflush(stdout);

  std::cout << "Laplace Area: " << laplace.CalculateArea()
            << std::endl;
}

extern double tau_s_0_init;
extern double K_init;
extern double nu_init;
void real_life_shape() {
  PointSet hookePoints;
  loadPointSet("../test.dat", &hookePoints);

  Laplace laplace;
  laplace.SetParameters(1.0377614536717459, 0.037851000904706);
  laplace.SetIntegrationStep(INTEGRATION_STEP);
  laplace.Solve();

  struct hooke_parameters fit_parameters;
  double error;
  printf("Started fitting\n");
  fitHookeShape(&laplace, &hookePoints, &fit_parameters, &error);
  printf("Finished fitting\n");
}

void test_generate_shape() {
  std::random_device r;
  std::mt19937 generator(r());
  std::uniform_real_distribution<> p_dist(0., 5.);
  std::uniform_real_distribution<> rho_dist(0., 5.);
  std::uniform_real_distribution<> K_dist(0., 50.0);
  std::uniform_real_distribution<> nu_dist(-1.0, 1.0);
  std::uniform_real_distribution<> tau_s_0_dist(0., 0.95);


  int hooke_valid = 0;
  int laplace_valid = 0;
  for(int i = 0; i < 10000; i++) {
    printf("\rGenerating shape %d/%d (%d:%d valid)",i, 10000, laplace_valid, hooke_valid);
    fflush(stdout);
    double K = K_dist(generator);
    double nu = nu_dist(generator);
    double tau_s_0 = tau_s_0_dist(generator);
    double p_L = p_dist(generator);
    double rho = rho_dist(generator);

    while (!isInParameterTriangle(p_L, rho)) {
      p_L = p_dist(generator);
      rho = rho_dist(generator);
    }

    Laplace laplace;
    laplace.SetParameters(p_L, rho);
    laplace.SetIntegrationStep(INTEGRATION_STEP);
    laplace.Solve();
    if (!laplace.IsValid()) {
      continue;
    }
    laplace_valid++;

    PointSet laplace_p;
    generatePointSet(&laplace, &laplace_p, true);

    Hooke hooke;
    hooke.SetLaplace(&laplace);
    hooke.SetParameters(1., nu, K);
    hooke.SetIntegrationStep(INTEGRATION_STEP);
    hooke.SetInitialConditions(0.0, 0.0, 0.0, tau_s_0);
    bool hookeShooting = hooke.PressureShooting(true);
    if (!hookeShooting) {
      continue;
    }
    hooke.Interpolate(true, true);

    PointSet p;
    generatePointSet(&hooke, &p, true);
    if (!hooke.IsInvalid()) { hooke_valid++; }
  }
}

int main(int argc, char **argv) {
  printf("\n\n--- Testing started ---\n");
  // real_life_shape();

  laplace_fit();
  hooke_fit();
  kelvin_fit();

  printf("\n--- Testing finished successful! ---\n\n");
  return 0;
}

#endif
