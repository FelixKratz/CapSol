#include "config.hpp"
#ifndef TEST
#include "main.hpp"
// #include "bending.hpp"
#include "contact.hpp"

// int main(int argc, char **argv) {
//   Bending bending;
//   bending.SetIntegrationStep(INTEGRATION_STEP);
//
//   LaplaceParameters p_laplace(3.999, 0.0);
//   Laplace laplace;
//   laplace.SetParameters(p_laplace);
//   laplace.SetIntegrationStep(INTEGRATION_STEP);
//   laplace.Solve();
//   if (!laplace.IsValid()) {
//     printf("Error: Invalid Laplace Shape...\n");
//   }
//
//   BendingParameters p_bending(0.0, 0, 1.0, 0.01, 0.5, 3.999, 1.0, 0.0);
//   bending.SetLaplace(&laplace);
//   bending.SetParameters(p_bending);
//   bending.SetInitialConditions(0, 0, 0, 0, 1.0);
//   bending.Integrate();
//
//   printf("%d\n", bending.IsValid());
//   return 0;
// }

// int main(int argc, char **argv) {
//   SlipContact contact;
//   contact.SetIntegrationStep(INTEGRATION_STEP);

//   LaplaceParameters p_laplace_u(2.0, 0.0);
//   LaplaceParameters p_laplace_d(2.0, 0.0);

//   HookeParameters p_hooke_u(0.9, 0.00000);
//   HookeParameters p_hooke_d(0.9, 0.00000);


//   Laplace laplace_u;
//   laplace_u.SetParameters(p_laplace_u);
//   laplace_u.SetIntegrationStep(INTEGRATION_STEP);
//   laplace_u.Solve();
//   if (!laplace_u.IsValid()) {
//     printf("Error: Invalid Laplace Shape...\n");
//   }

//   Laplace laplace_d;
//   laplace_d.SetParameters(p_laplace_d);
//   laplace_d.SetIntegrationStep(INTEGRATION_STEP);
//   laplace_d.Solve();
//   if (!laplace_d.IsValid()) {
//     printf("Error: Invalid Laplace Shape...\n");
//   }

//   contact.SetLaplace(&laplace_u, &laplace_d);
//   contact.MarkSymmetric();
//   contact.SetIntegrationStep(INTEGRATION_STEP);
//   contact.SetParameters(SlipContactParameters(p_hooke_u, p_hooke_d, 0.01, 1.0, 0.5, 0.5, 0.0));
//   contact.SetInitialConditions(0, 0, 0, 0, 0.5, 0.5);
//   bool success = contact.PressureShooting(true);
//   contact.GetUpperVolume();
//   printf("%d\n", success);
//   return 0;
// }

int main(int argc, char **argv) {
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
  int count = 33;
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

  if (!sequence.Solve(false)) {
    printf("\n[x] Error generating Kelvin Shape\n");
    exit(1);
  }

  // SlipContact contact;
  // contact.SetIntegrationStep(INTEGRATION_STEP);

  // LaplaceParameters p_laplace_u(3.99999, 0.0);
  // LaplaceParameters p_laplace_d(3.99999, 0.0);

  // HookeParameters p_hooke_u(0.33, 1.0);
  // HookeParameters p_hooke_d(0.33, 1.0);


  // Laplace laplace_u;
  // laplace_u.SetParameters(p_laplace_u);
  // laplace_u.SetIntegrationStep(INTEGRATION_STEP);
  // laplace_u.Solve();
  // if (!laplace_u.IsValid()) {
  //   printf("Error: Invalid Laplace Shape...\n");
  // }

  // Laplace laplace_d;
  // laplace_d.SetParameters(p_laplace_d);
  // laplace_d.SetIntegrationStep(INTEGRATION_STEP);
  // laplace_d.Solve();
  // if (!laplace_d.IsValid()) {
  //   printf("Error: Invalid Laplace Shape...\n");
  // }

  // contact.SetLaplace(&laplace_u, &laplace_d);
  // contact.MarkSymmetric();
  // contact.SetIntegrationStep(INTEGRATION_STEP);
  // contact.SetParameters(SlipContactParameters(p_hooke_u, p_hooke_d, 1.5, 1.0, 0.5, 0.5, 0.0));
  // contact.SetInitialConditions(0, 0, 0, 0, 0.5, 0.5);
  // bool success = contact.AngleShooting(true);
  // printf("%d\n", success);
  // return 0;
}

#endif
