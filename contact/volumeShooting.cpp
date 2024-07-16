#include "volumeShooting.hpp"
#include "contact.hpp"

bool ContactVolumeShooting::VolumeShooting(bool parallel) {
  // Perform the shooting for both the upper and the lower shape sequentially
  // and iterate to the fixpoint of both shapes reaching the boundary condition
  
  double target_volume_u = contact->hooke_u.undeformed->Volume();
  double target_volume_d = contact->hooke_d.undeformed->Volume();

  
  contact->SetApexStresses(1.0, 1.0);
  for(;;) {
    contact->parameters.Print();
    // Shoot for the upper shape boundary
    shootingForUpperShape = true;
    if (!Shoot(ShootingRange(target_volume_u, 0.5, 1., parallel, 1e-5))) {
      printf("Upper Shape fail\n");
      return false;
    }

    // Shoot for the lower shape boundary
    shootingForUpperShape = false;
    if (!Shoot(ShootingRange(target_volume_d, 0.5, 1., parallel, 1e-5))) {
      printf("Lower Shape fail\n");
      return false;
    }

    if (abs(contact->hooke_u.Volume() - target_volume_u) < 1e-4
        && abs(contact->hooke_d.Volume() - target_volume_d) < 1e-4)
      break;
  }

  return true;
}

double ContactVolumeShooting::ShootingSingleBoundaryDeviation(double target,
                                                            double parameter) {
  printf("Trying parameter: %f\n", parameter);
  SlipContact contactDummy;
  contactDummy.Clone(contact);
  if (shootingForUpperShape)
    contactDummy.SetApexStresses(parameter,
                                 contactDummy.hooke_d.GetApexStress());
  else 
    contactDummy.SetApexStresses(contactDummy.hooke_u.GetApexStress(),
                                 parameter);

  Hooke* hooke = shootingForUpperShape
                 ? &contactDummy.hooke_u
                 : &contactDummy.hooke_d;

  double deviation = DEVIATION_INVALID;
  if (contactDummy.PressureShooting(true))
    deviation = hooke->Volume() - target;

  printf("Deviation: %.10f\n", deviation);
  return deviation;
}

bool ContactVolumeShooting::ShootingResult(double optimal_parameter) {
  if (shootingForUpperShape)
    contact->SetApexStresses(optimal_parameter,
                             contact->hooke_d.GetApexStress());
  else
    contact->SetApexStresses(contact->hooke_u.GetApexStress(),
                             optimal_parameter);

  contact->PressureShooting(true);

  if (contact->hooke_u.IsInvalid() || contact->hooke_d.IsInvalid())
    return false;

  printf("Converged!\n");
  return true;
}

void ContactVolumeShooting::MarkAsShootingDummy() {
  ShootingMethod::MarkAsShootingDummy();
  contact->hooke_u.MarkAsShootingDummy();
  contact->hooke_d.MarkAsShootingDummy();
}
