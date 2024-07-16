#include "angleShooting.hpp"
#include "contact.hpp"

bool ContactAngleShooting::AngleShooting(bool parallel) {
  // Purely fluid contact needs no shooting...
  if (contact->parameters.hooke_u_parameters.k == 0.0
      && contact->parameters.hooke_d_parameters.k == 0.0) {
    printf("Purely fluid angle shooting is not implemented\n");
    exit(1);
  }

  contact->parameters.hooke_u_parameters.p = 0;
  contact->parameters.hooke_d_parameters.p = 0;
  contact->SetParameters(contact->parameters);
  shootingForUpperShape = true;
  if (!Shoot(ShootingRange(M_PI / 2.0,
                           0,
                           2.*contact->hooke_u.undeformed->GetPressure(),
                           parallel))) {
    contact->hooke_u.MarkInvalid();
    return false;
  }

  if (!contact->symmetric) {
    printf("Contact Needs to be symmetric for angle shooting routine!\n");
    exit(1);
  }

  return true;
}

double ContactAngleShooting::ShootingSingleBoundaryDeviation(double target,
                                                    double parameter) {
  SlipContact contactDummy;

  contactDummy.Clone(contact);

  if (shootingForUpperShape) {
    contactDummy.parameters.hooke_u_parameters.p = parameter;
    if (contact->symmetric) {
      contactDummy.parameters.hooke_d_parameters.p = parameter;
    }
  }
  else
    contactDummy.parameters.hooke_d_parameters.p = parameter;

  contactDummy.SetParameters(contactDummy.parameters);

  Hooke* hooke = shootingForUpperShape
                 ? &contactDummy.hooke_u
                 : &contactDummy.hooke_d;

  contactDummy.MarkAsPressureShootingDummy();
  contactDummy.Integrate();

  if (hooke->IsInvalid())
    return DEVIATION_INVALID;
  else
    return hooke->current.psi - target;
  return 0;
}

bool ContactAngleShooting::ShootingResult(double optimal_parameter) {
  if (shootingForUpperShape) {
    contact->parameters.hooke_u_parameters.p = optimal_parameter;
    if (contact->symmetric) {
      contact->parameters.hooke_d_parameters.p = optimal_parameter;
    }
  }
  else
    contact->parameters.hooke_d_parameters.p = optimal_parameter;

  contact->SetParameters(contact->parameters);
  contact->Integrate();

  if (contact->hooke_u.IsInvalid() || contact->hooke_d.IsInvalid()) return false;
  return true;
}

void ContactAngleShooting::MarkAsShootingDummy() {
  ShootingMethod::MarkAsShootingDummy();
  contact->hooke_u.MarkAsShootingDummy();
  contact->hooke_d.MarkAsShootingDummy();
}
