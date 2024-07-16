#include "pressureShooting.hpp"
#include "contact.hpp"

bool ContactPressureShooting::PressureShooting(bool parallel) {
  // Purely fluid contact needs no shooting...
  if (contact->parameters.hooke_u_parameters.k == 0.0
      && contact->parameters.hooke_d_parameters.k == 0.0) {
    contact->Integrate();
    return !contact->hooke_u.IsInvalid() && !contact->hooke_d.IsInvalid();
  }
  // Perform the shooting for both the upper and the lower shape sequentially
  // and iterate to the fixpoint of both shapes reaching the boundary condition
  
  int max_iterations = 15;

  // TODO: This is an arbitrary choice of the initial pressure for the lower
  // shape. This can cause problems if the pressure of the lower shape leads to
  // no solution being available for the initial integration of the upper
  // shape.
  contact->parameters.hooke_u_parameters.p = 0;
  contact->parameters.hooke_d_parameters.p = 0;
  contact->SetParameters(contact->parameters);
  for(;;) {
    if (max_iterations-- <= 0) return false;
    // Shoot for the upper shape boundary
    shootingForUpperShape = true;
    if (!Shoot(ShootingRange(0.5,
                             0,
                             2.*contact->hooke_u.undeformed->GetPressure(),
                             parallel))) {
      contact->hooke_u.MarkInvalid();
      return false;
    }

    if (!contact->symmetric) {
      // Shoot for the lower shape boundary
      shootingForUpperShape = false;
      // TODO: This assumes that both capillaries have width a, this does
      // not need to be the case. In general, we should search for the capillary
      // radius of the undeformed shape (which might be different for the lower
      // shape).
      if (!Shoot(ShootingRange(0.5,
                               0,
                               2.*contact->hooke_d.undeformed->GetPressure(),
                               parallel))) {
        contact->hooke_d.MarkInvalid();
        return false;
      }
    }

    if (abs(contact->hooke_u.current.r - 0.5) < 1e-8
        && (abs(contact->hooke_d.current.r - 0.5) < 1e-8
           || contact->symmetric))
      break;
  }

  return true;
}

double ContactPressureShooting::ShootingSingleBoundaryDeviation(double target,
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
    return hooke->current.r - target;
  return 0;
}

bool ContactPressureShooting::ShootingResult(double optimal_parameter) {
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

void ContactPressureShooting::MarkAsShootingDummy() {
  ShootingMethod::MarkAsShootingDummy();
  contact->hooke_u.MarkAsShootingDummy();
  contact->hooke_d.MarkAsShootingDummy();
}
