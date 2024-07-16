#include "pressureShooting.hpp"

// The actual pressure shooting implementation
bool Hooke::PressureShooting(bool parallel, double min, double max) {
  // The Shoot method is implemented in the generalized shooting method.
  if (!Shoot(ShootingRange(0.5, min, max, parallel))) {
    MarkInvalid();
    return false;
  }
  return true;
}

// Overload for the pressure shooting, which specifies a typical search range
// for the pressure.
bool Hooke::PressureShooting(bool parallel) {
  return PressureShooting(parallel, 0., 1.5*undeformed->GetPressure());
}


// Calculate the boundary deviation from the target and the parameter.
// This function is called from the generalized shooting method the Hooke
// object inherits from.
double Hooke::ShootingSingleBoundaryDeviation(double target, double parameter) {
  Hooke hookeDummy;
  hookeDummy.Clone(this);
  hookeDummy.SetPressure(parameter);
  hookeDummy.MarkAsShootingDummy();
  hookeDummy.Integrate();

  if (hookeDummy.IsInvalid())
    return DEVIATION_INVALID;
  else
    return hookeDummy.current.r - target;
}

// Calculate the resulting shape for the best parameter found in the generalized
// shooting method
bool Hooke::ShootingResult(double optimal_parameter) {
  SetPressure(optimal_parameter);
  Integrate();

  // We allow unphysical solutions iff the result is used as a dummy in
  // fitting or shooting.
  if ((IsUnphysical() && !fitting_dummy && !shooting_dummy)
     || abs(r.values.back() - 0.5) > 1e-2) {
    MarkInvalid();
    return false;
  }

  if (IsInvalid()) return false;

  return true;
}
