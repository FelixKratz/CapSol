#include "pressureShooting.hpp"

double Kelvin::ShootingSingleBoundaryDeviation(double target, double parameter) {
  Kelvin kelvinDummy;
  kelvinDummy.Clone(this);
  kelvinDummy.SetPressure(parameter);
  kelvinDummy.MarkAsShootingDummy();
  kelvinDummy.Integrate();

  if (kelvinDummy.IsInvalid())
    return DEVIATION_INVALID;
  else
    return kelvinDummy.current.r - target;
}
