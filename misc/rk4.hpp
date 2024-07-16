#pragma once
#include "operators.hpp"

template<class Variables>
class RK4 {
  private:
  virtual bool ShapeEquations(double s, Variables* current,
                                        Variables* derivative) = 0;

  Variables intermediate;
  Variables rk1, rk2, rk3, rk4;
  protected:
  Variables initial_conditions;
  Variables current;
  double integration_step = 1e5;

  bool RungeKutta(double s) {
    bool success = true;
    success &= ShapeEquations(s, &current, &rk1);
    intermediate = current + rk1*0.5*integration_step;
    success &= ShapeEquations(s + 0.5*integration_step, &intermediate, &rk2);
    intermediate = current + rk2*0.5*integration_step;
    success &= ShapeEquations(s + 0.5*integration_step, &intermediate, &rk3);
    intermediate = current + rk3*integration_step;
    success &= ShapeEquations(s + integration_step, &intermediate, &rk4);

    current += (rk1 + rk2*2.0 + rk3*2.0 + rk4)*integration_step*(1.0/6.0);
    return success;
  }
  public:

  void SetIntegrationStep(double step) { integration_step = step; }
};
