#pragma once
#include "pointSet.hpp"
#include "rk4.hpp"
#include "leastSquares.hpp"

#define LAPLACE_FIT_MAX_ITERATIONS 100
#define INITIAL_P_L 3.5
#define INITIAL_RHO 1.0

struct LaplaceParameters {
  double p_L = 0, rho = 0;
  int omega = 2;
  double capillary_radius = 0.5;

  LaplaceParameters(double p_L, double rho, int omega = 2,
                                            double capillary_radius = 0.5)
    : p_L(p_L), rho(rho), omega(omega), capillary_radius(capillary_radius) { };
  LaplaceParameters() {};
};

struct LaplaceFitParameters : Operators<LaplaceFitParameters> {
  double p_L, rho;
  LaplaceFitParameters() {};
  LaplaceFitParameters(LaplaceParameters p) : p_L(p.p_L), rho(p.rho) {}
};

struct LaplaceVariables : Operators<LaplaceVariables> {
  double r = 0, z = 0, psi = 0, s = 0;
};

class Laplace : virtual public Shape,
                virtual public LeastSquares<Laplace, LaplaceFitParameters>,
                public RK4<LaplaceVariables> {
  private:
  int integration_count = 0;
  LaplaceParameters parameters;

  bool ShapeEquations(double s, LaplaceVariables* current,
                                LaplaceVariables* derivative) override;

  void MarkInvalid() { invalid = true; Clear(); };

  void MoveDown() override;
  void Interpolate(bool sparse = false) override;

  public:
  DataSeries r_integration;
  DataSeries kappa_s_integration;
  DataSeries psi_integration;

  Laplace(double p = 0, double rho = 0, int omega = 2) {
    SetParameters(p, rho, omega);
  }

  void Clone(Laplace* _clone);
  void CloneInto(Laplace* into) override { into->Clone(this); };
  void Clear() override {
    Shape::Clear();
    if (fitting_dummy) return;

    r_integration.Clear();
    kappa_s_integration.Clear();
    psi_integration.Clear();
    integration_count = 0;
  }

  void Pop() override {
    Shape::Pop();
    if (fitting_dummy) return;
    if (r_integration.values.size() < 4) return;

    // If we pop a value from the integration, we need to pop four values from
    // the arrays tracing the intermediate integration steps (four intermediate
    // steps in RK4)
    while (r_integration.values.size() > integration_count - 4) {
      r_integration.Pop();
      kappa_s_integration.Pop();
      psi_integration.Pop();
    }
    integration_count -= 4;
  }

  LaplaceFitParameters GetFitParameters() override;
  LaplaceFitParameters SetParameters(LaplaceFitParameters* fitParameters,
                                     bool clampToFit) override;

  void SetParameters(double p, double rho, int omega = 2);
  void SetParameters(LaplaceParameters params) {
    SetParameters(params.p_L, params.rho, params.omega);
  };
  void SetCapillaryRadius(double r) { parameters.capillary_radius = r; }

  double GetPressure() { return parameters.p_L; }
  double GetDensity() { return parameters.rho; }
  int GetOmega() { return parameters.omega; }
  double CalculateArea();

  bool Solve(bool parallel = false) override;
  friend class SlipContact;
};

void SetLaplaceInitialGuess(double p_L, double rho);
