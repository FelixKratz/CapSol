#pragma once

#include "hooke.hpp"

// These are the additional parameters the Kelvin shapes need over the Hooke
// shape
struct KelvinParameters {
  double eta, dt;
};

// All of the functionality for the Kelvin shapes is derived from the Hooke
// class, we make surgical changes by replacing some of the virtual functions
// in the Hooke class with new functionality (i.e. the shape equations and 
// the function which calculates the boundary deviations for the shooting
// method).
class Kelvin : public Hooke {
  Kelvin* previous = nullptr;
  DataSeries lambda_s_detailed;
  DataSeries lambda_phi_detailed;
  DataSeries tau_s_detailed;
  DataSeries tau_phi_detailed;
  int integration_count = 0;
  int previousCount = 0;
  bool use_generalized_equation = false;

  KelvinParameters parameters_kelvin;
  double ShootingSingleBoundaryDeviation(double target, double parameter) override;
  bool ShapeEquations(double s, HookeVariables* current,
                                HookeVariables* derivative) override;

  bool ShapeEquationsGeneralized(double s, HookeVariables* current,
                                           HookeVariables* derivative);

  public:
    Kelvin() : Hooke() {};
    Kelvin(Kelvin* _clone) { Clone(_clone); };

    void Clone(Kelvin* _clone);
    virtual HookeFitParameters SetParameters(HookeFitParameters* p, bool clampToFit) override {
      printf("Unimplemented: %s\n", __FUNCTION__);
      return *p;
    };

    virtual void Clear() override {
      Hooke::Clear();
      lambda_s_detailed.Clear();
      lambda_phi_detailed.Clear();
      tau_s_detailed.Clear();
      tau_phi_detailed.Clear();
      integration_count = 0;
    }

    void SetParameters(double p, double nu, double k, double eta, double dt);
    void SetParameters(HookeParameters params) override {
      printf("This should never ever be called..\n");
      exit(1);
    }
    void SetParameters(double p, double nu, double k) override {
      SetParameters(p, nu, k, parameters_kelvin.eta, parameters_kelvin.dt);
    };

    double GetViscosity();
    double GetTimeStep();

    double J_dot_s_s(double t);
    double J_dot_phi_phi(double t);
    double J_dot_s_phi(double t);
    double J_dot_phi_s(double t);

    void SetPreviousShape(Kelvin *data);
};
