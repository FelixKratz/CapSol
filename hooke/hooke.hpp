#pragma once

#include "shape.hpp"
#include "laplace.hpp"
#include "config.hpp"
#include "shooting.hpp"
#include "leastSquares.hpp"

#define INITIAL_K 2.5
#define INITIAL_NU 0.5
#define INITIAL_STRESS 0.95


// NOTE: p and tau_s_0 are *not* independent and can not be controlled both
struct HookeParameters {
  double p, g, nu, k, tau_s_0;
  HookeParameters(double g, double k, double tau_s_0 = 1, double p = 0)
    : p(p), nu((k - g )/(k + g)), k(k), tau_s_0(tau_s_0) { };
  HookeParameters() {};
};

struct HookeFitParameters : Operators<HookeFitParameters> {
  double tau_s_0, g, k;
  HookeFitParameters() {};
  HookeFitParameters(HookeParameters p)
    : tau_s_0(p.tau_s_0), g(p.g), k(p.k) {};
};

// These are the independent variables of the differential equations, i.e.
// all those variables we have a differential equation for. Those go to the RK4
// method and are used to integrate the shape.
struct HookeVariables : Operators<HookeVariables> {
  double sd = 0, r = 0, z = 0, psi = 0, tau_s = 0;
};

struct HookeConstitutiveQuantities {
  double lambda_s, lambda_phi, tau_phi, kappa_s, kappa_phi;
};

class Hooke : virtual public Shape,
              public LeastSquares<Hooke, HookeFitParameters>,
              public RK4<HookeVariables>,
              public ShootingMethod<SHOOTING_INTERVALS> {
  int integration_count = 0;
  protected:
  // This field should never be used as a parameter and is only used if the
  // shape equations are non-dimensionalized by a different surface tension
  // than the surface tension of the reference shape e.g. the lower shape of
  // the elastic contact problem. This ratio is defined as \gamma_{reference} /
  // \gamma_{non-dim} and divides the reference surface tension by the tension
  // used to non-dimensionalize the equations.
  double gamma_ratio = 1.0;
  bool use_splines = false;
  bool flip_gravity = false;
  bool enable_wrinkling = true;

  HookeParameters parameters;
  HookeConstitutiveQuantities constitutive_parameters;

  virtual bool ShapeEquations(double s, HookeVariables* current,
                                        HookeVariables* derivative) override;

  void MoveDown() override;
  virtual bool ShootingResult(double optimal_parameter) override;
  virtual double ShootingSingleBoundaryDeviation(double target, double parameter) override;

  public:
  DataSeries tau_s, lambda_s, lambda_phi,
             kappa_phi, kappa_s, tau_phi, ls_times_rs;

  Hooke() : Shape() {};
  Hooke(Hooke* _clone) : Hooke() { Clone(_clone); };
  ~Hooke();

  virtual bool PressureShooting(bool parallel);
  virtual bool PressureShooting(bool parallel, double min, double max);
  void Clone(Hooke* _clone);
  void CloneInto(Hooke* into) override { into->Clone(this); };

  virtual void Clear() override;
  void ResetAdditionalShapeVectors();
  virtual void Push(double s);

  virtual void Integrate(double s_0 = 0., bool reset = true);
  virtual bool Solve(bool async) override { return PressureShooting(async); };

  void Interpolate(bool sparse) override { Interpolate(sparse, false); };
  void Interpolate(bool sparse, bool deformed);
  void MarkInvalid() { invalid = true; Clear(); }
  void SetLaplace(Laplace* laplace) { Clear(); undeformed = laplace; }
  virtual void SetParameters(double p, double nu, double k);
  virtual void SetParameters(HookeParameters params) {
    Clear();
    parameters = params;
  }

  virtual HookeFitParameters GetFitParameters() override;
  virtual HookeFitParameters SetParameters(HookeFitParameters* fitParameters,
                                           bool clampToFit) override;


  void SetPressure(double p) { SetParameters(p, parameters.nu, parameters.k); }
  void SetInitialConditions(double r, double z, double psi, double tau_s);
  void SetApexStress(double tau_s) { SetInitialConditions(0., 0., 0., tau_s); }

  double GetCapillaryStress() { return tau_s.values[tau_s.values.size() - 1]; }
  double GetPressure() { return parameters.p; }
  double GetPoisson() { return parameters.nu; }
  double GetCompression() { return parameters.k; }
  double GetEH0() { return GetCompression() * (2.0 * (1.0 - GetPoisson())); }
  double GetApexStress() { return parameters.tau_s_0; }

  Laplace* undeformed = nullptr;

  bool is_in_wrinkling_domain = false;
  double wrinkling_start = 0.0;
  double wrinkling_stop = 0.0;

  std::pair<double, double> GetXYPairFromS(double s);

  friend class SlipContact;
  friend class ContactPressureShooting;
  friend class ContactVolumeShooting;
  friend class ContactAngleShooting;
};

void SetHookeInitialGuess(double K_init, double nu_init, double tau_s_0_init);
