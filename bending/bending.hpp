#pragma once
#include "hooke.hpp"
#include "pressureShooting.hpp"
#include "volumeShooting.hpp"

// This file will integrate the shape for a shell with bending stiffness
// pressed against a solid wall with an external force. The shell shall be a
// half-sphere initially, clamped at the top (i.e. \Psi(L) = \pi / 2).

struct BendingParameters {
  // The length of the contacting region in *deformed* arc length coordinates
  // NOTE: Only the contact_length can be specified, the force may be
  // calculated at the splitting point
  double contact_length = 0, force = 0;
  double k, e_b, nu, p;
  double tau_s_0, q_plus = 0;
  double gamma_contact = 0, gamma_wall = GAMMA_SCALE;

  BendingParameters() {};
  BendingParameters(double contact_length, double force, double k, double e_b, double nu, double p, double tau_s_0, double q_plus)
    : contact_length(contact_length),
      force(force),
      k(k),
      e_b(e_b),
      nu(nu),
      p(p),
      tau_s_0(tau_s_0),
      q_plus(q_plus) {};

  void Print() {
    for (int i = 0; i < sizeof(BendingParameters) / sizeof(double); i++) {
      printf("%f ", *((double*)this + i));
    }
    printf("\n");
  }
};

struct BendingConstitutiveParameters {
  double lambda_s, lambda_phi, tau_phi, K_s, K_phi, m_phi, kappa_s, kappa_phi;
};

// NOTE: These are the integration variables for the free regions of the
// solution.
struct BendingVariables : Operators<BendingVariables> {
  double s = 0, r = 0, z = 0, psi = 0, m_s = 0, tau_s = 0, q = 0;
};

class Bending : virtual public Shape, public RK4<BendingVariables> {
  private:
  int integration_count = 0;
  bool contact_region = true;
  Laplace* undeformed = NULL;
  BendingParameters parameters;
  BendingConstitutiveParameters constitutive_parameters;
  void Clone(Bending* clone);

  virtual bool ShapeEquations(double s_0, BendingVariables* current, BendingVariables* derivative) override;
  virtual bool SplittingPoint(double s_0, BendingVariables* current);

  public:
  Bending() {};

  virtual void Clear() override { Shape::Clear(); integration_count = 0; };
  void SetLaplace(Laplace* laplace) { undeformed = laplace; };
  void SetParameters(BendingParameters params);
  void SetInitialConditions(double s, double r, double z, double psi, double tau_s_0);
  void SetApexStress(double tau_s_0);

  void MarkInvalid() { invalid = true; Clear(); };
  void Integrate();

  virtual void Interpolate(bool sparse) override { };
  virtual void MoveDown() override { };
};
