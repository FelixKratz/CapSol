#include "hooke.hpp"
#include "pressureShooting.hpp"
#include "volumeShooting.hpp"
#include "angleShooting.hpp"

struct SlipContactParameters {
  HookeParameters hooke_u_parameters;
  HookeParameters hooke_d_parameters;

  // The length of the contacting region in *deformed* arc length coordinates
  // NOTE: Only the contact_length can be specified, the force may be
  // calculated at the splitting point
  double contact_length;
  double force = 0;

  // The ratio of surface tensions
  // (gamma_i_d + gamma_o_d) / (gamma_i_u + gamma_o_u)
  double tension_ratio;

  // This gives (gamma_i_u) / (gamma_i_u + gamma_o_u)
  double tension_inner_u;

  // This gives (gamma_i_d) / (gamma_i_d + gamma_o_d)
  double tension_inner_d;

  // This is the possible adhesion tension (gamma_ud) / (gamma_i_u + gamma_o_u)
  double tension_ud;

  SlipContactParameters(HookeParameters p_hooke_u, HookeParameters p_hooke_d,
                                                   double contact_length,
                                                   double tension_ratio,
                                                   double tension_inner_u,
                                                   double tension_inner_d,
                                                   double tension_ud = 0.)
    : hooke_u_parameters(p_hooke_u),
      hooke_d_parameters(p_hooke_d),
      contact_length(contact_length),
      tension_ratio(tension_ratio),
      tension_inner_u(tension_inner_u),
      tension_inner_d(tension_inner_d),
      tension_ud(tension_ud) { };
  SlipContactParameters() {};

  void Print() {
    for (int i = 0; i < sizeof(SlipContactParameters) / sizeof(double); i++) {
      printf("%f ", *((double*)this + i));
    }
    printf("\n");
  }
};

// NOTE: These are the integration variables for the free *slip* case.
// In the no-slip case, only tau_s_u + tau_s_d appear
struct SlipContactVariables : Operators<SlipContactVariables> {
  double s = 0, r = 0, z = 0, psi = 0, tau_s_u = 0, tau_s_d = 0, s_0_d = 0.;
};

class SlipContact : public RK4<SlipContactVariables> {
  private:
  // This class is declared as a friend to Hooke, such that we are allowed to
  // access its private members, here the book keeping is performed
  Hooke hooke_u, hooke_d;
  ContactPressureShooting pressure_cannon;
  ContactVolumeShooting volume_cannon;
  ContactAngleShooting angle_cannon;

  // HACK: Calculate symmetric configurations faster (no double shooting)
  bool symmetric = false;

  SlipContactParameters parameters;
  void Clone(SlipContact* clone);
  void IntegrateLiquid();

  virtual bool ShapeEquations(double s_0, SlipContactVariables* current,
                                          SlipContactVariables* derivative) override;

  virtual bool SplittingPoint(double s_0, SlipContactVariables* current);
  void SyncAndPush(double s);

  public:
  SlipContact() : pressure_cannon(this), volume_cannon(this), angle_cannon(this) {};

  void SetLaplace(Laplace* laplace_u, Laplace* laplace_d);
  void SetParameters(SlipContactParameters params);
  void SetInitialConditions(double s, double r, double z, double psi, double tau_s_0_u, double tau_s_0_d);
  void SetApexStresses(double tau_s_0_u, double tau_s_0_d);

  void MarkAsPressureShootingDummy() { pressure_cannon.MarkAsShootingDummy(); }
  void MarkSymmetric() { symmetric = true; }
  void Integrate();

  // Searches for the pressure required to reach the capillary at s_0 = L_0
  bool PressureShooting(bool parallel);

  // Searches for the apex stresses to achieve the reference shape volume
  // and the pressure required to reach the capillary
  bool VolumeShooting(bool parallel);

  // Searches for the pressure required to reach Psi = pi at s_0 = L_0/2
  bool AngleShooting(bool parallel);
  double GetUpperVolume() { return hooke_u.Volume(); }
  double GetLowerVolume() { return hooke_d.Volume(); }

  friend class ContactAPI;
  friend class ContactUnitCellAPI;
  friend class ContactPressureShooting;
  friend class ContactVolumeShooting;
  friend class ContactAngleShooting;
};
