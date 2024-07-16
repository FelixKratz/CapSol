#pragma once
#include "shooting.hpp"

class SlipContact;

class ContactVolumeShooting : public ShootingMethod<SHOOTING_INTERVALS> {
  private:
  bool shootingForUpperShape = false;
  SlipContact* contact;

  public:
  ContactVolumeShooting(SlipContact* contact) : contact(contact) {};
  virtual void MarkAsShootingDummy() override;

  bool VolumeShooting(bool parallel);
  virtual bool ShootingResult(double optimal_parameter) override;
  virtual double ShootingSingleBoundaryDeviation(double target,
                                                 double parameter) override;
};
