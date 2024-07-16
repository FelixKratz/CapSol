#pragma once
#include "shooting.hpp"

class SlipContact;

class ContactAngleShooting : public ShootingMethod<SHOOTING_INTERVALS> {
  private:
  bool shootingForUpperShape = false;
  SlipContact* contact;

  public:
  ContactAngleShooting(SlipContact* contact) : contact(contact) {};
  virtual void MarkAsShootingDummy() override;

  bool AngleShooting(bool parallel);
  virtual bool ShootingResult(double optimal_parameter) override;
  virtual double ShootingSingleBoundaryDeviation(double target,
                                                 double parameter) override;
};
