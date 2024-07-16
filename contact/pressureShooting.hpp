#pragma once
#include "shooting.hpp"

class SlipContact;

class ContactPressureShooting : public ShootingMethod<SHOOTING_INTERVALS> {
  private:
  bool shootingForUpperShape = false;
  SlipContact* contact;

  public:
  ContactPressureShooting(SlipContact* contact) : contact(contact) {};
  virtual void MarkAsShootingDummy() override;

  bool PressureShooting(bool parallel);
  virtual bool ShootingResult(double optimal_parameter) override;
  virtual double ShootingSingleBoundaryDeviation(double target,
                                                 double parameter) override;
};
