#pragma once
#define _USE_MATH_DEFINES
#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <iomanip>
#include <time.h>
#include <vector>
#include <string>
#include <sstream>
#include <math.h>
#include <algorithm>
#include "dataSeries.hpp"
#include "config.hpp"

class Residual {
  public:
  double Norm() { return (delta_r * delta_r + delta_z * delta_z); }

  double delta_r;
  double delta_z;
  double s;
};


typedef std::vector< std::pair<double,double> > PointSet;

class Shape {
  protected:
  bool invalid = true;
  bool unphysical = false;

  public:
  virtual ~Shape() { if (acc) gsl_interp_accel_free(acc); }

  DataSeries s, sd, r, z, psi;
  gsl_interp_accel* acc = nullptr;
  uint32_t size = 0;
  double L0;

  virtual void Clear() {
    s.Clear();
    sd.Clear();
    r.Clear();
    z.Clear();
    psi.Clear();
    size = 0;
  };

  void Push(double s_i, double sd_i, double r_i, double z_i, double psi_i) {
    s.Push(s_i);
    sd.Push(sd_i);
    r.Push(r_i);
    z.Push(z_i);
    psi.Push(psi_i);
    size++;
  }
  
  virtual void Pop() {
    if (size < 1) return;

    s.Pop();
    sd.Pop();
    r.Pop();
    z.Pop();
    psi.Pop();
    size--;
  }

  bool IsValid() { return !invalid && !unphysical; }
  bool IsInvalid() { return invalid; }
  bool IsUnphysical() { return unphysical; }

  void MarkValid() { invalid = false; }
  void MarkUnphysical() { unphysical = true; }

  virtual void Interpolate(bool sparse) = 0;
  virtual void MoveDown() = 0;

  double MSError(PointSet *points) {
    int size = (int)points->size();
    double ms = 0.0;
    for (int i = 0; i < size; i++) {
      ms += ResidualForPoint((*points)[i]).Norm();
    }
    return ms / (double)size;
  }

  double RMSError(PointSet *points) {
    return sqrt(MSError(points)*((double)size)/((double) size + 1.));
  }

  double SquareDistance(double s_ref, std::pair<double, double> point) {
    double dz = point.second - gsl_spline_eval(z.spline, s_ref, acc);
    double dr = point.first - gsl_spline_eval(r.spline, s_ref, acc);
    return dr*dr + dz*dz;
  }

  double AbsoluteDistance(double s_ref, std::pair<double, double> point) {
    return sqrt(SquareDistance(s_ref, point));
  }

  Residual ResidualForPoint(std::pair<double,double> point) {
    double pos_plus_dx, pos_minus_dx;
    double pos = 0.5 * L0, step = 0.5 * L0;

    for (;;) {
      step /= 2.0;
      pos_plus_dx = SquareDistance(pos + step, point);
      pos_minus_dx = SquareDistance(pos - step, point);

      if (pos_plus_dx < pos_minus_dx) pos += step;
      else if (pos_plus_dx > pos_minus_dx) pos -= step;
      else break;
    }

    return {
      .delta_r = point.first - gsl_spline_eval(r.spline, pos, acc),
      .delta_z = point.second - gsl_spline_eval(z.spline, pos, acc),
      .s = pos
    };
  }

  double Volume() {
    if (z.values.size() < 1 || sd.values.size() < 1) return 0;

    // Prepare the integrand for the volume calculation
    DataSeries r2SinPsi;
    for (int i = 0; i < sd.values.size(); i++) {
      r2SinPsi.Push(r.values[i] * r.values[i] * sin(psi.values[i]));
    }

    r2SinPsi.Interpolate(&sd);

    double volume = gsl_spline_eval_integ(r2SinPsi.spline,
                                          0,
                                          sd.values.back(),
                                          acc              ) * M_PI;

    return volume;
  }
};
