#pragma once

#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>

#include "config.hpp"
#include "shape.hpp"
#include "kelvinFit.hpp"
#include "shapeEquations.hpp"

#include <algorithm>
#include <time.h>
#include <cmath>
#include <utility>
#include <math.h>
#include <chrono>

template <typename T>
std::string to_string(T value) {
  std::ostringstream buffer;
  buffer << value;
  return buffer.str();
}

struct laplace_parameters {
  double p_L;
  double rho;
};

struct hooke_parameters {
  double p_L;
  double rho;
  double K;
  double Psi;
  double tau_s_0;
};

struct kelvin_parameters {
  double p_L;
  double rho;
  double K;
  double Psi;
  double dtau_s_0;
  double eta;
  double dt;
};

bool fitLaplaceShape(PointSet* laplacePointSet, Laplace* laplaceShape, struct laplace_parameters* fit_parameters, double* err = NULL);
bool fitLaplaceShape(PointSet* laplacePointSet, struct laplace_parameters* fit_parameters, double* err = NULL);
bool fitLaplaceShapes(std::vector<PointSet*> laplacePointSets, Laplace* laplace, struct laplace_parameters* fit_parameters, double* err);
bool fitLaplaceShapes(std::vector<PointSet*> laplacePointSets, struct laplace_parameters* fit_parameters, double* err);


bool fitHookeShape(Laplace* laplaceShape, PointSet* hookePointSet, struct hooke_parameters* fit_parameters, double* err = NULL);
bool fitHookeShapes(Laplace* laplace, std::vector<PointSet*> hookePointSet, struct hooke_parameters* fit_parameters, double* err);


bool fitKelvinSequence(std::vector<PointSet>* kelvinPointSet, Laplace* laplaceShape, struct kelvin_parameters* fit_parameters, int count, int periods = 1, double* err = NULL);
