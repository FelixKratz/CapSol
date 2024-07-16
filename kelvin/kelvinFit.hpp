#pragma once
#include "kelvinSequence.hpp"
#include "pointSet.hpp"

#define ARMA_WARN_LEVEL 1
#include <armadillo>

#define INITIAL_STRESS_AMPLITUDE 0.2
#define INITIAL_ETA 0.25

bool KelvinFit(KelvinSequence* _sequence, std::vector<PointSet>* _points, int count, double* err = NULL);
