#pragma once
// Disable the warning when a system is singular and an approx solution is attempted
#define ARMA_WARN_LEVEL 1
#include <armadillo>
#include "config.hpp"
#include "shape.hpp"
#include "operators.hpp"

template<class Object>
void* solve_proc(void* context) {
  return (void*)((Object*)context)->Solve(false);
}

template<class Object, class Parameters>
class LeastSquares : virtual public Shape {
  bool LeastSquaresIteration(std::vector<PointSet*> pointSets, double* dE, double* err) {
    Parameters p = GetFitParameters();

    int pointCount = 0;
    for (auto set : pointSets) { pointCount += set->size(); }

    int d = sizeof(p) / sizeof(double);
    bool converged = false;

    arma::mat jacobian(2*pointCount, d);
    arma::vec residual(2*pointCount);
    arma::vec update(d);
    arma::vec parameters(d);

    for (int i = 0; i < d; i++) parameters[i] = p[i];

    std::pair<double,double>* point;
    Residual error, error_right, error_left;
    double err_prev = 0;
    for (auto set : pointSets) err_prev += MSError(set);
    err_prev = sqrt(err_prev / pointSets.size());

    double err_next = 1e10;

    #ifdef WATCH_LEAST_SQUARES
      std::cout << std::endl << "_____ FIT INTERATION _____" << std::endl;
      for (int i = 0; i < d; i++) { std::cout << parameters[i] << std::endl; }
      std::cout << std::endl;
    #endif

    bool valid_derivatives[d];
    pthread_t threads[d];
    Object object_local[d];
    double dx[d];

    for (int i = 0; i < d; i++) {
      dx[i] = DX_FIT;

      // Calculate the varied parameter set for each fitting parameter
      Parameters parameters_local;
      for (int j = 0; j < d; j++) {
        parameters_local[j] = parameters[j] + (i == j) * dx[i];
      }

      // Clone the fit object into a local dummy and set the varied parameter
      // set, used to calculate the derivatives
      CloneInto(&object_local[i]);
      Parameters valid_parameters
                      = object_local[i].SetParameters(&parameters_local, true);

      dx[i] = valid_parameters[i] - parameters[i];

      // The parameter space could be constrained here, such that we need to
      // sample the negative direction to calculate a gradient.
      if (dx[i] < 1e-10) {
        dx[i] = -DX_FIT;
        for (int j = 0; j < d; j++) {
          parameters_local[j] = parameters[j] + (i == j) * dx[i];
          object_local[i].SetParameters(&parameters_local, true);
        }
      }
      object_local[i].MarkAsFittingDummy();

      // Spawn a thread to generate the solution
      pthread_create(&threads[i], NULL, solve_proc<Object>, &object_local[i]);
    }

    bool valid = true;
    for (int i = 0; i < d; i++) {
      // Join the worker threads and wait until they are finished calculating
      // their shapes
      pthread_join(threads[i], (void**)&valid_derivatives[i]);
      valid &= valid_derivatives[i];
    }

    if (!valid) {
      #ifdef WATCH_LEAST_SQUARES
        std::cout << "One or more derivatives became invalid... Trying to continue." << std::endl;
      #endif

      valid = true;
      for (int i = 0; i < d; i++) {
        if (!valid_derivatives[i]) {
          // If there is no solution in a given direction we simply set the
          // original parameter set to calculate the derivatives
          // object_local[i].SetParameters(&p, false);
          // valid &= object_local[i].Solve(true);
          dx[i] = 0;
        }
      }
    }

    // if (!valid) {
    //   printf("Previous solution invalid...should be impossible, abort\n");
    //   return true;
    // }

    // Calculate the Jacobian Matrix and error vector
    int set_count = 0;
    for (int j = 0; j < pointSets.size(); j++) {
      PointSet* points = pointSets[j];
      for (int y = 0; y < points->size(); y++) {
        point = &(*points)[y];
        error = ResidualForPoint(*point);
        residual[set_count + y] = - error.delta_r;
        residual[set_count + y + pointCount] = - error.delta_z;
        for (int i = 0; i < d; i++) {
          double offset = dx[i];

          // The offset can actually be zero (very unlikely), because we might
          // need to revert both the left and the right step due to convergence
          // problems.
          if (!offset) {
            jacobian(set_count + y, i) = 0;
            jacobian(set_count + y + pointCount, i) = 0;
          } else {
            error_right = object_local[i].ResidualForPoint(*point);
            error_left = error;

            double dr = error_right.delta_r - error_left.delta_r;
            double dz = error_right.delta_z - error_left.delta_z;

            jacobian(set_count + y, i) = dr / offset;
            jacobian(set_count + y + pointCount, i) = dz / offset;
          }
        }
      }
      set_count += points->size();
    }

    // Invert the Jacobian and calculate the parameter update

    arma::solve(update, jacobian, residual);
    // Alternatively (but a bit slower and less precise):
    // update = arma::pinv(jacobian) * residual;

    // We limit the maxiumum step size to 1.0, where we calculate a new update
    // to mitigate huge parameter jumps in shallow parameter regions
    double norm = arma::norm(update, 2);
    if (norm > 1.0) update /= norm;

    double xi = 0.5;

    // Backtrace through the line connecting the new parameter point with the
    // current parameter point, until a solution is found that actually reduces
    // the error.

    for (;;) {
      Parameters local;
      for (int j = 0; j < d; j++) local[j] = parameters[j] + update[j];
      SetParameters(&local, true);
      bool success = Solve(true);

      if (success) {
        double err_next = 0;
        for (auto set : pointSets) err_next += MSError(set);
        err_next = sqrt(err_next / pointSets.size());
        *dE = fabs(err_next - err_prev);
        if (norm < EPS_NEWTON_HOOKE) {
          if (err_next <= err_prev) {
            converged = true;
          }
          else {
            SetParameters(&p, true);
            if (Solve(true)) {
              converged = true;
              err_next = err_prev;
            }
          }
          break;
        } else if (err_next <= err_prev) break;
      } else if (norm < EPS_NEWTON_HOOKE) return true;

      #ifdef WATCH_LEAST_SQUARES
        std::cout << std::setw(20) << "shift" << std::endl;
        for (int i = 0; i < d; i++) std::cout << std::setw(20) << update[i];
        std::cout << std::endl;
      #endif

      // The new parameter set did not give a better fit -- or a solution did
      // not exist => scale the update vector and try again...
      update *= xi;
      norm = arma::norm(update, 2);
    }

    if (!converged) {
      // Create the parameter set for the update
      Parameters local;
      for (int j = 0; j < d; j++) local[j] = parameters[j] + update[j];

      // Try to set the updated parameter set -- constraints are applied in
      // the SetParameters function
      SetParameters(&local, true);
      Solve(true);
    }

    if (err) *err = err_next;
    return converged;
  }

  protected:
  bool fitting_dummy = false;

  public:
  bool Fit(std::vector<PointSet*> pointSets, double* error_out) {
    if (error_out) *error_out = 1e10;
    bool converged = false;

    // This keeps track of the last fitting error improvements
    int errorCache = 5;
    double prevErrors[errorCache];
    for (int i = 0; i < errorCache; i++) prevErrors[i] = 1e10;

    MarkAsFittingDummy();
    int counter = 0;
    while (!converged) {
      converged = LeastSquaresIteration(pointSets,
                                        &prevErrors[counter % errorCache],
                                        error_out);
      if (counter++ > 50) {
        #ifdef WATCH_LEAST_SQUARES
          std::cout << "Limiting Iterations... Converged." << std::endl;
        #endif
        converged = true;
        break;
      }

      // If the slope of the error landscape is too shallow we abort the fit
      double dE = 0;
      for (int i = 0; i < errorCache; i++) { dE += prevErrors[i]; }
      if (dE < 1e-8) {
        #ifdef WATCH_LEAST_SQUARES
          std::cout << "Error not changing any more... Converged." << std::endl;
        #endif
        converged = true;
      }
    }

    if (converged && !invalid) {
      UnmarkAsFittingDummy();
      Interpolate(false);
    }

    return converged;
  }

  void MarkAsFittingDummy() { fitting_dummy = true; }
  void UnmarkAsFittingDummy() { fitting_dummy = false; }

  virtual void CloneInto(Object* clone) = 0;
  virtual bool Solve(bool parallel = false) = 0;
  virtual Parameters SetParameters(Parameters* fit_parameters, bool clampToFit) = 0;
  virtual Parameters GetFitParameters() = 0;
};
