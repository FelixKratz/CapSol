#pragma once
#include <iostream>
#include <cmath>
#include "config.hpp"

#define RANGE_INVALID -2
#define RANGE_FULL -1

#define DEVIATION_INVALID -1e10

template<int m> class ShootingMethod;
template<int m> struct shooting_context {
  double target, parameter;
  double deviation;
  ShootingMethod<m>* method;
};

template<int m>
void* shooting_proc(void* context) {
  struct shooting_context<m>* shooting = (struct shooting_context<m>*)context;
  shooting->deviation =
        shooting->method->ShootingSingleBoundaryDeviation(shooting->target,
                                                          shooting->parameter);
  return NULL;
}

struct ShootingRange {
  bool parallel;
  int shifts;
  double target;
  double min;
  double max;
  int min_id;
  int max_id;
  bool range_changed;
  bool machine_precision_reached;
  double minimum_precision;

  ShootingRange(double target, double minimum, double maximum, bool parallel = false, double minimum_precision = 1e-14, int shifts = 0)
    : parallel(parallel),
      shifts(shifts),
      target(target),
      min(minimum),
      max(maximum),
      min_id(0),
      max_id(-1),
      range_changed(false),
      machine_precision_reached(false),
      minimum_precision(minimum_precision) { }
};

template<int m>
class ShootingMethod {
  private:
  double deviation_min;
  double parameter_best;

  inline bool SearchForValidRange(ShootingRange* range) {
    int firstValidIndex = RANGE_INVALID;
    int lastValidIndex = RANGE_FULL;
  
    for (int i = 0; i < m - 1; i++) {
      if (deviations[i] > DEVIATION_INVALID
          && firstValidIndex == RANGE_INVALID)
        firstValidIndex = i - 1;
      else if (deviations[i] <= DEVIATION_INVALID
               && firstValidIndex != RANGE_INVALID
               && lastValidIndex == RANGE_FULL) {
        lastValidIndex = i;
      }
    }
    if (firstValidIndex == RANGE_INVALID && lastValidIndex == RANGE_FULL) {
      firstValidIndex = RANGE_INVALID;
      lastValidIndex = RANGE_INVALID;
      return false;
    }
  
    // Check Right border
    if (lastValidIndex > 0 && firstValidIndex < 0) {
      firstValidIndex = lastValidIndex - 1;
    }
    // Check Left Border
    if (lastValidIndex < 0 && firstValidIndex > 0) {
      lastValidIndex = firstValidIndex + 1;
    }
  
    // Set Full Range if no invalid regions found
    if (lastValidIndex == RANGE_FULL) lastValidIndex = m - 1;
    if (firstValidIndex == RANGE_INVALID) firstValidIndex = RANGE_FULL;
    range->range_changed = false;

    if ((lastValidIndex - firstValidIndex <= m - 2)
        && (lastValidIndex - firstValidIndex > 0)) {
      range->min = parameters[firstValidIndex < 0 ? 0 : firstValidIndex];
      range->max = parameters[lastValidIndex > m - 1 ? m - 1 : lastValidIndex];
      range->range_changed = true;
    }

    range->min_id = firstValidIndex;
    range->max_id = lastValidIndex;
    return true;
  }

  inline void SearchForZeroCrossing(ShootingRange* currentRange) {
    for (int i = 0; i < m - 1; i++) {
      if (((deviations[i] <= 0 && deviations[i+1] >= 0)
            || (deviations[i] >= 0 && deviations[i+1] <= 0))
          && (deviations[i] > DEVIATION_INVALID
              && deviations[i+1] > DEVIATION_INVALID)) {

        currentRange->min_id = i;
        currentRange->max_id = i+1;
        currentRange->min = parameters[i];
        currentRange->max = parameters[i+1];
        
        currentRange->range_changed = true;
        break;
      }
    }

    // If the range has not changed, we did not find a zero crossing, the zero
    // crossing might be in between the first invalid point and the first valid
    // point or between the last valid point and the neighbouring invalid
    // point. We prioritize the lower deviation next to the invalid point.
    if (!currentRange->range_changed) {
      int first_valid_id = 0;
      int last_valid_id = m - 1;
      for (int i = 0; i < m - 1; i++) {
        if ((deviations[i] <= DEVIATION_INVALID)
          && (deviations[i+1] > DEVIATION_INVALID)) {
          first_valid_id = i + 1;
          break;
        }
      }
      for (int i = m - 2; i >= 0; i--) {
        if ((deviations[i+1] <= DEVIATION_INVALID)
            && (deviations[i] > DEVIATION_INVALID)) {
          last_valid_id = i;
          break;
        }
      }

      if (((first_valid_id > 0 && last_valid_id < m - 1)
          && (abs(deviations[first_valid_id]) < abs(deviations[last_valid_id])))
          || (first_valid_id > 0 && last_valid_id >= m - 1)) {
        currentRange->min_id = first_valid_id - 1;
        currentRange->max_id = first_valid_id;
        currentRange->min = parameters[currentRange->min_id];
        currentRange->max = parameters[currentRange->max_id];
        
        currentRange->range_changed = true;
      } else if ((first_valid_id > 0 && last_valid_id < m - 1)
                || ((last_valid_id < m - 1) && first_valid_id <= 0)) {
        currentRange->min_id = last_valid_id;
        currentRange->max_id = last_valid_id + 1;
        currentRange->min = parameters[currentRange->min_id];
        currentRange->max = parameters[currentRange->max_id];
        
        currentRange->range_changed = true;
      }
    }

    // HACK: This detects numerical elimination and will be true if the
    // interval is wider than the numerical resolution with respect to the
    // unit scale and false if it is smaller.
    currentRange->machine_precision_reached
                      = (1.0 - (currentRange->max - currentRange->min)) == 1.0;
  }

  inline bool ShiftRange(ShootingRange* currentRange) {
    // Set an arbitrary limit to shifting to not get stuck.
    if (currentRange->shifts < 10) {
    // Shifting only makes sense if the borders of the shooting region exist.
      if (fabs(deviations[currentRange->max_id - 1])
          < fabs(deviations[currentRange->min_id + 1])) {
        // Right Shifting
        ShootingRange rightRange(currentRange->target,
                                 parameters[currentRange->max_id - 1],
                                 2.*parameters[currentRange->max_id - 1]
                                 - currentRange->min,
                                 currentRange->parallel,
                                 currentRange->minimum_precision,
                                 currentRange->shifts + 1);

        if (Shoot(rightRange)) return true;
      } else {
        // Left Shifting
        ShootingRange leftRange(currentRange->target,
                                2.*parameters[currentRange->min_id + 1]
                                - currentRange->max,
                                parameters[currentRange->min_id + 1],
                                currentRange->parallel,
                                currentRange->minimum_precision,
                                currentRange->shifts + 1);

        if (Shoot(leftRange)) return true;
      }
    }
    return false;
  }

  inline bool SecantMethod(ShootingRange* range) {
    double x_n_minus_2 = range->min;
    double x_n_minus_1 = range->max;
    double f_x_n_minus_2 = deviations[range->min_id];
    double f_x_n_minus_1 = deviations[range->max_id];

    double difference, norm;
    if (f_x_n_minus_1 == DEVIATION_INVALID
        || f_x_n_minus_2 == DEVIATION_INVALID) return false;

    int max_iterations = 30;
    for(;;) {
      difference = x_n_minus_2 * f_x_n_minus_1 - x_n_minus_1 * f_x_n_minus_2;
      norm = f_x_n_minus_1 - f_x_n_minus_2;

      x_n_minus_2 = x_n_minus_1;
      f_x_n_minus_2 = f_x_n_minus_1;
      x_n_minus_1 = difference / norm;
      f_x_n_minus_1 = ShootingSingleBoundaryDeviation(range->target,
                                                      x_n_minus_1);

      if ((max_iterations-- < 0) || (f_x_n_minus_1 == DEVIATION_INVALID))
        return false;

      double abs_f_x_n_minus_1 = fabs(f_x_n_minus_1);

      if (abs_f_x_n_minus_1 < deviation_min) {
        deviation_min = abs_f_x_n_minus_1;
        parameter_best = x_n_minus_1;
      } else if (deviation_min < range->minimum_precision) {
        // In this case we have hit the machine precision and can not
        // improve further, but are already better than the minimum
        // required precision
        break;
      } else {
        // The error is not decreasing and it is not already sufficiently
        // small, the secant method search has failed.
        return false;
      }
    }

    return true;
  }

  inline void ShootingBoundaryDeviations(double target, bool parallel) {
    if (parallel) {
      pthread_t threads[num_sections];
      struct shooting_context<m> contexts[num_sections];

      for (int i = 0; i < num_sections; i++) {
        contexts[i].parameter = parameters[i];
        contexts[i].target = target;
        contexts[i].method = this;

        pthread_create(&threads[i], NULL, shooting_proc<m>, &contexts[i]);
      }
      for (int i = 0; i < num_sections; i++) {
        pthread_join(threads[i], NULL);
        deviations[i] = contexts[i].deviation;
      }
    } else {
      for (int i = 0; i < num_sections; i++) {
        deviations[i] = ShootingSingleBoundaryDeviation(target,
                                                        parameters[i]);
      }
    }
  }

  protected:
  bool shooting_dummy = false;
  double deviations[m];
  double parameters[m];
  static const int num_sections = m;

  public:
  virtual void MarkAsShootingDummy() { shooting_dummy = true; }
  inline virtual bool ShootingResult(double optimal_parameter) = 0;
  inline virtual double ShootingSingleBoundaryDeviation(double target, double parameter) = 0;

  inline bool Shoot(ShootingRange range) {
    static const int intervals = m;
    double h = (range.max - range.min) / ((double)intervals - 1.0);
    bool disableSecantMethod = false;
    int counter = 0;

    for(;;) {
      if (counter++ > 1000) return false;

      range.range_changed = false;
      range.machine_precision_reached = false;

      h = (range.max - range.min) / ((double)intervals - 1.0);

      // The maximum precision we are able to achieve is prescribed by the
      // length of the mantissa, which is 53bit wide. This makes
      // 53log10(2) significant digits available to our calculation.
      // The maxiumum precision is thus given by
      // 10^(-53log10(2)) * range.min ~ 1e-16 * range.min
      if (h < 1e-16*range.min) break;

      for (int i = 0; i < intervals; i++) parameters[i] = range.min + i * h;

      // This is implemented in the parent class
      ShootingBoundaryDeviations(range.target, range.parallel);

      // Find minimal boundary deviation
      deviation_min = fabs(deviations[0]);
      parameter_best = parameters[0];
      for (int i = 1; i < intervals; i++) {
        double absolute_deviation = fabs(deviations[i]);
        if (absolute_deviation < deviation_min) {
          deviation_min = absolute_deviation;
          parameter_best = parameters[i];
        }
      }

      if (deviation_min < range.minimum_precision) break;

      SearchForZeroCrossing(&range);
      if (range.machine_precision_reached) break;

      if (!range.range_changed) {
        // Search for a valid subrange in the shooting interval, i.e. where
        // solutions to the problem exist.
        if (!SearchForValidRange(&range)) return false;

        // If the current range has changed and we did not zoom too many times
        // already, we go back to the top of the loop and evaluate the new
        // range again to hopefully find the boundary crossing (it might be
        // between the first/last valid point and the neighbouring invalid
        // point).
        if (range.range_changed) continue;

        // When we are in a valid subrange and still do not find a
        // zero-crossing, it might be to the left/right of our current
        // interval. We check this by shifting the target range in the
        // opposite direction of the absolute deviation gradient.
        // TODO: Fix cyclic range shifts by disallowing shifting to the
        // previous range
        if (ShiftRange(&range)) return true;

        // Everything failed, we have still not found a zero-crossing.
        return false;
      }

      // The multisection method will bring us in the vincinity
      // of the zero crossing, from there we try to converge even faster
      if (!disableSecantMethod) {
        // Now that we have found a reasonable range we
        // return to a secant method because of its superior convergence.
        if (!SecantMethod(&range)) {
          // If the secant method fails (because there is an invalid shape,
          // or because it did not reach the target precision)
          // we simply fall back to the multisection search
          disableSecantMethod = true;
          continue;
        }

        break;
      }
    }

    return ShootingResult(parameter_best);
  }
};
