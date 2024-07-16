#include "pointSet.hpp"
#include <utility>
#include <iostream>
#include "shape.hpp"

bool loadPointSet(std::string filename, PointSet* _pointSet) {
  std::ifstream inputFile(filename);
  if (inputFile.is_open()) {
    while(!inputFile.eof()) {
      double r = 0, z = 0;
      inputFile >> r;
      inputFile >> z;
      if (r == 0 && z == 0) continue;
      std::pair<double, double> xyPair;
      _pointSet->push_back(xyPair);
      _pointSet->back().first = r;
      _pointSet->back().second = z;
    }   
  }
  else {
    return false;
  }
  return true;
}

void addNoiseToPointSet(PointSet* _pointSet, double amount, char mode) {
  std::default_random_engine gen;
  std::normal_distribution<double> g_dist(0.,amount);
  std::uniform_real_distribution<double> u_dist(0., amount);

  double rand_r = 0;
  double rand_z = 0;

  for (int j = 0; j < _pointSet->size(); j++) {
    switch (mode)
    {
      case 'g':
        rand_r = g_dist(gen);
        rand_z = g_dist(gen);
        break;
      case 'u':
        rand_r = u_dist(gen);
        rand_z = u_dist(gen);
        break;
      default:
        std::cout << "Invalid noise mode " << mode << std::endl;
        exit(1);
        break;
    }

    _pointSet->at(j).first += rand_r;
    _pointSet->at(j).second += rand_z;
  }
}

void generatePointSet(Shape* _data, PointSet* _pointSet, bool deformed) {
  // NOTE: the step count for the point set is fixed to 250.

  int stepCount = 250;
  _pointSet->clear();
  double step = (deformed ? _data->sd.values.back()
                          : _data->s.values.back()) / (double(stepCount) - 1.);
  if (deformed) {
    if (_data->sd.values.size() == 0
        || !_data->r.spline
        || !_data->acc
        || !_data->z.spline) {
      printf("Error: Shape has no deformed parametrization\n");
      return;
    }
  }

  // for (int i = 1; i < _data->size; i++) {
  //   double big_step = 100.0 * INTEGRATION_STEP;
  //   bool r_negative = _data->r.values[i] < 0;
  //   bool sd_negative = _data->sd.values[i] < 0;
  //   bool r_big_step = abs(_data->r.values[i - 1] - _data->r.values[i]) > big_step;
  //   bool z_big_step = abs(_data->z.values[i - 1] - _data->z.values[i]) > big_step;
  //   bool sd_big_step = abs(_data->sd.values[i - 1] - _data->sd.values[i]) > big_step;
  //   if (r_negative || sd_negative || r_big_step || z_big_step || sd_big_step) {
  //     printf("Impossible: %d -> %d %d %d %d %d!\n", i, r_negative, sd_negative, r_big_step, z_big_step, sd_big_step);
  //     printf("%d : %f %f %f\n", i, _data->sd.values[i], _data->r.values[i], _data->z.values[i]);
  //     exit(1);
  //   }
  // }

  for (int i = 0; i < stepCount - 1; i += 1) {
    _pointSet->push_back(
        std::make_pair(gsl_spline_eval(_data->r.spline, step*i, _data->acc),
                       gsl_spline_eval(_data->z.spline, step*i, _data->acc)));

    // if (abs(_pointSet->back().first) > 20 || abs(_pointSet->back().second) > 20) {
    //   printf("We have a problem here.. %d %zu %zu %zu\n", _data->IsValid(), _data->r.values.size(), _data->z.values.size(), _data->sd.values.size());
    //   for (int i = _data->sd.values.size() - 40; i < _data->sd.values.size(); i++) {
    //     printf("%d/%zu %f %f %f\n", i, _data->sd.values.size(), _data->sd.values[i], _data->r.values[i], _data->z.values[i]);
    //   }
    //   exit(1);
    // }
  }

  // The last point must always be (0.5, 0), this is the point that connects to
  // the capillary
  _pointSet->push_back(std::make_pair(0.5, 0.));
}
