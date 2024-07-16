#include "kelvinFit.hpp"

double stress_amplitude_init = INITIAL_STRESS_AMPLITUDE;
double eta_init = INITIAL_ETA;

void* fit_kelvin_proc(void* context) {
  return (void*)((KelvinSequence*) context)->Solve();
}

bool KelvinFitIteration(KelvinSequence* _sequence, std::vector<PointSet>* _points, int count, double* dE, double* err) {
  // The points of all shapes in the sequence need to be considered here
  int N = 0;
  for (int i = 0; i < count; i++) {
    N += (*_points)[i].size();
  }

  const int nparams = 3;
  double parameters[] = { _sequence->shapes[0].GetCompression(),
                          _sequence->shapes[0].GetViscosity(),
                          _sequence->stressAmplitude             };

  arma::mat _J(2*N, nparams);
  arma::vec _F(2*N);
  arma::vec _D(nparams);
  arma::vec _PI(nparams);

  std::pair<double,double> *p;
  Residual error, error_x1, error_x2;
  double xi = 0.5;
  bool continue_linesearch = true;
  bool converged = false;
  double err_prev = 0, err_next = 0;

  for (int i = 0; i < nparams; i++) _PI[i] = parameters[i];

  KelvinSequence c[2*nparams];
  for (int i = 0; i < 2*nparams; i++) {
    c[i].SetCount(count);
  } 

  #ifdef WATCH_KELVIN_FITTING
    std::cout << std::setw(20) << "PI(K,eta, dtau_s_0) ->";
    for (int i = 0; i < nparams; i++) std::cout << std::setw(20) << std::setprecision(10) << _PI[i];
    std::cout << " Error: " << err_prev << std::endl;
  #endif

  for (int j = 0; j < count; j++) {
    err_prev += _sequence->shapes[j].MSError(&(*_points)[j]);
  }
  err_prev = sqrt(err_prev / (double)count);

  bool valid_derivatives[2*nparams];

  #pragma omp parallel for
  for (int i = 0; i < 2*nparams; i++) {
    double local[nparams];
    double dx = (i % 2 == 0) ? DX_FIT : - DX_FIT;
    #ifdef WATCH_KELVIN_FITTING
      std::cout << "." << std::flush;
    #endif

    for (int j = 0; j < nparams; j++) {
      local[j] = _PI[j] + (i / 2 == j ? dx : 0.0);
    }

    std::vector<double> apexStresses;
    for (int j = 0; j < count; j++) {
      apexStresses.push_back(1. + local[2]*sin(((double)j) * _sequence->dt[j]));
    }

    c[i].SetLaplace(_sequence->shapes[0].undeformed);
    c[i].SetParameters(1., _sequence->shapes[0].GetPoisson(), local[0], local[1], &_sequence->dt);
    c[i].SetApexStresses(&apexStresses);
    valid_derivatives[i] = c[i].Solve();
  }

  bool valid = true;
  for (int i = 0; i < 2*nparams; i++) {
    valid &= valid_derivatives[i];
  }

  if (!valid) {
    #ifdef WATCH_KELVIN_FITTING
      std::cout << "Derivatives became invalid... Converged." << std::endl;
    #endif
    if (err) *err = err_prev;
    for (int i = 0; i < 2*nparams; i++) {
      if (!valid_derivatives[i]) {
        std::vector<double> apexStresses;
        for (int j = 0; j < count; j++) {
          apexStresses.push_back(1. + INITIAL_STRESS_AMPLITUDE*sin(((double)j) * _sequence->dt[j]));
        }

        c[i].SetLaplace(_sequence->shapes[0].undeformed);
        c[i].SetParameters(1., _sequence->shapes[0].GetPoisson(), INITIAL_K, INITIAL_NU, &_sequence->dt);
        c[i].SetApexStresses(&apexStresses);
        c[i].Solve(true);
      }
    }
    return true;
  }

  #ifdef WATCH_KELVIN_FITTING
    std::cout << std::endl << std::flush;
  #endif

  int offset = 0;
  for (int j = 0; j < count; j++) {
    Kelvin* data = &_sequence->shapes[j];
    PointSet* points = &(*_points)[j];
    /* set up jacobian and residual */
    for (int y = 0; y < points->size(); y++) {
      p = &(*points)[y];
      error = data->ResidualForPoint(*p);
      _F[offset + y] = - error.delta_r;
      _F[offset + y + N] = - error.delta_z;
      for (int x = 0; x < 2*nparams; x+=2) {
        error_x1 = c[x].shapes[j].ResidualForPoint(*p);
        error_x2 = c[x + 1].shapes[j].ResidualForPoint(*p);
        _J(offset + y,x/2) = (error_x1.delta_r - error_x2.delta_r) / (2.*DX_FIT);
        _J(offset + y + N, x/2) = (error_x1.delta_z - error_x2.delta_z) / (2.*DX_FIT);
      }
    }
    offset += points->size();
  }

  arma::solve(_D, _J, _F);
  // Alternatively:
  // _D = arma::pinv(_J) * _F;

  while (arma::norm(_D, 2) > 1.0) {
    _D *= xi;
  }

  int line_limiter = 0;
  while (continue_linesearch) {
    line_limiter++;
    if(line_limiter > 100) {
      continue_linesearch = false;
      converged = true;
    }

    double local[nparams];
    for (int j = 0; j < nparams; j++) {
      local[j] = _PI[j] + _D[j];
    }

    if (local[0] < 0) local[0] = 1e-3;
    if (local[1] < 0) local[1] = 1e-3;

    std::vector<double> apexStresses;
    for (int j = 0; j < count; j++) {
      apexStresses.push_back(1. + local[2]*sin(((double)j) * _sequence->dt[j]));
    }

    _sequence->SetCount(count);
    _sequence->SetLaplace(_sequence->shapes[0].undeformed);
    _sequence->SetParameters(1,
                             _sequence->shapes[0].GetPoisson(),
                             local[0],
                             local[1],
                             &_sequence->dt);

    _sequence->SetApexStresses(&apexStresses);
    _sequence->stressAmplitude = local[2];

    if (_sequence->Solve(true)) { 
      err_next = 0;
      for (int i = 0; i < count; i++) {
        err_next += _sequence->shapes[i].MSError(&(*_points)[i]);
      }
      err_next = sqrt(err_next / (double)count);

      #ifdef WATCH_KELVIN_FITTING
        printf("Err: next: %.8f prev: %.8f\n", err_next, err_prev);
      #endif
      if (arma::norm(_D, 2) < EPS_NEWTON_HOOKE) {
        if (err_next <= err_prev) {
          continue_linesearch = false;
          converged = true;
        } else {
          _sequence->SetCount(count);
          _sequence->SetLaplace(_sequence->shapes[0].undeformed);
          _sequence->SetParameters(1,
                                   _sequence->shapes[0].GetPoisson(),
                                   parameters[0],
                                   parameters[1],
                                   &_sequence->dt);

          _sequence->SetApexStresses(&apexStresses);
          _sequence->stressAmplitude = parameters[2];

          if (_sequence->Solve(true)) {
            converged = true;
            continue_linesearch = false;
          }
        }
        break;
      }
    }

    if (err_next <= err_prev) { 
      break;
    }
    // if (arma::norm(_D, 2) < EPS_NEWTON_HOOKE) {
    //   continue_linesearch = false;
    // }

    #ifdef WATCH_KELVIN_FITTING
      std::cout << std::setw(20) << "shift" << std::endl;
      std::cout << std::setw(20) << _D[0] << std::setw(20) << _D[1] << std::setw(20) << _D[1] << std::endl;
    #endif

    _D *= xi;
  } /* end of linesearch */

  if (!converged) {
    _PI += _D;

    if (_PI[0] < 0) _PI[0] = 1e-4;
    if (_PI[1] < 0) _PI[1] = 1e-4;
    std::vector<double> apexStresses;
    for (int j = 0; j < count; j++) {
      apexStresses.push_back(1. + _PI[2]*sin(((double)j) * _sequence->dt[j]));
    }

    _sequence->SetCount(count);
    _sequence->SetLaplace(_sequence->shapes[0].undeformed);
    _sequence->SetParameters(1,
                             _sequence->shapes[0].GetPoisson(),
                             _PI[0],
                             _PI[1],
                             &_sequence->dt);

    _sequence->SetApexStresses(&apexStresses);
    _sequence->stressAmplitude = _PI[2];

    _sequence->Solve(true);
  }

  if (err) *err = err_next;
  /* clean up */
  p = NULL;
  return converged;
}

bool KelvinFit(KelvinSequence* _sequence, std::vector<PointSet>* _points, int count, double* err) {
  bool converged = false;
  int errorCache = 5;
  double prevErrors[errorCache];
  for (int i = 0; i < errorCache; i++) {
    prevErrors[i] = 1e10;
  }


  int counter = 0;
  while (!converged) {
    converged = KelvinFitIteration(_sequence, _points, count, &prevErrors[counter % errorCache], err);
    ++counter;
    if (counter > 50) {
      #ifdef WATCH_KELVIN_FITTING
        std::cout << "Limiting Iterations... Converged." << std::endl;
      #endif
      return true;
    }
    double dE = 0;
    for (int i = 0; i < errorCache; i++) {
      dE += prevErrors[i];
    }
    if (dE < 1e-5) {
      #ifdef WATCH_KELVIN_FITTING
        std::cout << "Error not changing any more... Converged." << std::endl;
      #endif
      converged = true;
    }
  }

  return converged;
}
