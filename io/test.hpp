#ifdef TEST

// Laplace Benchmark time (M1 Pro): 137ms
#define LAPLACE_TARGET_PRECISION 1.414e-17 

// Hooke Benchmark time (M1 Pro): 182ms
#define HOOKE_TARGET_PRECISION 1.447e-09

// Kelvin Benchmark time (M1 Pro): 7778ms
#define KELVIN_TARGET_PRECISION 1.06e-5

#include "interface.hpp"
#include "../hooke/pressureShooting.hpp"
#include "../kelvin/pressureShooting.hpp"
#include "hookeSolver.hpp"
#include <chrono>
#endif
