#pragma once

#define EPS_NEWTON_LAPLACE 1e-8
#define EPS_NEWTON_HOOKE 1e-8
#define INTEGRATION_STEP 1e-3
#define GAMMA_SCALE 1.0
#define SHOOTING_INTERVALS 40

// This should not be changed to be much lower, since it will lead to severe
// numerical floating point precision problems. Especially if the shape is
// more insensitive to the parameter changes. The shape changes might be
// several orders of magnitude smaller. When the derivatives in the parameter
// landscape become numerically zero *before* we multiply with 1./DX_FIT we
// have a problem.
#define DX_FIT 1e-5
