#include <cmath>

static inline bool isInParameterTriangle(double p, double rho) {
    if (4.7* (p - 4.0) > rho)
        return false;
    if (std::pow(p, 2.0010582695119425) * 0.156 < rho)
        return false;
    return true;
}
