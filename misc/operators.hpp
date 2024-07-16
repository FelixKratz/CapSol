#pragma once
// These are the operators all Variable template arguments should inherit from
template<class Variables>
struct Operators {
  Variables operator+(const Variables& v) {
    Variables ret;
    for (int i = 0; i < sizeof(Variables) / sizeof(double); i++) {
      *((double*)&ret + i) = *((double*)this + i) + *((double*)&v + i);
    }
    return ret;
  }

  void operator+=(const Variables& v) {
    for (int i = 0; i < sizeof(Variables) / sizeof(double); i++) {
      *((double*)this + i) += *((double*)&v + i);
    }
  }

  Variables operator*(const double scale) {
    Variables ret;
    for (int i = 0; i < sizeof(Variables) / sizeof(double); i++) {
      *((double*)&ret + i) = *((double*)this + i) * scale;
    }
    return ret;
  }

  double& operator[](const int index) {
    return *((double*)this + index);
  }
};
