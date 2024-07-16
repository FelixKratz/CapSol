#include <vector>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp.h>

class DataSeries {
  private:
  void DestroySpline() {
    if (spline) {
      gsl_spline_free(spline);
      spline = nullptr;
    }
  }

  public:
  std::vector<double> values;
  gsl_spline* spline = nullptr;

  ~DataSeries() { Clear(); }

  void Push(double value) {
    values.push_back(value);

    // Always invalidate the spline when new values are pushed.
    // It should never be expected behaviour to work with an outdated spline...
    DestroySpline();
  }

  void Pop() {
    values.pop_back();
    DestroySpline();
  }

  void Clear() {
    values.clear();
    DestroySpline();
  }

  void Interpolate(DataSeries* x) {
    if (x->values.size() == 0 || values.size() == 0) return;
    DestroySpline();

    spline = gsl_spline_alloc(gsl_interp_cspline, values.size());
    gsl_spline_init(spline, &(x->values[0]), &(values[0]), values.size());
  }
};

