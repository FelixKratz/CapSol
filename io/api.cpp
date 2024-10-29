#ifdef API
#include "api.hpp"

typedef std::pair<PointSet, PointSet> ShapePoints;
typedef std::vector<ShapePoints> Shapes;

std::random_device r;
std::mt19937 generator(r());
std::uniform_real_distribution<> p_dist(0.5, 3.95);
std::uniform_real_distribution<> rho_dist(0., 5.);
std::uniform_real_distribution<> K_dist(0.1, 5.0);
std::uniform_real_distribution<> nu_dist(-0.9, 0.9);
std::uniform_real_distribution<> tau_s_0_dist(0.01, 0.9);

class LaplaceAPI {
  public:
  double p_L, rho;
  double error, area, volume;
  bool valid;
  py::array shape;

  LaplaceAPI() {
    do {
      randomShape();
    } while (!valid);
  };

  LaplaceAPI(double p_L, double rho) {
    setParameters(p_L, rho);
  };

  LaplaceAPI(double p_L, double rho, double omega) {
    setParameters(p_L, rho, omega);
  };

  LaplaceAPI(ShapePoints _shape) {
    setShape(&_shape);
  };

  LaplaceAPI(ShapePoints _shape, double p_L_init, double rho_init) {
    SetLaplaceInitialGuess(p_L_init, rho_init);
    setShape(&_shape);
  };


  std::string repr() {
    std::stringstream stream;
    stream << "p_L=" << p_L << ", rho=" << rho << std::endl
           << "A=" << area << ", V=" << volume << std::endl
           << "valid=" << (valid ? "true" : "false") << ", e=" << error << std::endl;
    return stream.str();
  }

  void randomShape() {
    double p = p_dist(generator);
    double r = rho_dist(generator);
    error = 0.;
    area = 0.;
    volume = 0.;
    valid = false;

    while (!isInParameterTriangle(p, r)) {
      p = p_dist(generator);
      r = rho_dist(generator);
    }

    setParameters(p, r);
  }

  void setShape(ShapePoints* _shape) {
    std::vector<PointSet*> pointSets;
    pointSets.push_back(&_shape->first);
    pointSets.push_back(&_shape->second);
    
    error = 1e10;
    laplace_parameters prms;
    fitLaplaceShapes(pointSets, &prms, &error);
    setParameters(prms.p_L, prms.rho);
  };

  bool setParameters(double _p_L, double _rho, double _omega = 2) {
    p_L = _p_L;
    rho = _rho;
    Laplace laplace;
    laplace.SetParameters(p_L, rho, _omega);
    laplace.SetIntegrationStep(INTEGRATION_STEP);
    laplace.Solve();
    
    if (!laplace.IsValid()) {
      valid = false;
      return false;
    }

    volume = laplace.Volume();

    PointSet p;
    generatePointSet(&laplace, &p, false);
    ShapePoints s = { p, p };
    shape = py::cast(s);

    valid = true;
    return true;
  };
};

class HookeAPI {
  public:
  double K, nu, tau_s_0, p_L, rho, p_a;
  double error, volume;
  std::vector<double> tau_s, tau_phi;
  std::vector<double> lambda_s, lambda_phi;

  bool valid;
  py::array shape;
  std::pair<double, double> wrinkling_start;
  std::pair<double, double> wrinkling_stop;

  HookeAPI(double _p_L, double _rho) {
    p_L = _p_L;
    rho = _rho;
    int counter = 0;
    do {
      randomShape();
    } while(!valid && counter++ < 100);
  };

  HookeAPI(double _p_L, double _rho, double _nu) {
    p_L = _p_L;
    rho = _rho;
    nu = _nu;
    int counter = 0;
    do {
      randomShape(false);
    } while(!valid && counter++ < 100);
  };

  HookeAPI(double K, double nu, double tau_s_0, double p_L, double rho) {
    setParameters(K, nu, tau_s_0, p_L, rho);
  };

  HookeAPI(ShapePoints shape, double p_L, double rho) {
    setShape(&shape, p_L, rho);
  };

  HookeAPI(ShapePoints shape, double K_init, double nu_init, double tau_s_0_init, double p_L, double rho) {
    SetHookeInitialGuess(K_init, nu_init, tau_s_0_init);
    setShape(&shape, p_L, rho);
  };

  std::string repr() {
    std::stringstream stream;
    stream << "p_L=" << p_L << ", rho=" << rho << std::endl
           << "K=" << K << ", nu=" << nu << std::endl
           << "tau_s_0=" << tau_s_0 << std::endl
           << "V=" << volume << std::endl
           << "wrinkling=[[" << wrinkling_start.first
                             << ", " << wrinkling_start.second
                             << "],["  << wrinkling_stop.first
                             << ", " << wrinkling_stop.second << "]]" << std::endl
           << "valid=" << (valid ? "true" : "false") << ", e=" << error << std::endl;
    return stream.str();
  }

  void randomShape(bool generate_nu = true) {
    K = K_dist(generator);
    if (generate_nu)
      nu = nu_dist(generator);

    tau_s_0 = tau_s_0_dist(generator);

    valid = false;
    error = 0.;
    volume = 0.;

    setParameters(K, nu, tau_s_0, p_L, rho);
  }

  bool setParameters(double _K, double _nu, double _tau_s_0, double _p_L, double _rho) {
    K = _K;
    nu = _nu;
    tau_s_0 = _tau_s_0;
    p_L = _p_L;
    rho = _rho;
    Laplace laplace;
    laplace.SetParameters(p_L, rho);
    laplace.SetIntegrationStep(INTEGRATION_STEP);
    laplace.Solve();
    if (!laplace.IsValid()) {
      valid = false;
      return false;
    }

    Hooke hooke;
    hooke.SetLaplace(&laplace);
    hooke.SetParameters(1., nu, K);
    hooke.SetIntegrationStep(INTEGRATION_STEP);
    hooke.SetInitialConditions(0.0, 0.0, 0.0, tau_s_0);
    bool hookeShooting = hooke.PressureShooting(true);
    if (!hookeShooting) {
      valid = false;
      return false;
    }
    hooke.Interpolate(true, true);

    PointSet p;
    generatePointSet(&hooke, &p, true);
    ShapePoints s = { p , p };
    shape = py::cast(s);
    valid = true;

    // wrinkling_start = hooke.GetXYPairFromS(hooke.wrinkling_start);
    // wrinkling_stop = hooke.GetXYPairFromS(hooke.wrinkling_stop);
    tau_s = hooke.tau_s.values;
    tau_phi = hooke.tau_phi.values;

    lambda_s = hooke.lambda_s.values;
    lambda_phi = hooke.lambda_phi.values;

    p_a = hooke.GetPressure();
    volume = hooke.Volume();
    return true;
  };

  void setShape(ShapePoints* _shape, double _p_L, double _rho) {
    p_L = _p_L;
    rho = _rho;
    Laplace laplace(p_L, rho);
    laplace.SetIntegrationStep(INTEGRATION_STEP);
    laplace.Solve();
    if (!laplace.IsValid()) {
      valid = false;
      return;
    }

    std::vector<PointSet*> pointSets;
    pointSets.push_back(&_shape->first);
    pointSets.push_back(&_shape->second);

    error = 1e10;
    hooke_parameters prms;
    fitHookeShapes(&laplace, pointSets, &prms, &error);
    setParameters(prms.K, prms.Psi, prms.tau_s_0, _p_L, _rho);
  };
};


class KelvinAPI {
  public:
    py::array shapes;
    py::array apexStresses;
    double K, nu, eta, p_L, rho;
    double error;
    bool valid;
    int count;

    KelvinAPI(std::vector<double> stresses, double _K, double _nu, double _eta, double _p_L, double _rho, int periods) {
      setParameters(&stresses, _K, _nu, _eta, _p_L, _rho, periods);
    }

    KelvinAPI(double amplitude, int periods, int count, double _K, double _nu, double _eta, double _p_L, double _rho) {
      std::vector<double> stresses;
      double dt = ((double)periods) * 2.* M_PI / ((double)count - 1.);
      for (int i = 0; i < count; i++) {
        stresses.push_back(1.0 + amplitude*sin(((double)i) * dt));
      }

      setParameters(&stresses, _K, _nu, _eta, _p_L, _rho, periods);
    }

    KelvinAPI(Shapes shapes, double _p_L, double _rho, int periods) {
      setShapes(&shapes, _p_L, _rho, periods);
    }

    std::string repr() {
      std::stringstream stream;
      stream << "p_L=" << p_L << ", rho=" << rho << std::endl
             << "K=" << K << ", nu=" << nu << ", eta=" << eta << std::endl
             << "valid=" << (valid ? "true" : "false") << ", e=" << error
             << std::endl;
      return stream.str();
    }

    bool setParameters(std::vector<double>* stresses, double _K, double _nu, double _eta, double _p_L, double _rho, int periods) {
      K = _K;
      nu = _nu;
      eta = _eta;
      p_L = _p_L;
      rho = _rho;

      apexStresses = py::cast(stresses);

      Laplace laplace(p_L, rho);
      laplace.SetIntegrationStep(INTEGRATION_STEP);
      laplace.Solve();
      if (!laplace.IsValid()) {
        valid = false;
        return false;
      }

      KelvinSequence sequence(stresses->size());
      sequence.SetLaplace(&laplace);
      double dt = periods * 2.*M_PI / ((double)stresses->size() - 1.);
      sequence.SetParameters(1., nu, K, eta, dt);
      sequence.SetApexStresses(stresses);

      if (!sequence.Solve(true)) {
        valid = false;
        return false;
      }

      sequence.Interpolate();
      Shapes _shapes;
      for (int i = 0; i < sequence.count; i++) {
        PointSet p;
        generatePointSet(&sequence.shapes[i], &p, true);
        _shapes.push_back({p, p});
      }

      count = _shapes.size();
      shapes = py::cast(_shapes);
      valid = true;

      return true;
    }

    void setShapes(Shapes* shapes, double _p_L, double _rho, int periods) {

      p_L = _p_L;
      rho = _rho;

      Laplace laplace(p_L, rho);
      laplace.SetIntegrationStep(INTEGRATION_STEP);
      laplace.Solve();

      if (!laplace.IsValid()) {
        printf("[x] Laplace Invalid\n");
        valid = false;
        return;
      } 

      std::vector<PointSet> leftPointSet;
      std::vector<PointSet> rightPointSet;
      for (int i = 0; i < shapes->size(); i++) {
        leftPointSet.push_back((*shapes)[i].first);
        rightPointSet.push_back((*shapes)[i].second);
      }

      struct kelvin_parameters fit_parameters_l;
      double err_l = 1e10;
      bool success = fitKelvinSequence(&leftPointSet,
                                       &laplace,
                                       &fit_parameters_l,
                                       leftPointSet.size(),
                                       periods,
                                       &err_l              );

      std::vector<double> apexStresses;
      double dt = ((double)periods) * 2.* M_PI / ((double)shapes->size() - 1.);
      for (int i = 0; i < shapes->size(); i++) {
        apexStresses.push_back(1.0 + fit_parameters_l.dtau_s_0*sin(((double)i) * dt));
      }

      setParameters(&apexStresses, fit_parameters_l.K,
                                   fit_parameters_l.Psi,
                                   fit_parameters_l.eta,
                                   fit_parameters_l.p_L,
                                   fit_parameters_l.rho,
                                   periods);

      valid = true;
      error = err_l;
    }
};

class ContactAPI {
  public:
  Laplace laplace_u, laplace_d;
  SlipContact contact;
  std::vector<double> r_u, z_u;
  std::vector<double> r_d, z_d;
  std::vector<double> tau_s_u, tau_s_d, tau_phi_u, tau_phi_d;
  std::vector<double> lambda_s_u, lambda_s_d, lambda_phi_u, lambda_phi_d;
  double f;
  double V_u;
  double V_d;
  double p_a_u;
  double p_a_d;

  bool valid;

  // PressureShooting
  ContactAPI(double p_L_u, double rho_u, double p_L_d, double rho_d,
             double nu_u, double K_u, double tau_s_0_u,
             double nu_d, double K_d, double tau_s_0_d,
             double contact_length, double tension_ratio,
             double tension_inner_u, double tension_inner_d,
             double tension_ud) {

    LaplaceParameters p_laplace_u(p_L_u, rho_u);
    // The laplace paramters of the lower shape must be scaled such that they
    // are in the proper dimensionless units for integration (i.e. we get
    // parameters with \gamma^u as a tension scale and the solver needs
    // parameters in terms of \gamma^d as a tension scale). Dividing by the
    // tension ratio is sufficient.
    LaplaceParameters p_laplace_d(p_L_d / tension_ratio, rho_d / tension_ratio);

    double g_u = K_u * (1. - nu_u) / (1. + nu_u);
    double g_d = K_d * (1. - nu_d) / (1. + nu_d);
    HookeParameters p_hooke_u(g_u, K_u, tau_s_0_u, 0);
    HookeParameters p_hooke_d(g_d, K_d, tau_s_0_d, 0);

    laplace_u.SetParameters(p_laplace_u);
    laplace_u.SetIntegrationStep(INTEGRATION_STEP);
    laplace_u.Solve();
    if (!laplace_u.IsValid()) {
      printf("Error: Invalid Laplace Shape...\n");
    }

    laplace_d.SetParameters(p_laplace_d);
    laplace_d.SetIntegrationStep(INTEGRATION_STEP);
    laplace_d.Solve();
    if (!laplace_d.IsValid()) {
      printf("Error: Invalid Laplace Shape...\n");
    }

    contact.SetLaplace(&laplace_u, &laplace_d);
    contact.SetIntegrationStep(INTEGRATION_STEP);
    contact.SetParameters(SlipContactParameters(p_hooke_u, p_hooke_d, contact_length, tension_ratio, tension_inner_u, tension_inner_d, tension_ud));
    contact.SetInitialConditions(0, 0, 0, 0, p_hooke_u.tau_s_0, p_hooke_d.tau_s_0);

    if (K_u == K_d
        && nu_u == nu_d
        && tau_s_0_u == tau_s_0_d
        && p_L_u == p_L_d
        && rho_u == rho_d
        && tension_inner_u == tension_inner_d) {
      contact.MarkSymmetric();
    }
    valid = contact.PressureShooting(true);

    r_u = contact.hooke_u.r.values;
    z_u = contact.hooke_u.z.values;

    r_d = contact.hooke_d.r.values;
    z_d = contact.hooke_d.z.values;

    tau_s_u = contact.hooke_u.tau_s.values;
    tau_s_d = contact.hooke_d.tau_s.values;

    tau_phi_u = contact.hooke_u.tau_phi.values;
    tau_phi_d = contact.hooke_d.tau_phi.values;

    lambda_s_u = contact.hooke_u.lambda_s.values;
    lambda_s_d = contact.hooke_d.lambda_s.values;

    lambda_phi_u = contact.hooke_u.lambda_phi.values;
    lambda_phi_d = contact.hooke_d.lambda_phi.values;
    f = contact.parameters.force;

    V_u = contact.hooke_u.Volume();
    V_d = contact.hooke_d.Volume();
    p_a_u = contact.hooke_u.GetPressure();
    p_a_d = contact.hooke_d.GetPressure();
  };

  // VolumeShooting
  ContactAPI(double p_L_u, double rho_u, double p_L_d, double rho_d,
             double nu_u, double K_u,
             double nu_d, double K_d,
             double contact_length, double tension_ratio,
             double tension_inner_u, double tension_inner_d,
             double tension_ud) {

    LaplaceParameters p_laplace_u(p_L_u, rho_u);
    // The laplace paramters of the lower shape must be scaled such that they
    // are in the proper dimensionless units for integration (i.e. we get
    // parameters with \gamma^u as a tension scale and the solver needs
    // parameters in terms of \gamma^d as a tension scale). Dividing by the
    // tension ratio is sufficient.
    LaplaceParameters p_laplace_d(p_L_d / tension_ratio, rho_d / tension_ratio);

    double g_u = K_u * (1. - nu_u) / (1. + nu_u);
    double g_d = K_d * (1. - nu_d) / (1. + nu_d);
    HookeParameters p_hooke_u(g_u, K_u, 0.5, 0);
    HookeParameters p_hooke_d(g_d, K_d, 0.5, 0);

    laplace_u.SetParameters(p_laplace_u);
    laplace_u.SetIntegrationStep(INTEGRATION_STEP);
    laplace_u.Solve();
    if (!laplace_u.IsValid()) {
      printf("Error: Invalid Laplace Shape...\n");
      valid = false;
      return;
    }

    laplace_d.SetParameters(p_laplace_d);
    laplace_d.SetIntegrationStep(INTEGRATION_STEP);
    laplace_d.Solve();
    if (!laplace_d.IsValid()) {
      printf("Error: Invalid Laplace Shape...\n");
      valid = false;
      return;
    }

    contact.SetLaplace(&laplace_u, &laplace_d);
    contact.SetIntegrationStep(INTEGRATION_STEP);
    contact.SetParameters(SlipContactParameters(p_hooke_u, p_hooke_d, contact_length, tension_ratio, tension_inner_u, tension_inner_d, tension_ud));
    contact.SetInitialConditions(0, 0, 0, 0, p_hooke_u.tau_s_0, p_hooke_d.tau_s_0);

    if (K_u == K_d
        && nu_u == nu_d
        && p_L_u == p_L_d
        && rho_u == rho_d
        && tension_inner_u == tension_inner_d) {
      contact.MarkSymmetric();
    }
    valid = contact.VolumeShooting(true);

    r_u = contact.hooke_u.r.values;
    z_u = contact.hooke_u.z.values;

    r_d = contact.hooke_d.r.values;
    z_d = contact.hooke_d.z.values;

    tau_s_u = contact.hooke_u.tau_s.values;
    tau_s_d = contact.hooke_d.tau_s.values;

    tau_phi_u = contact.hooke_u.tau_phi.values;
    tau_phi_d = contact.hooke_d.tau_phi.values;

    lambda_s_u = contact.hooke_u.lambda_s.values;
    lambda_s_d = contact.hooke_d.lambda_s.values;

    lambda_phi_u = contact.hooke_u.lambda_phi.values;
    lambda_phi_d = contact.hooke_d.lambda_phi.values;
    f = contact.parameters.force;

    V_u = contact.hooke_u.Volume();
    V_d = contact.hooke_d.Volume();
    p_a_u = contact.hooke_u.GetPressure();
    p_a_d = contact.hooke_d.GetPressure();
  };
};

class ContactUnitCellAPI {
  public:
  Laplace laplace_u, laplace_d;
  SlipContact contact;
  std::vector<double> r, z;
  std::vector<double> tau_s, tau_phi;
  std::vector<double> lambda_s, lambda_phi;
  double f;
  double V;
  double p_a;

  bool valid;

  // PressureShooting
  ContactUnitCellAPI(double nu, double K, double tau_s_0,
                     double contact_length, double tension_ratio,
                     double tension_inner, double tension_ud) {

    LaplaceParameters p_laplace_u(3.99999, 0.0);
    double g = K * (1. - nu) / (1. + nu);
    HookeParameters p_hooke_u(g, K, tau_s_0, 0);

    laplace_u.SetParameters(p_laplace_u);
    laplace_u.SetIntegrationStep(1e-4);
    laplace_u.Solve();
    if (!laplace_u.IsValid()) {
      printf("Error: Invalid Laplace Shape...\n");
      valid = false;
      return;
    }

    contact.SetLaplace(&laplace_u, &laplace_u);
    contact.SetIntegrationStep(1e-4);
    contact.SetParameters(SlipContactParameters(p_hooke_u, p_hooke_u, contact_length, tension_ratio, tension_inner, tension_inner, tension_ud));
    contact.SetInitialConditions(0, 0, 0, 0, p_hooke_u.tau_s_0, p_hooke_u.tau_s_0);

    contact.MarkSymmetric();
    valid = contact.AngleShooting(true);

    r = contact.hooke_u.r.values;
    z = contact.hooke_u.z.values;

    tau_s = contact.hooke_u.tau_s.values;
    tau_phi = contact.hooke_u.tau_phi.values;

    lambda_s = contact.hooke_u.lambda_s.values;
    lambda_phi = contact.hooke_u.lambda_phi.values;

    f = contact.parameters.force;
    V = contact.hooke_u.Volume();
    p_a = contact.hooke_u.GetPressure();
  };
};

class BendingAPI {
  public:
  std::vector<double> r, z;

  BendingAPI(double p_L, double L, double K, double E_B, double p_a, double tau_s_0, double m_s_0, double q_plus) {
    Bending bending;
    bending.SetIntegrationStep(INTEGRATION_STEP);

    LaplaceParameters p_laplace(p_L, 0.0);
    Laplace laplace;
    laplace.SetParameters(p_laplace);
    laplace.SetIntegrationStep(INTEGRATION_STEP);
    laplace.Solve();
    if (!laplace.IsValid()) {
      printf("Error: Invalid Laplace Shape...\n");
    }

    BendingParameters p_bending(L, 0, K, E_B, 0.5, p_a, tau_s_0, q_plus);
    bending.SetLaplace(&laplace);
    bending.SetParameters(p_bending);
    bending.SetInitialConditions(0, 0, 0, m_s_0, tau_s_0);
    bending.Integrate();

    r = bending.r.values;
    z = bending.z.values;
  }
};


PYBIND11_MODULE(CapSol, m) {
    m.doc() = "Library for fitting drops and viscoelastic capsules";
    py::class_<BendingAPI>(m, "Bending")
      .def(py::init<double, double, double, double, double, double, double, double>())
      .def_readwrite("r", &BendingAPI::r)
      .def_readwrite("z", &BendingAPI::z);

    py::class_<ContactUnitCellAPI>(m, "ContactUnitCell")
      .def(py::init<double, double, double,
                    double, double, double, double>())
      .def_readwrite("r", &ContactUnitCellAPI::r)
      .def_readwrite("z", &ContactUnitCellAPI::z)
      .def_readwrite("tau_s", &ContactUnitCellAPI::tau_s)
      .def_readwrite("tau_phi", &ContactUnitCellAPI::tau_phi)
      .def_readwrite("lambda_s", &ContactUnitCellAPI::lambda_s)
      .def_readwrite("lambda_phi", &ContactUnitCellAPI::lambda_phi)
      .def_readwrite("f", &ContactUnitCellAPI::f)
      .def_readwrite("V", &ContactUnitCellAPI::V)
      .def_readwrite("p_a", &ContactUnitCellAPI::p_a)
      .def_readwrite("valid", &ContactUnitCellAPI::valid);

    py::class_<ContactAPI>(m, "Contact")
      .def(py::init<double, double, double, double, double,
                    double, double, double, double, double,
                    double, double, double, double, double>())
      .def(py::init<double, double, double, double, double,
                    double, double, double, double, double,
                    double, double, double>())
      .def_readwrite("r_u", &ContactAPI::r_u)
      .def_readwrite("z_u", &ContactAPI::z_u)
      .def_readwrite("r_d", &ContactAPI::r_d)
      .def_readwrite("z_d", &ContactAPI::z_d)
      .def_readwrite("tau_s_u", &ContactAPI::tau_s_u)
      .def_readwrite("tau_s_d", &ContactAPI::tau_s_d)
      .def_readwrite("tau_phi_u", &ContactAPI::tau_phi_u)
      .def_readwrite("tau_phi_d", &ContactAPI::tau_phi_d)
      .def_readwrite("lambda_s_u", &ContactAPI::lambda_s_u)
      .def_readwrite("lambda_s_d", &ContactAPI::lambda_s_d)
      .def_readwrite("lambda_phi_u", &ContactAPI::lambda_phi_u)
      .def_readwrite("lambda_phi_d", &ContactAPI::lambda_phi_d)
      .def_readwrite("f", &ContactAPI::f)
      .def_readwrite("V_u", &ContactAPI::V_u)
      .def_readwrite("V_d", &ContactAPI::V_d)
      .def_readwrite("p_a_u", &ContactAPI::p_a_u)
      .def_readwrite("p_a_d", &ContactAPI::p_a_d)
      .def_readwrite("valid", &ContactAPI::valid);

    py::class_<LaplaceAPI>(m, "Laplace")
      .def(py::init<>())
      .def(py::init<double, double>())
      .def(py::init<double, double, double>())
      .def(py::init<ShapePoints>())
      .def(py::init<ShapePoints, double, double>())
      .def("__repr__", &::LaplaceAPI::repr)
      .def_readwrite("p_L", &LaplaceAPI::p_L)
      .def_readwrite("rho", &LaplaceAPI::rho)
      .def_readwrite("shape", &LaplaceAPI::shape)
      .def_readwrite("area", &LaplaceAPI::area)
      .def_readwrite("volume", &LaplaceAPI::volume)
      .def_readwrite("valid", &LaplaceAPI::valid)
      .def_readwrite("error", &LaplaceAPI::error)
      .def_readwrite("area", &LaplaceAPI::area)
      .def_readwrite("volume", &LaplaceAPI::volume);

    py::class_<HookeAPI>(m, "Hooke")
      .def(py::init<double, double, double, double, double>())
      .def(py::init<ShapePoints, double, double>())
      .def(py::init<ShapePoints, double, double, double, double, double>())
      .def(py::init<double, double>())
      .def(py::init<double, double, double>())
      .def("__repr__", &::HookeAPI::repr)
      .def_readwrite("K", &HookeAPI::K)
      .def_readwrite("nu", &HookeAPI::nu)
      .def_readwrite("tau_s_0", &HookeAPI::tau_s_0)
      .def_readwrite("p_a", &HookeAPI::p_a)
      .def_readwrite("shape", &HookeAPI::shape)
      .def_readwrite("valid", &HookeAPI::valid)
      .def_readwrite("error", &HookeAPI::error)
      .def_readwrite("tau_s", &HookeAPI::tau_s)
      .def_readwrite("tau_phi", &HookeAPI::tau_phi)
      .def_readwrite("lambda_s", &HookeAPI::lambda_s)
      .def_readwrite("lambda_phi", &HookeAPI::lambda_phi)
      .def_readwrite("volume", &HookeAPI::volume)
      .def_readwrite("wrinkling_start", &HookeAPI::wrinkling_start)
      .def_readwrite("wrinkling_stop", &HookeAPI::wrinkling_stop);

    py::class_<KelvinAPI>(m, "Kelvin")
      .def(py::init<std::vector<double>, double, double, double, double, double, int>())
      .def(py::init<double, int, int, double, double, double, double, double>())
      .def(py::init<Shapes, double, double, int>())
      .def("__repr__", &::KelvinAPI::repr)
      .def_readwrite("K", &KelvinAPI::K)
      .def_readwrite("eta", &KelvinAPI::eta)
      .def_readwrite("nu", &KelvinAPI::nu)
      .def_readwrite("tau_s_0", &KelvinAPI::apexStresses)
      .def_readwrite("shapes", &KelvinAPI::shapes)
      .def_readwrite("valid", &KelvinAPI::valid)
      .def_readwrite("error", &KelvinAPI::error)
      .def_readwrite("count", &KelvinAPI::count);
}

#endif
