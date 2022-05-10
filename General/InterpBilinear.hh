#ifndef KinKal_InterpBilinear_hh
#define KinKal_InterpBilinear_hh
// class to do a 2D interpolation, which is utilized by the CylBFieldMap interface.
#include<vector>

namespace KinKal {
  class InterpBilinear {
    public:
      using VEC = std::vector<double>;
      using MAT = std::vector<VEC>;
      InterpBilinear(const VEC &x1v, const VEC &x2v, const MAT &ym)
        : m_(x1v.size()), n_(x2v.size()), dx1_(x1v[1]-x1v[0]), dx2_(x2v[1]-x2v[0]), x1(x1v), x2(x2v), y(ym) {}
      double interp(double x1p, double x2p) const;
      // the gradient vector of the interpolating function, e.g. dy/dx_i
      VEC gradient(double x1p, double x2p) const;
    private:
      const int m_, n_;
      const double dx1_, dx2_;
      const VEC &x1, &x2;
      const MAT &y;
  };
}
#endif