#include "KinKal/General/InterpBilinear.hh"

#include <algorithm>

namespace KinKal {
  using VEC = std::vector<double>;
  using MAT = std::vector<std::vector<double>>;

  double InterpBilinear::interp(double x1p, double x2p) const {
    // interpolate values using bilinear interpolation formula at a point x1p, x2p
    int i, j;
    double yy, t, u;
    // find i & j of point at lower left corner of square that encloses x1p,x2p
    i = (int)((x1p-x1[0])/dx1_);
    j = (int)((x2p-x2[0])/dx2_);
    // make sure in available index range
    i = std::max(0, std::min(m_-2, i));
    j = std::max(0, std::min(n_-2, j));
    // calculate fractional location of x1p,x2p in square of points
    t = (x1p-x1[i])/dx1_;
    u = (x2p-x2[j])/dx2_;
    // bilinear interpolation 
    yy = (1.-t)*(1.-u)*y[i][j] + t*(1.-u)*y[i+1][j] + (1.-t)*u*y[i][j+1] + t*u*y[i+1][j+1];
    return yy;
  }

  VEC InterpBilinear::gradient(double x1p, double x2p) const {
    // calculate derivative of interpolation function at a point x1p, x2p
    int i, j;
    double t, u, dtx1, dux2;
    VEC dyxi(2);
    // find i & j of point at lower left corner of square that encloses x1p,x2p
    i = (int)((x1p-x1[0])/dx1_);
    j = (int)((x2p-x2[0])/dx2_);
    // make sure in available index range
    i = std::max(0, std::min(m_-2, i));
    j = std::max(0, std::min(n_-2, j));
    // calculate fractional location of x1p,x2p in square of points
    t = (x1p-x1[i])/dx1_;
    u = (x2p-x2[j])/dx2_;
    // calculate numerical derivatives
    dtx1 = 1./(x1[i+1]-x1[i]);
    dux2 = 1./(x2[j+1]-x2[j]);
    dyxi[0] = dtx1 * ((1.-u)*(y[i+1][j]-y[i][j]) + u*(y[i+1][j+1]-y[i][j+1]));
    dyxi[1] = dux2 * ((1.-t)*(y[i][j+1]-y[i][j]) + t*(y[i+1][j+1]-y[i+1][j]));
    return dyxi; 
  }
}