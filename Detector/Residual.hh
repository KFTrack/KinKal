#ifndef KinKal_Residual_hh
#define KinKal_Residual_hh
//
//  A Residual describes the 1-dimensional tension between an individual sensor measurement and a prediction of a particle trajectory.
//  It does not have fixed units, but must have consistent units between value, variance and derivatives.
//  Its value may depend on any aspect of the measurement, reduced to a single dimension.
//
#include <ostream>
#include "KinKal/General/Vectors.hh"

namespace KinKal {
  class Residual {
    public:
      // accessors
      double value() const { return value_; } // residual value
      double variance() const  { return var_; } // variance on value, based on measurement uncertainty
      DVEC const& dRdP() const { return dRdP_; } // derivative of this residual WRT parameters
      Residual(double value, double var, DVEC const& dRdP) : value_(value), var_(var), dRdP_(dRdP){}
      Residual() : value_(0.0), var_(-1.0) {}
    private:
      double value_;  // value for this residual
      double var_; // estimated variance of the residual due to sensor measurement uncertainty ONLY
      DVEC dRdP_; // derivative of this residual WRT the trajectory parameters, evaluated at the reference parameters
  };
  std::ostream& operator <<(std::ostream& ost, Residual const& res);
}
#endif
