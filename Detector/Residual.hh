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
      double measurementVariance() const  { return mvar_; }
      double parameterVariance() const  { return pvar_; }
      double variance() const { return mvar_ + pvar_; }
      DVEC const& dRdP() const { return dRdP_; } // derivative of this residual WRT parameters
      double chisq() const { return (value_*value_)/variance(); }
      double chi() const { return value_/sqrt(variance()); }
      double pull() const { return chi(); }
      Residual(double value, double mvar, double pvar, DVEC const& dRdP) : value_(value), mvar_(mvar), pvar_(pvar), dRdP_(dRdP){}
      Residual() : value_(0.0), mvar_(-1.0), pvar_(-1.0) {}
    private:
      double value_;  // value for this residual
      double mvar_; // estimated variance due to measurement uncertainty
      double pvar_; // estimated variance due to parameter uncertainty
      DVEC dRdP_; // derivative of this residual WRT the trajectory parameters, evaluated at the reference parameters
  };
  std::ostream& operator <<(std::ostream& ost, Residual const& res);
}
#endif
