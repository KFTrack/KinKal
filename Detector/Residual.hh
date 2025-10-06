#ifndef KinKal_Residual_hh
#define KinKal_Residual_hh
//
//  A Residual describes the 1-dimensional tension between an individual sensor measurement and a prediction of a particle trajectory.
//  It does not have fixed units, but must have consistent units between value, variance and derivatives.
//  Its value may depend on any aspect of the measurement, reduced to a single dimension.
//
#include <ostream>
#include "KinKal/General/Vectors.hh"
#include "KinKal/General/Weights.hh"

namespace KinKal {
  class Residual {
    public:
      // accessors
      double value() const { return value_; } // residual value
      double measurementVariance() const  { return mvar_; }
      double parameterVariance() const  { return pvar_; }
      double variance() const { return mvar_ + pvar_; }
      DVEC const& dRdP() const { return dRdP_; } // derivative of this residual WRT parameters
      bool active() const { return active_; }
      double chisq() const { return active() ? (value_*value_)/variance() : 0.0; }
      double chi() const { return active() ? value_/sqrt(variance()): 0.0; }
      double pull() const { return chi(); }
      unsigned nDOF() const { return active() ? 1 : 0; }
      // calculate the weight WRT some parameters implied by this residual.  Optionally scale the variance
      Weights weight(DVEC const& params, double varscale=1.0) const;
      Residual(double value, double mvar, double pvar, DVEC const& dRdP,bool active=true) :
        value_(value), mvar_(mvar), pvar_(pvar), dRdP_(dRdP), active_(active){}
      Residual() : value_(0.0), mvar_(-1.0), pvar_(-1.0), active_(false) {}
    private:
      double value_;  // value for this residual
      double mvar_; // estimated variance due to measurement uncertainty
      double pvar_; // estimated variance due to parameter uncertainty
      DVEC dRdP_; // derivative of this residual WRT the trajectory parameters, evaluated at the reference parameters
      bool active_;
  };
  std::ostream& operator <<(std::ostream& ost, Residual const& res);
}
#endif
