#ifndef KinKal_Residual_hh
#define KinKal_Residual_hh
//
//  A Residual describes the 1-dimensional tension between an individual sensor measurement and a prediction of a particle trajectory
//  It does not have fixed units, but must have consistent units between value, variance and derivatives.
//  It is based on the concept of TPOCA (time point of closest approach) between a measurement and a prediction
//  The residual value may depend on any aspect of the measurement, reduced to a single dimension.
//  used as part of the kinematic kalman fit
//
#include "KinKal/TPocaBase.hh"
#include <ostream>
#include <string>

namespace KinKal {
  template <size_t DDIM> class Residual {
    public:
      enum rdim {unknown=-1,dtime=0,distance}; // residual dimension
      static std::string dimensionName(rdim dim) {
	switch(dim) {
	  case unknown: default:
	    return std::string("Unknown");
	    break;
	  case dtime:
	    return std::string("Time");
	    break;
	  case distance:
	    return std::string("Distance");
	    break;
	}
      }
      typedef ROOT::Math::SVector<double,DDIM> DVEC; // data vector
      // accessors
      rdim dimension() const { return dim_; }
      TPocaBase const& tPoca() const { return tpoca_; }
      double time() const { return tpoca_.particleToca(); }
      double value() const { return value_; }
      double variance() const  { return var_; }
      DVEC const& dRdP() const { return dRdP_; }
      Residual(rdim dim, TPocaBase const& tpoca, double value, double var, DVEC const& dRdP) : dim_(dim), tpoca_(tpoca), value_(value), var_(var), dRdP_(dRdP) {}
      Residual() : dim_(unknown), value_(0.0), var_(-1.0) {}
    private:
      rdim dim_; // dimension of this residual
      TPocaBase tpoca_; // TPOCA associated with this residual
      double value_;  // value for this residual
      double var_; // estimated variance of the residual due to sensor measurement uncertainty ONLY
      DVEC dRdP_; // derivative of residual WRT the reference parameters
  };

  template <size_t DDIM> std::ostream& operator <<(std::ostream& ost, Residual<DDIM> const& res) {
    ost << " residual dimension " << Residual<DDIM>::dimensionName(res.dimension()) << " value " << res.value() << " variance " << res.variance() << " time " << res.time() << " dRdP " << res.dRdP();
    return ost;
  }
}
#endif
