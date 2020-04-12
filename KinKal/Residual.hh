#ifndef KinKal_Residual_hh
#define KinKal_Residual_hh
//
//  A Residual describes the 1-dimensional tension between an individual sensor measurement and a prediction of a particle trajectory
//  It does not have fixed units, but must have consistent units between value, variance and derivatives.
//  It is based implicitly on the concept of TPOCA (time point of closest approach) between a measurement and a prediction
//  The residual value may depend on space, time, or both
//  used as part of the kinematic kalman fit
//
#include <ostream>
namespace KinKal {
  template <class KTRAJ> class Residual {
    public:
      typedef typename KTRAJ::PDER PDER; // forward derivative type from the particle trajectory
      // accessors
      float time() const { return time_; } 
      float value() const { return value_; }
      float variance() const  { return var_; }
      PDER const& dRdP() const { return dRdP_; }
      Residual(float time, float value, float var, PDER const& dRdP) : time_(time), value_(value), var_(var), dRdP_(dRdP) {}
      Residual() : time_(0.0), value_(0.0), var_(-1.0) {}
    private:
      float time_; // particle time associated with this residual
      float value_; // residual value
      float var_; // estimated variance of the residual due to sensor measurement uncertainty ONLY
      PDER dRdP_; // derivative of residual WRT the reference parameters
  };

  template <class KTRAJ> std::ostream& operator <<(std::ostream& ost, Residual<KTRAJ> const& res) {
    ost << " residual " << res.value() << " variance " << res.variance() << " time " << res.time() << " dRdP " << res.dRdP();
    return ost;
  }
}
#endif
