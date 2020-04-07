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
      float residTime() const { return tresid_; } 
      float resid() const { return resid_; }
      float residVar() const  { return rvar_; }
      PDER const& dRdP() const { return dRdP_; }
      Residual(float tresid, float rvec, float rvar, PDER const& dRdP) : tresid_(tresid), resid_(rvec), rvar_(rvar), dRdP_(dRdP) {}
      Residual() : tresid_(0.0), resid_(0.0), rvar_(-1.0) {}
    private:
      float tresid_; // particle time at which this residual was calculated
      float resid_; // residual value
      float rvar_; // estimated variance of the residual due to sensor measurement uncertainty ONLY
      PDER dRdP_; // derivative of residual WRT the reference parameters
  };

  template <class KTRAJ> std::ostream& operator <<(std::ostream& ost, Residual<KTRAJ> const& res) {
    ost << " residual " << res.resid() << " variance " << res.residVar() << " time " << res.residTime() << " dRdP " << res.dRdP();
    return ost;
  }
}
#endif
