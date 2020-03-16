#ifndef KinKal_Residual_hh
#define KinKal_Residual_hh
//
//  A Residual describes the tension between an individual sensor measurement and a prediction
//  It does not have fixed units, but must have consistent units between value and variance
//  The relationship between the residual value and time and space are explicit
//  used as part of the kinematic kalman fit
//
#include <ostream>
namespace KinKal {
  class Residual {
    public:
      // accessors
      float resid() const { return resid_; }
      float residVar() const  { return rvar_; }
      float dRdD() const { return drdd_; } 
      float dRdT() const { return drdt_; } 
      Residual(float rvec, float rvar, float drdd, float drdt) : resid_(rvec), rvar_(rvar), drdd_(drdd), drdt_(drdt) {}
      Residual() : drdd_(0.0), drdt_(0.0) {}
    private:
      float resid_; // residual value
      float rvar_; // estimated variance of the residual
      float drdd_; // how the residual value directly depends on distance
      float drdt_; // how the residual value directly depends on time
  };
  std::ostream& operator <<(std::ostream& ost, Residual const& res);
}
#endif
