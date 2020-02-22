#ifndef KinKal_Residual_hh
#define KinKal_Residual_hh
//
//  A Residual describes the tension between an individual sensor measurement and a prediction
//  It does not have fixed units, but is consistent units between value and covariance
//  used as part of the kinematic kalman fit
//
#include "Math/SVector.h"
#include "Math/SMatrix.h"
namespace KinKal {
  class Residual {
    public:
      typedef ROOT::Math::SVector<double,1> RVec; // residual vector type
      typedef ROOT::Math::SMatrix<double,1,1,ROOT::Math::MatRepSym<double,1> > RCov;  // associated covariance type
      // accessors
      RVec const& residual() const { return rvec_; }
      RCov const& covariance() const  { return rcov_; }
      float dRdt() const { return drdt_; } 
      // construct from matrix data
      Residual(RVec const& rvec, RCov const& rcov, double drdt) : rvec_(rvec), rcov_(rcov), drdt_(drdt) {}
      // construct from scalars
      Residual(double r, double ms, double drdt) : rvec_(r), rcov_(ms), drdt_(drdt) {}
      Residual() : drdt_(0.0) {}
    private:
      RVec rvec_; // residual value
      RCov rcov_; // covariance matrix on residual (from measurement errors)
      float drdt_; // how the residual value directly depends on the time component of TDOCA
  };
}
#endif
