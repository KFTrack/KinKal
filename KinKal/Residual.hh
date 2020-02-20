#ifndef KinKal_Residual_hh
#define KinKal_Residual_hh
#include "Math/SVector.h"
#include "Math/SMatrix.h"
namespace KinKal {
  template <unsigned HDIM=1> class Residual {
    public:
      typedef ROOT::Math::SVector<double,HDIM> RVec; // residual vector type
      typedef ROOT::Math::SMatrix<double,HDIM,HDIM,ROOT::Math::MatRepSym<double,HDIM> > RCov;  // associated covariance type
      // accessors
      RVec const& residual() const { return rvec_; }
      RCov const& covariance() const  { return rcov_; }
      // construct from data
      Residual(RVec const& rvec, RCov const& rcov) : rvec_(rvec), rcov_(rcov) {}
    private:
      RVec rvec_; // residual value
      RCov rcov_; // covariance matrix on residual
  };
}
#endif
