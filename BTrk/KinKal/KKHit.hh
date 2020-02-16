#ifndef KinKal_KKHit_hh
#define KinKal_KKHit_hh
//
//  class to use information from a hit in the Kinematic fit.
//  Used as part of the kinematic Kalman fit
//
#include "BTrk/KinKal/KKWeight.hh"
#include "BTrk/KinKal/TPOCA.hh"
#include "BTrk/KinKal/TLine.hh"

namespace KinKal {

  template <class KTRAJ, unsigned HDIM=1> class KKHit : public KKWeight<KTRAJ> {
    public:

      typedef ROOT::Math::SVector<double,HDIM> RVec; // residual vector
      typedef ROOT::Math::SMatrix<double,HDIM,HDIM,ROOT::Math::MatRepSym<double,HDIM> > RMat;  // associated measurement errors

      virtual unsigned nCons() const override { return NDIM; }
      // access to the TTraj describing this measurement
      virtual TTraj const& traj() const = 0;
      // interpet the TPOCA result as a residual
      virtual void resid(TPOCABase const& tpoca, Resid& resid) const = 0;
      

    private:
      TLine traj_; // trajectory representing this hit
      TPOCA<KTRAJ,TLine> tpoca_; // POCA between the reference trajectory and the line returned by the TrkHit
  };
}
#endif
