#ifndef KinKal_KKWireHit_hh
#define KinKal_KKWireHit_hh
//
//  class to use information from a hit in the Kinematic fit.
//  Used as part of the kinematic Kalman fit
//
#include "KinKal/KKWeight.hh"
#include "KinKal/TPOCA.hh"
#include "KinKal/TLine.hh"

namespace KinKal {

  template class WireHit : public TrkHit<1> {
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
      TPOCA<KTRAJ,TLINE> tpoca_; // POCA between the reference
  };
}
#endif
