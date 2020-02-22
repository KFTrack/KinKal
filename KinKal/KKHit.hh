#ifndef KinKal_KKHit_hh
#define KinKal_KKHit_hh
//
//  class to use information from a hit in the Kinematic fit.
//  Used as part of the kinematic Kalman fit
//
#include "KinKal/KKWeight.hh"
#include "KinKal/THit.hh"
#include "KinKal/TPOCA.hh"
#include "KinKal/TLine.hh"
#include "KinKal/Residual.hh"

namespace KinKal {
  class TTraj;
  template <class KTRAJ> class KKHit : public KKWeight<KTRAJ> {
    public:
      virtual unsigned nDOF() const override { return thit_.nDOF(); }
      THit const& hit() const { return thit_; }
      virtual void update(KTRAJ const& ref)  override { 
	KKEffect<KTRAJ>::reftraj_ = &ref;
//	thit_.update(ref); FIXME!
      }
      virtual double time() const override { return tdpoca_.poca0().T(); } // time on the main trajectory
      virtual ~KKHit(){}
      // construct from a hit and reference trajectory
      KKHit(THit const& thit, KTRAJ const& reftraj);
    private:
      THit const& thit_; // hit used for this constraint
      TDPOCA<KTRAJ,TLine> tdpoca_; // POCA between the reference trajectory and the sensor trajectory
      Residual rresid_; // residual between the hit and reference
  };

  template<class KTRAJ> KKHit<KTRAJ>::KKHit(THit const& thit, KTRAJ const& reftraj) : KKWeight<KTRAJ>(reftraj), thit_(thit) , tdpoca_(reftraj,thit_.sensorTraj()) {
    // translate tdpoca into a residual
    thit_.resid(tdpoca_,rresid_);
    // compute the total derivative (space and time combined) of the residual WRT ktraj parameters
    auto dRdP = rresid_.dRdD()*tdpoca_.dDdP() + rresid_.dRdT()*tdpoca_.dTdP();
    // convert this to a Nx1 matrix (for root)
    ROOT::Math::SMatrix<double,KTRAJ::NParams(),1> dRdPM;
    dRdPM.Place_in_row(dRdP,0,0);
    // convert the variance into a 1X1 matrix
    ROOT::Math::SMatrix<double, 1,1, ROOT::Math::MatRepSym<double,1> > RVarM;
    RVarM(0,0) = rresid_.residVar();
    // expand these into the weight matrix
    KKWeight<KTRAJ>::weight_.weightMat() = ROOT::Math::Similarity(dRdPM,RVarM);
    // reference weight vector from reference parameters
    auto refvec = KKWeight<KTRAJ>::weight().weightMat()*reftraj.params().parameters();
    // translate residual value into weight vector WRT the reference parameters
    auto delta = dRdP*rresid_.resid()*rresid_.residVar();
    // add change WRT reference; sign convention reflects resid = measure - prediction
    KKWeight<KTRAJ>::weight_.weightVec() = refvec + delta; 
  }
}
#endif
