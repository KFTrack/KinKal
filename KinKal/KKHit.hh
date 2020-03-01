#ifndef KinKal_KKHit_hh
#define KinKal_KKHit_hh
//
//  class to use information from a hit in the Kinematic fit.
//  Used as part of the kinematic Kalman fit
//
#include "KinKal/KKWeight.hh"
#include "KinKal/PKTraj.hh"
#include "KinKal/THit.hh"
#include "KinKal/TPoca.hh"
#include "KinKal/TLine.hh"
#include "KinKal/Residual.hh"

namespace KinKal {
  class TTraj;
  template <class KTRAJ> class KKHit : public KKWeight<KTRAJ> {
    public:
      typedef KKEff<KTRAJ> KKEFF;
      typedef PKTraj<KTRAJ> PKTRAJ;
      typedef TDPoca<PKTRAJ,TLine> TDPOCA;
      typedef typename KTRAJ::PDATA PDATA; // forward derivative type
      typedef typename KTRAJ::PDer PDer; // forward derivative type
      virtual unsigned nDOF() const override { return thit_.isActive() ? thit_.nDOF() : 0; }
      THit const& hit() const { return thit_; }
      virtual bool update(PKTRAJ const& ref)  override;
      virtual bool isActive() const override { return thit_.isActive(); }
      virtual double time() const override { return tdpoca_.poca0().T(); } // time on the main trajectory
      virtual double chisq(PDATA const& pars) const override;
      virtual ~KKHit(){}
      // construct from a hit and reference trajectory
      KKHit(THit const& thit, KTRAJ const& reftraj);
      KKHit(THit const& thit, PKTRAJ const& reftraj);
      Residual const& resid() const { return rresid_; }
      TDPOCA const& poca() const { return tdpoca_; }
      // residual derivatives WRT local trajectory parameters
      PDer dRdP() const { return rresid_.dRdD()*tdpoca_.dDdP() + rresid_.dRdT()*tdpoca_.dTdP(); }
    private:
    // helpers
      void setWeight();
      THit const& thit_; // hit used for this constraint
      TDPOCA tdpoca_; // POCA between the reference trajectory and the sensor trajectory
      Residual rresid_; // residual between the hit and reference
  };

  template<class KTRAJ> KKHit<KTRAJ>::KKHit(THit const& thit, PKTRAJ const& reftraj) : thit_(thit) ,
    tdpoca_(reftraj,thit_.sensorTraj()) {
    KKEFF::setRefTraj(reftraj.nearestPiece(tdpoca_.poca0().T()));
    thit_.resid(tdpoca_,rresid_);
    setWeight();
  }
  
  template<class KTRAJ> bool KKHit<KTRAJ>::update(PKTRAJ const& ref) {
  // update the hit state to the previous TPOCA
    thit_.update(tdpoca_);
    tdpoca_ = TDPOCA(ref,thit_.sensorTraj());
    KKEFF::setRefTraj(ref);
    thit_.resid(tdpoca_,rresid_);
    setWeight();
    KKEFF::setStatus(TDir::forwards,KKEFF::unprocessed);
    KKEFF::setStatus(TDir::backwards,KKEFF::unprocessed);
    return true;
  }

  template<class KTRAJ> void KKHit<KTRAJ>::setWeight() {
    // convert derivatives to a Nx1 matrix (for root)
    ROOT::Math::SMatrix<double,KTRAJ::NParams(),1> dRdPM;
    auto drdp = dRdP();
    dRdPM.Place_in_col(drdp,0,0);
    // convert the variance into a 1X1 matrix
    ROOT::Math::SMatrix<double, 1,1, ROOT::Math::MatRepSym<double,1> > RVarM;
    RVarM(0,0) = 1.0/rresid_.residVar();
    // expand these into the weight matrix
    KKWeight<KTRAJ>::weight_.weightMat() = ROOT::Math::Similarity(dRdPM,RVarM);
    // reference weight vector from reference parameters
    auto refvec = KKWeight<KTRAJ>::weight().weightMat()*KKEFF::referenceTraj().params().parameters();
    // translate residual value into weight vector WRT the reference parameters
    auto delta = drdp*rresid_.resid()/rresid_.residVar();
    // add change WRT reference; sign convention reflects resid = measurement - prediction
    KKWeight<KTRAJ>::weight_.weightVec() = refvec + delta;
  }

  template<class KTRAJ> double KKHit<KTRAJ>::chisq(PDATA const& pars) const {
    // compute the difference between these parameters and the reference parameters
    auto dpvec = pars.parameters() - KKEFF::referenceTraj().params().parameters();
    // use the differnce to 'correct' the reference residual to be WRT these parameters
    double newres = resid().resid() - ROOT::Math::Dot(dpvec,dRdP()); // check sign FIXME!
    // project the parameter covariance into a residual space variance (adding the intrinsic variance)
    double rvar = ROOT::Math::Similarity(dRdP(),pars.covariance()) + resid().residVar();
    // chisquared is the residual squared divided by the variance
    double chisq = newres*newres/rvar;
    return chisq;
  }

}
#endif
