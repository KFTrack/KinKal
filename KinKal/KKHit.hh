#ifndef KinKal_KKHit_hh
#define KinKal_KKHit_hh
//
//  class to use information from a hit in the Kinematic fit.
//  Used as part of the kinematic Kalman fit
//
#include "KinKal/KKWEff.hh"
#include "KinKal/PKTraj.hh"
#include "KinKal/THit.hh"
#include "KinKal/TPoca.hh"
#include "KinKal/TLine.hh"
#include "KinKal/Residual.hh"

namespace KinKal {
  class TTraj;
  template <class KTRAJ> class KKHit : public KKWEff<KTRAJ> {
    public:
      typedef KKEff<KTRAJ> KKEFF;
      typedef PKTraj<KTRAJ> PKTRAJ;
      typedef TDPoca<PKTRAJ,TLine> TDPOCA;
      typedef typename KTRAJ::PDATA PDATA; // forward derivative type
      typedef typename KTRAJ::PDER PDER; // forward derivative type
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
      Residual const& refResid() const { return rresid_; }
      TDPOCA const& poca() const { return tdpoca_; }
      // residual derivatives WRT local trajectory parameters
      PDER dRdP() const { return rresid_.dRdD()*tdpoca_.dDdP() + rresid_.dRdT()*tdpoca_.dTdP(); }
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
    KKEffBase::updateStatus();
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
    KKWEff<KTRAJ>::wdata_.weightMat() = ROOT::Math::Similarity(dRdPM,RVarM);
    KKWEff<KTRAJ>::wdata_.setStatus(PDATA::valid);
    // reference weight vector from reference parameters
    auto refvec = KKWEff<KTRAJ>::wData().weightMat()*KKEFF::refTraj().params().parameters();
    // translate residual value into weight vector WRT the reference parameters
    auto delta = drdp*rresid_.resid()/rresid_.residVar();
    // add change WRT reference; sign convention reflects resid = measurement - prediction
    KKWEff<KTRAJ>::wdata_.weightVec() = refvec + delta;
  }

  template<class KTRAJ> double KKHit<KTRAJ>::chisq(PDATA const& pdata) const {
    // compute the difference between these parameters and the reference parameters
    typename PDATA::DVEC dpvec = pdata.parameters() - KKEFF::refTraj().params().parameters();
    // use the differnce to 'correct' the reference residual to be WRT these parameters
    double newres = refResid().resid() - ROOT::Math::Dot(dpvec,dRdP()); 
    // project the parameter covariance into a residual space variance (adding the intrinsic variance)
    double rvar = ROOT::Math::Similarity(dRdP(),pdata.covariance()) + refResid().residVar();
    // chisquared is the residual squared divided by the variance
    double chisq = newres*newres/rvar;
    return chisq;
  }

}
#endif
