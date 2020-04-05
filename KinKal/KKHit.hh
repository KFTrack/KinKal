#ifndef KinKal_KKHit_hh
#define KinKal_KKHit_hh
//
//  class to use information from a hit in the Kinematic fit.
//  Used as part of the kinematic Kalman fit
//
#include "KinKal/KKWEff.hh"
#include "KinKal/PKTraj.hh"
#include "KinKal/THit.hh"
#include "KinKal/TPocaBase.hh"
#include "KinKal/Residual.hh"
#include <ostream>
#include <memory>

namespace KinKal {
  class TTraj;
  template <class KTRAJ> class KKHit : public KKWEff<KTRAJ> {
    public:
      typedef KKEff<KTRAJ> KKEFF;
      typedef KKWEff<KTRAJ> KKWEFF;
      typedef PKTraj<KTRAJ> PKTRAJ;
      typedef THit<KTRAJ> THIT;
      typedef std::shared_ptr<THIT> THITPTR;
      typedef typename KTRAJ::PDATA PDATA; // forward derivative type
      typedef typename KTRAJ::PDATA::DVEC DVEC; // forward derivative type
      typedef typename KTRAJ::PDER PDER; // forward derivative type
      virtual unsigned nDOF() const override { return thit_->isActive() ? thit_->nDOF() : 0; }
      virtual bool update(PKTRAJ const& ref)  override;
      virtual bool isActive() const override { return thit_->isActive(); }
      virtual double time() const override { return thit_->poca().particleToca(); } // time on the particle trajectory
      virtual double chisq(PDATA const& pars) const override;
      virtual void print(std::ostream& ost=std::cout,int detail=0) const override;
      virtual ~KKHit(){}
      // construct from a hit and reference trajectory
      KKHit(THITPTR const& thit, PKTRAJ const& reftraj);
      Residual const& refResid() const { return rresid_; }
      // residual derivatives WRT local trajectory parameters
      PDER dRdP() const { return rresid_.dRdD()*thit_->dDdP() + rresid_.dRdT()*thit_->dTdP(); }
      // accessors
      THITPTR const& tHit() const { return thit_; }
      TPocaBase const& poca() const { return thit_->poca(); }
    private:
    // helpers
      void setWeight();
      THITPTR thit_ ; // hit used for this constraint
      DVEC ref_; // reference parameters
      Residual rresid_; // residual between the hit and reference
  };

  template<class KTRAJ> KKHit<KTRAJ>::KKHit(THITPTR const& thit, PKTRAJ const& reftraj) : thit_(thit) {
    update(reftraj);
  }
  
  template<class KTRAJ> bool KKHit<KTRAJ>::update(PKTRAJ const& ref) {
    thit_->update(ref);
    ref_ = ref.nearestPiece(thit_->poca().particleToca()).params().parameters();
    thit_->resid(rresid_);
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
    // add annealing temperature to weight FIXME!
    RVarM(0,0) = 1.0/rresid_.residVar();
    // expand these into the weight matrix
    KKWEFF::wdata_.weightMat() = ROOT::Math::Similarity(dRdPM,RVarM);
    KKWEFF::wdata_.setStatus(PDATA::valid);
    // reference weight vector from reference parameters
    KKWEFF::wdata_.weightVec() = KKWEFF::wData().weightMat()*ref_;
    // translate residual value into weight vector WRT the reference parameters
    // and add change WRT reference; sign convention reflects resid = measurement - prediction
    KKWEFF::wdata_.weightVec() += drdp*rresid_.resid()/rresid_.residVar();
  }

  template<class KTRAJ> double KKHit<KTRAJ>::chisq(PDATA const& pdata) const {
    // compute the difference between these parameters and the reference parameters
    DVEC dpvec = pdata.parameters() - ref_; 
    // use the differnce to 'correct' the reference residual to be WRT these parameters
    double newres = refResid().resid() - ROOT::Math::Dot(dpvec,dRdP()); 
    // project the parameter covariance into a residual space variance (adding the intrinsic variance)
    double rvar = ROOT::Math::Similarity(dRdP(),pdata.covariance()) + refResid().residVar();
    // chisquared is the residual squared divided by the variance
    double chisq = newres*newres/rvar;
    return chisq;
  }

  template <class KTRAJ> void KKHit<KTRAJ>::print(std::ostream& ost, int detail) const {
    ost << "KKHit " << static_cast<KKEff<KTRAJ> const&>(*this) << 
      " doca " << poca().doca() << " resid " << refResid() << std::endl;
    if(detail > 0){
      thit_->print(ost,detail);    
      ost << "Reference " << ref_ << std::endl;

    }
  }

  template <class KTRAJ> std::ostream& operator <<(std::ostream& ost, KKHit<KTRAJ> const& kkhit) {
    kkhit.print(ost,0);
    return ost;
  }

}
#endif
