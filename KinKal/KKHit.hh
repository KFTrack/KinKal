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
  template <class KTRAJ> class KKHit : public KKWEff<KTRAJ> {
    public:
      typedef KKEff<KTRAJ> KKEFF;
      typedef KKWEff<KTRAJ> KKWEFF;
      typedef PKTraj<KTRAJ> PKTRAJ;
      typedef THit<KTRAJ> THIT;
      typedef Residual<KTRAJ> RESIDUAL;
      typedef std::shared_ptr<THIT> THITPTR;
      typedef typename KTRAJ::PDATA PDATA; // forward derivative type
      typedef TData<PDATA::PDim()> TDATA;
      typedef typename KTRAJ::PDATA::DVEC DVEC; // forward derivative type
      typedef typename KTRAJ::PDER PDER; // forward derivative type
      virtual unsigned nDOF() const override { return thit_->isActive() ? thit_->nDOF() : 0; }
      virtual bool update(PKTRAJ const& ref)  override;
      virtual bool isActive() const override { return thit_->isActive(); }
      virtual float time() const override { return rresid_.residTime(); } // time on the particle trajectory
      virtual float chisq(PDATA const& pars) const override;
      virtual void print(std::ostream& ost=std::cout,int detail=0) const override;
      virtual ~KKHit(){}
      // construct from a hit and reference trajectory
      KKHit(THITPTR const& thit, PKTRAJ const& reftraj);
      // interface for reduced residual needed here FIXME
      // accessors
      THITPTR const& tHit() const { return thit_; }
      RESIDUAL const& refResid() const { return rresid_; }
      DVEC const& refParams() const { return ref_; }
    private:
      THITPTR thit_ ; // hit used for this constraint
      DVEC ref_; // reference parameters
      RESIDUAL rresid_; // residuals for this reference and hit
  };

  template<class KTRAJ> KKHit<KTRAJ>::KKHit(THITPTR const& thit, PKTRAJ const& reftraj) : thit_(thit) {
    update(reftraj);
  }
  
  template<class KTRAJ> bool KKHit<KTRAJ>::update(PKTRAJ const& pktraj) {
    // compute residual and derivatives from hit using reference parameters
    thit_->resid(pktraj, rresid_);
    ref_ = pktraj.nearestPiece(rresid_.residTime()).params().parameters();
    // convert derivatives to a Nx1 matrix (for root)
    ROOT::Math::SMatrix<double,KTRAJ::NParams(),1> dRdPM;
    dRdPM.Place_in_col(rresid_.dRdP(),0,0);
    // convert the variance into a 1X1 matrix
    ROOT::Math::SMatrix<double, 1,1, ROOT::Math::MatRepSym<double,1> > RVarM;
    // add annealing temperature to weight FIXME!
    RVarM(0,0) = 1.0/rresid_.residVar();
    // expand these into the weight matrix
    KKWEFF::wdata_.weightMat() = ROOT::Math::Similarity(dRdPM,RVarM);
    KKWEFF::wdata_.setStatus(TDATA::valid);
    // reference weight vector from reference parameters
    KKWEFF::wdata_.weightVec() = KKWEFF::wData().weightMat()*ref_;
    // translate residual value into weight vector WRT the reference parameters
    // and add change WRT reference; sign convention reflects resid = measurement - prediction
    KKWEFF::wdata_.weightVec() += rresid_.dRdP()*rresid_.resid()/rresid_.residVar();
    KKEffBase::updateStatus();
    return true;
  }

  template<class KTRAJ> float KKHit<KTRAJ>::chisq(PDATA const& pdata) const {
    // compute the difference between these parameters and the reference parameters
    DVEC dpvec = pdata.parameters() - ref_; 
    // use the differnce to 'correct' the reference residual to be WRT these parameters
    float newres = rresid_.resid() - ROOT::Math::Dot(dpvec,rresid_.dRdP()); 
    // project the parameter covariance into a residual space variance (adding the intrinsic variance)
    float rvar = ROOT::Math::Similarity(rresid_.dRdP(),pdata.covariance()) + rresid_.residVar();
    // chisquared is the residual squared divided by the variance
    float chisq = newres*newres/rvar;
    return chisq;
  }

  template <class KTRAJ> void KKHit<KTRAJ>::print(std::ostream& ost, int detail) const {
    ost << "KKHit " << static_cast<KKEff<KTRAJ> const&>(*this) << " resid " << refResid() << std::endl;
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
