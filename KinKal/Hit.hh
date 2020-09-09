#ifndef KinKal_Hit_hh
#define KinKal_Hit_hh
//
//  class to use information from a hit in the Kinematic fit.
//  Used as part of the kinematic Kalman fit
//
#include "KinKal/Effect.hh"
#include "KinKal/ParticleTrajectory.hh"
#include "KinKal/DetectorHit.hh"
#include "KinKal/Residual.hh"
#include <ostream>
#include <memory>

namespace KinKal {
  template <class KTRAJ> class Hit : public Effect<KTRAJ> {
    public:
      using KKEFF = Effect<KTRAJ>;
      using PKTRAJ = ParticleTrajectory<KTRAJ>;
      using THIT = DetectorHit<KTRAJ>;
      using THITPTR = std::shared_ptr<THIT>;
      
      virtual unsigned nDOF() const override { return thit_->isActive() ? thit_->nDOF() : 0; }
      virtual double fitChi() const override; 
      virtual double chisq(Parameters const& pdata) const override{ double chival = chi(pdata); return chival*chival; } 
      virtual void update(PKTRAJ const& pktraj)  override;
      virtual void update(PKTRAJ const& pktraj, MetaIterConfig const& miconfig) override;
      virtual void process(FitData& kkdata,TimeDir tdir) override;
      virtual bool isActive() const override { return thit_->isActive(); }
      virtual double time() const override { return rresid_.time(); } // time on the particle trajectory
      virtual void print(std::ostream& ost=std::cout,int detail=0) const override;
      virtual ~Hit(){}
      // local functions
      void updateCache(PKTRAJ const& pktraj);
      // construct from a hit and reference trajectory
      Hit(THITPTR const& thit, PKTRAJ const& reftraj,double precision=1e-6);
      // interface for reduced residual
      double chi(Parameters const& pdata) const;
      // accessors
      THITPTR const& tHit() const { return thit_; }
      Residual const& refResid() const { return rresid_; }
      Parameters const& refParams() const { return ref_; }
      Weights const& weightCache() const { return wcache_; }
      double precision() const { return precision_; }
      // compute the reduced residual
    private:
      THITPTR thit_ ; // hit used for this constraint
      Parameters ref_; // reference parameters
      Weights wcache_; // sum of processing weights in opposite directions, excluding this hit's information. used to compute chisquared and reduced residuals
      Weights hiteff_; // wdata representation of this effect's constraint/measurement
      Residual rresid_; // residuals for this reference and hit
      double vscale_; // variance factor due to annealing 'temperature'
      double precision_; // precision used in TPOCA calcuation
  };

  template<class KTRAJ> Hit<KTRAJ>::Hit(THITPTR const& thit, PKTRAJ const& reftraj,double precision) : thit_(thit), vscale_(1.0), precision_(precision) {
    update(reftraj);
  }
 
  template<class KTRAJ> void Hit<KTRAJ>::process(FitData& kkdata,TimeDir tdir) {
    // direction is irrelevant for adding information
    if(this->isActive()){
      // cache the processing weights, adding both processing directions
      wcache_ += kkdata.wData();
      // add this effect's information
      kkdata.append(hiteff_);
    }
    KKEFF::setStatus(tdir,KKEFF::processed);
  }

  template<class KTRAJ> void Hit<KTRAJ>::update(PKTRAJ const& pktraj) {
    // compute residual and derivatives from hit using reference parameters
    thit_->resid(pktraj, rresid_, precision_);
    updateCache(pktraj);
  }

  template<class KTRAJ> void Hit<KTRAJ>::update(PKTRAJ const& pktraj, MetaIterConfig const& miconfig) {
    // reset the annealing temp
    vscale_ = miconfig.varianceScale();
    precision_ = miconfig.tprec_;
    // update the hit internal state; this can depend on specific configuration parameters
    if(miconfig.updatehits_)
      thit_->update(pktraj,miconfig, rresid_);
    else
      thit_->resid(pktraj, rresid_, precision_);
    // update the state of this object
    updateCache(pktraj);
  }

  template<class KTRAJ> void Hit<KTRAJ>::updateCache(PKTRAJ const& pktraj) {
    // reset the processing cache
    wcache_ = Weights();
    // scale resid variance by temp normalization
    double tvar = rresid_.variance()*vscale_; 
    ref_ = pktraj.nearestPiece(rresid_.time()).params();
    // convert derivatives to a Nx1 matrix (for root)
    ROOT::Math::SMatrix<double,NParams(),1> dRdPM;
    dRdPM.Place_in_col(rresid_.dRdP(),0,0);
    // convert the variance into a 1X1 matrix
    ROOT::Math::SMatrix<double, 1,1, ROOT::Math::MatRepSym<double,1>> RVarM;
    // weight by inverse variance
    RVarM(0,0) = 1.0/tvar;
    // expand these into the weight matrix
    hiteff_.weightMat() = ROOT::Math::Similarity(dRdPM,RVarM);
    // translate residual value into weight vector WRT the reference parameters
    // sign convention reflects resid = measurement - prediction
    hiteff_.weightVec() = hiteff_.weightMat()*ref_.parameters() + rresid_.dRdP()*rresid_.value()/tvar;
    KKEFF::updateStatus();
  }

  template<class KTRAJ> double Hit<KTRAJ>::fitChi() const {
    double retval(0.0);
    if(this->isActive() && KKEFF::wasProcessed(TimeDir::forwards) && KKEFF::wasProcessed(TimeDir::backwards)) {
    // Invert the cache to get unbiased parameters at this hit
      Parameters unbiased(wcache_);
      retval = chi(unbiased);
    }
    return retval;
  }

  template<class KTRAJ> double Hit<KTRAJ>::chi(Parameters const& pdata) const {
    double retval(0.0);
    if(this->isActive()) {
      // compute the difference between these parameters and the reference parameters
      DVEC dpvec = pdata.parameters() - ref_.parameters(); 
      // use the differnce to 'correct' the reference residual to be WRT these parameters
      double uresid = rresid_.value() - ROOT::Math::Dot(dpvec,rresid_.dRdP());
      // project the parameter covariance into a residual space variance
      double rvar = ROOT::Math::Similarity(rresid_.dRdP(),pdata.covariance());
      // add the measurement variance, scaled by the current temperature normalization
      rvar +=  rresid_.variance()*vscale_;
      // chi is the ratio of these
      retval = uresid/sqrt(rvar);
    }
    return retval;
  }

  template <class KTRAJ> void Hit<KTRAJ>::print(std::ostream& ost, int detail) const {
    ost << "Hit " << static_cast<Effect<KTRAJ> const&>(*this) << " resid " << refResid()  << std::endl;
    if(detail > 0){
      thit_->print(ost,detail);    
      ost << "Reference " << ref_ << std::endl;
    }
  }

  template <class KTRAJ> std::ostream& operator <<(std::ostream& ost, Hit<KTRAJ> const& kkhit) {
    kkhit.print(ost,0);
    return ost;
  }

}
#endif
