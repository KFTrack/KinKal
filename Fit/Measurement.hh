#ifndef KinKal_Measurement_hh
#define KinKal_Measurement_hh
//
//  class wrapping a detector measurement for the Kinematic fit.
//  Used as part of the kinematic Kalman fit
//
#include "KinKal/Fit/Effect.hh"
#include "KinKal/Trajectory/ParticleTrajectory.hh"
#include "KinKal/Detector/DetectorHit.hh"
#include <ostream>
#include <memory>

namespace KinKal {
  template <class KTRAJ> class Measurement : public Effect<KTRAJ> {
    public:
      using KKEFF = Effect<KTRAJ>;
      using PKTRAJ = ParticleTrajectory<KTRAJ>;
      using DHIT = DetectorHit<KTRAJ>;
      using DHITPTR = std::shared_ptr<DHIT>;
      
      unsigned nDOF() const override { return hit_->isActive() ? hit_->nDOF() : 0; }
      double fitChi() const override; 
      double chisq(Parameters const& pdata) const override;
      void update(PKTRAJ const& pktraj) override;
      void update(PKTRAJ const& pktraj, MetaIterConfig const& miconfig) override;
      void process(FitState& kkdata,TimeDir tdir) override;
      bool isActive() const override { return hit_->isActive(); }
      double time() const override { return hit_->time(); }
      void print(std::ostream& ost=std::cout,int detail=0) const override;
      virtual ~Measurement(){}
      // local functions
      // construct from a hit and reference trajectory
      Measurement(DHITPTR const& hit, PKTRAJ const& reftraj,double precision=1e-6);
      // interface for reduced residual
      double chi(Parameters const& pdata) const;
      // accessors
      DHITPTR const& detectorHit() const { return hit_; }
      Weights const& weightCache() const { return wcache_; }
      Weights const& measurementWeight() const { return mwt_; }
      double precision() const { return precision_; }
      // compute the reduced residual
    private:
      DHITPTR hit_ ; // hit used for this constraint
      Weights wcache_; // sum of processing weights in opposite directions, excluding this hit's information. used to compute chisquared and reduced residuals
      Weights mwt_; // wdata representation of this effect's constraint/measurement
      double vscale_; // variance factor due to annealing 'temperature'
      double precision_; // precision used in TCA calcuation
  };

  template<class KTRAJ> Measurement<KTRAJ>::Measurement(DHITPTR const& hit, PKTRAJ const& reftraj,double precision) : hit_(hit), vscale_(1.0), precision_(precision) {
    update(reftraj);
  }
 
  template<class KTRAJ> void Measurement<KTRAJ>::process(FitState& kkdata,TimeDir tdir) {
    // direction is irrelevant for processing measurements 
    if(this->isActive()){
      // cache the processing weights, adding both processing directions
      wcache_ += kkdata.wData();
      // add this effect's information
      kkdata.append(mwt_);
    }
    KKEFF::setStatus(tdir,KKEFF::processed);
  }

  template<class KTRAJ> void Measurement<KTRAJ>::update(PKTRAJ const& pktraj) {
    // reset the processing cache
    wcache_ = Weights();
    // update the hit
    hit_->update(pktraj);
    // get the weight from the hit 
    mwt_ = hit_->weight();
    // scale weight for the temp
    mwt_ *= 1.0/vscale_;
    // ready for processing!
    KKEFF::updateStatus();
  }

  template<class KTRAJ> void Measurement<KTRAJ>::update(PKTRAJ const& pktraj, MetaIterConfig const& miconfig) {
    // reset the annealing temp and measurement precision
    vscale_ = miconfig.varianceScale();
    precision_ = miconfig.tprec_;
    // update the measurement's internal state; this depends on the configuration parameters
    if(miconfig.updatehits_)hit_->update(pktraj,miconfig );
    // update the state of this object
    update(pktraj);
  }

  template<class KTRAJ> double Measurement<KTRAJ>::chisq(Parameters const& pdata) const {
    double retval(0.0);
    if(this->isActive()) {
      double chi = hit_->chi(pdata);
      // correct for current variance scaling
      retval = chi*chi/vscale_;
    }
    return retval;
  }

  template<class KTRAJ> double Measurement<KTRAJ>::fitChi() const {
    double retval(0.0);
    if(this->isActive() && KKEFF::wasProcessed(TimeDir::forwards) && KKEFF::wasProcessed(TimeDir::backwards)) {
    // Invert the cache to get unbiased parameters at this hit
      Parameters unbiased(wcache_);
      retval = hit_->chi(unbiased)/vscale_;
    }
    return retval;
  }

  template <class KTRAJ> void Measurement<KTRAJ>::print(std::ostream& ost, int detail) const {
    ost << "Measurement " << static_cast<Effect<KTRAJ> const&>(*this) << std::endl;
    if(detail > 0){
      hit_->print(ost,detail);    
      ost << " Measurement Weight " << mwt_ << std::endl;
    }
  }

  template <class KTRAJ> std::ostream& operator <<(std::ostream& ost, Measurement<KTRAJ> const& kkhit) {
    kkhit.print(ost,0);
    return ost;
  }

}
#endif
