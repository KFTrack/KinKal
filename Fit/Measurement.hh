#ifndef KinKal_Measurement_hh
#define KinKal_Measurement_hh
//
//  class represeting a constraint on the fit parameters due external information (typically a measurement).
//  Used as part of the kinematic Kalman fit
//
#include "KinKal/Fit/Effect.hh"
#include "KinKal/Trajectory/ParticleTrajectory.hh"
#include "KinKal/Detector/Hit.hh"
#include <ostream>
#include <memory>

namespace KinKal {
  template <class KTRAJ> class Measurement : public Effect<KTRAJ> {
    public:
      using KKEFF = Effect<KTRAJ>;
      using PKTRAJ = ParticleTrajectory<KTRAJ>;
      using HIT = Hit<KTRAJ>;
      using HITPTR = std::shared_ptr<HIT>;

      Chisq chisq(Parameters const& pdata) const override;
      Chisq chisq() const override;
      void update(PKTRAJ const& pktraj) override;
      void update(PKTRAJ const& pktraj, MetaIterConfig const& miconfig) override;
      void update(Config const& config) override {}
      void process(FitState& kkdata,TimeDir tdir) override;
      bool active() const override { return hit_->active(); }
      double time() const override { return hit_->time(); }
      void print(std::ostream& ost=std::cout,int detail=0) const override;
      virtual ~Measurement(){}
      // local functions
      // construct from a hit and reference trajectory
      Measurement(HITPTR const& hit, PKTRAJ const& reftraj);
      // the unbiased parameters are the fit parameters not including the information content of this effect
      Parameters unbiasedParameters() const;
      // access the contents
      HITPTR const& hit() const { return hit_; }
      Weights const& weightCache() const { return wcache_; }
      Weights const& hitWeight() const { return hitwt_; }
    private:
      HITPTR hit_ ; // hit used for this constraint
      Weights wcache_; // sum of processing weights in opposite directions, excluding this hit's information. used to compute unbiased parameters and chisquared
      Weights hitwt_; // weight representation of the hits constraint
      double vscale_; // variance factor due to annealing 'temperature'
  };

  template<class KTRAJ> Measurement<KTRAJ>::Measurement(HITPTR const& hit, PKTRAJ const& reftraj) : hit_(hit), vscale_(1.0) {
    update(reftraj);
  }

  template<class KTRAJ> void Measurement<KTRAJ>::process(FitState& kkdata,TimeDir tdir) {
    // direction is irrelevant for processing hits
    if(this->active()){
      // cache the processing weights, adding both processing directions
      wcache_ += kkdata.wData();
      // add this effect's information
      kkdata.append(hitwt_);
    }
    KKEFF::setState(tdir,KKEFF::processed);
  }

  template<class KTRAJ> void Measurement<KTRAJ>::update(PKTRAJ const& pktraj) {
    // reset the processing cache
    wcache_ = Weights();
    // update the hit
    hit_->update(pktraj);
    // get the weight from the hit
    hitwt_ = hit_->weight();
    // scale weight for the temp
    hitwt_ *= 1.0/vscale_;
    // ready for processing!
    KKEFF::updateState();
  }

  template<class KTRAJ> void Measurement<KTRAJ>::update(PKTRAJ const& pktraj, MetaIterConfig const& miconfig) {
    // reset the annealing temp and hit precision
    vscale_ = miconfig.varianceScale();
    // update the hit's internal state; the actual update depends on the hit
    hit_->update(pktraj,miconfig );
    // update the state of this object
    update(pktraj);
  }

  template<class KTRAJ> Chisq Measurement<KTRAJ>::chisq(Parameters const& pdata) const {
    return hit_->chisq(pdata);
  }

  template<class KTRAJ> Chisq Measurement<KTRAJ>::chisq() const {
    if(this->active()) {
      return hit_->chisq(unbiasedParameters());
    } else
      return Chisq();
  }

  template<class KTRAJ> Parameters Measurement<KTRAJ>::unbiasedParameters() const {
    // this function can't be called on an unprocessed effect
    if( !KKEFF::wasProcessed(TimeDir::forwards) || !KKEFF::wasProcessed(TimeDir::backwards))
      throw  std::invalid_argument("Can't compute unbiased parameters for unprocessed constraint");
    // Invert the cache to get unbiased parameters at this constraint
    return Parameters(wcache_);
  }

  template <class KTRAJ> void Measurement<KTRAJ>::print(std::ostream& ost, int detail) const {
    ost << "Measurement " << static_cast<Effect<KTRAJ> const&>(*this) << std::endl;
    if(detail > 0){
      hit_->print(ost,detail);
      ost << " Measurement Weight " << hitwt_ << std::endl;
    }
  }

  template <class KTRAJ> std::ostream& operator <<(std::ostream& ost, Measurement<KTRAJ> const& kkhit) {
    kkhit.print(ost,0);
    return ost;
  }

}
#endif
