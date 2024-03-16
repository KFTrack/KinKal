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
      using PTRAJ = ParticleTrajectory<KTRAJ>;
      using HIT = Hit<KTRAJ>;
      using HITPTR = std::shared_ptr<HIT>;
      // Effect Interface
      double time() const override { return hit_->time(); }
      bool active() const override { return hit_->active(); }
      void process(FitState& kkdata,TimeDir tdir) override;
      void updateState(MetaIterConfig const& miconfig,bool first) override;
      void updateConfig(Config const& config) override {}
      void updateReference(PTRAJ const& ptraj) override;
      void append(PTRAJ& fit,TimeDir tdir) override;
      Chisq chisq(Parameters const& pdata) const override;
      void print(std::ostream& ost=std::cout,int detail=0) const override;
      virtual ~Measurement(){}
      // local functions
      // construct from a hit and reference trajectory
      Measurement(HITPTR const& hit,PTRAJ const& ptraj);
      // access the underlying hit
      HITPTR const& hit() const { return hit_; }
    private:
      HITPTR hit_ ; // hit used for this measurement
  };

  template<class KTRAJ> Measurement<KTRAJ>::Measurement(HITPTR const& hit,PTRAJ const& ptraj) : hit_(hit) {
    this->updateReference(ptraj);
  }

  template<class KTRAJ> void Measurement<KTRAJ>::process(FitState& kkdata,TimeDir tdir) {
    // add this effect's information. direction is irrelevant for processing hits
    if(this->active())kkdata.append(hit_->weight());
  }

  template<class KTRAJ> void Measurement<KTRAJ>::updateState(MetaIterConfig const& miconfig,bool first) {
    // update the hit's internal state; the actual update depends on the hit
    hit_->updateState(miconfig,first);
  }

  template<class KTRAJ> void Measurement<KTRAJ>::append(PTRAJ& ptraj,TimeDir tdir) {
  }

  template<class KTRAJ> void Measurement<KTRAJ>::updateReference(PTRAJ const& ptraj) {
    hit_->updateReference(ptraj.nearestTraj(hit_->time()));
  }

  template<class KTRAJ> Chisq Measurement<KTRAJ>::chisq(Parameters const& pdata) const {
    return hit_->chisq(pdata);
  }

  template<class KTRAJ> void Measurement<KTRAJ>::print(std::ostream& ost, int detail) const {
    ost << "Measurement " << static_cast<Effect<KTRAJ> const&>(*this) << std::endl;
    if(detail > 0){
      hit_->print(ost,detail);
    }
  }

  template <class KTRAJ> std::ostream& operator <<(std::ostream& ost, Measurement<KTRAJ> const& measurement) {
    measurement.print(ost,0);
    return ost;
  }

}
#endif
