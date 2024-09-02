#ifndef KinKal_Material_hh
#define KinKal_Material_hh
//
// Class to describe effect of a particle passing through discrete material on the fit (ie material transport)
// This effect adds no information content, just noise, and is KKEFF::processed in params space
//
#include "KinKal/Fit/Effect.hh"
#include "KinKal/Detector/ElementXing.hh"
#include "KinKal/General/TimeDir.hh"
#include <iostream>
#include <stdexcept>
#include <array>
#include <ostream>

namespace KinKal {
  template<class KTRAJ> class Material : public Effect<KTRAJ> {
    public:
      using KKEFF = Effect<KTRAJ>;
      using PTRAJ = ParticleTrajectory<KTRAJ>;
      using EXING = ElementXing<KTRAJ>;
      using EXINGPTR = std::shared_ptr<EXING>;
      double time() const override { return exingptr_->time();}
      bool active() const override { return  exingptr_->active(); }
      void process(FitState& kkdata,TimeDir tdir) override;
      void updateState(MetaIterConfig const& miconfig,bool first) override;
      void updateConfig(Config const& config) override {}
      void append(PTRAJ& fit,TimeDir tdir) override;
      void extrapolate(PTRAJ& fit,TimeDir tdir) override;
      void updateReference(PTRAJ const& ptraj) override;
      Chisq chisq(Parameters const& pdata) const override { return Chisq();}
      void print(std::ostream& ost=std::cout,int detail=0) const override;
      virtual ~Material(){}
      // create from the material and a trajectory
      Material(EXINGPTR const& exingptr, PTRAJ const& ptraj);
      // accessors
      auto const& elementXing() const { return *exingptr_; }
      auto const& elementXingPtr() const { return exingptr_; }
      auto const& referenceTrajectory() const { return exingptr_->referenceTrajectory(); }
    private:
      EXINGPTR exingptr_; // element crossing for this effect
      Weights nextwt_; // cache of weight forwards of this effect.
  };

  template<class KTRAJ> Material<KTRAJ>::Material(EXINGPTR const& exingptr, PTRAJ const& ptraj) : exingptr_(exingptr) {
    updateReference(ptraj);
  }

  template<class KTRAJ> void Material<KTRAJ>::process(FitState& kkdata,TimeDir tdir) {
    // The assymetry in the cache processing WRT fit update implements the convention that the cache is forwards in time of this effect, and insures the effect is only included once
    if(exingptr_->active()){
      if(tdir == TimeDir::forwards) {
        kkdata.append(exingptr_->params(),tdir);
        nextwt_ += kkdata.wData();
      } else {
        nextwt_ += kkdata.wData();
        kkdata.append(exingptr_->params(),tdir);
      }
    }
  }

  template<class KTRAJ> void Material<KTRAJ>::updateState(MetaIterConfig const& miconfig,bool first) {
    // update the ElementXing
    exingptr_->updateState(miconfig,first);
    // reset the cached weights
    nextwt_ = Weights();
  }

  template<class KTRAJ> void Material<KTRAJ>::append(PTRAJ& ptraj,TimeDir tdir) {
    if(exingptr_->active()){
      // create a trajectory piece from the cached weight
      double etime = this->time();
      // make sure this effect is appendable
      if( (tdir == TimeDir::forwards && etime < ptraj.back().range().begin()) ||
          (tdir == TimeDir::backwards && etime > ptraj.front().range().end()) )
        throw std::invalid_argument("New piece overlaps existing");
      KTRAJ newpiece = (tdir == TimeDir::forwards) ? ptraj.back() : ptraj.front();
      // make sure the range includes the transit time
      newpiece.range() = (tdir == TimeDir::forwards) ? TimeRange(etime,std::max(ptraj.range().end(),etime+exingptr_->transitTime())) :
        TimeRange(std::min(ptraj.range().begin(),etime-exingptr_->transitTime()),etime);
      // invert the cached weight to define the parameters
      newpiece.params() = Parameters(nextwt_);
      if( tdir == TimeDir::forwards){
        ptraj.append(newpiece);
      } else {
        // Since the cache was forwards of this site, we have to apply the effect of this material xing to the parameters.
        auto temp = exingptr_->params();
        temp.parameters() *= -1; // reverse the sign of the parameter change
        newpiece.params() += temp;
        ptraj.prepend(newpiece);
      }
    }
  }

  template<class KTRAJ> void Material<KTRAJ>::extrapolate(PTRAJ& ptraj,TimeDir tdir) {
  // make sure the traj can be extrapolated
    if((tdir == TimeDir::forwards && ptraj.back().range().begin() > time()) ||
        (tdir == TimeDir::backwards && ptraj.front().range().end() < time()) )
      throw std::invalid_argument("Material: Can't append piece");
    auto ktrajptr = ptraj.nearestTraj(time());
    // convert parameters at this material to momentum space representation
    auto pstate = ktrajptr->stateEstimate(time());
    auto bnom = ktrajptr->bnom();
    // change in momentum due to the associated element Xing
    SVEC3 dmom;
    SMAT dmomvar;
    exingptr_->momentumChange(dmom,dmomvar);
    // expand these to the full particle state, with null position and cross-covariance
    DMAT dstatevar; dstatevar.Place_at(dmomvar,3,3); // lower right corner for momentum
    auto psvar = pstate.stateCovariance() + dstatevar;
    SVEC6 dstate; dstate.Place_at(dmom,3);
    // update the pstate accordingly. momentum vector change depends on the time direction, but not the covariance
    if( tdir == TimeDir::forwards){
      // add the material effect
      auto psvec = pstate.state() + dstate;
      ParticleStateEstimate newpstate(psvec,psvar,time(),pstate.mass(),pstate.charge());
      TimeRange range(time(),ptraj.range().end());
      KTRAJ newpiece(newpstate, bnom, range);
      ptraj.append(newpiece);
    } else {
      auto psvec = pstate.state() - dstate;
      ParticleStateEstimate newpstate(psvec,psvar,time(),pstate.mass(),pstate.charge());
      TimeRange range(ptraj.range().begin(),time());
      KTRAJ newpiece(newpstate, bnom, range);
      ptraj.prepend(newpiece);
    }
  }

  template<class KTRAJ> void Material<KTRAJ>::updateReference(PTRAJ const& ptraj) {
    exingptr_->updateReference(ptraj);
  }

  template<class KTRAJ> void Material<KTRAJ>::print(std::ostream& ost,int detail) const {
    ost << "Material " << static_cast<Effect<KTRAJ>const&>(*this);
    ost << " ElementXing ";
    exingptr_->print(ost,detail);
    if(detail >3){
      ost << " forward cache ";
      nextwt_.print(ost,detail);
    }
  }

  template <class KTRAJ> std::ostream& operator <<(std::ostream& ost, Material<KTRAJ> const& kkmat) {
    kkmat.print(ost,0);
    return ost;
  }
}
#endif
