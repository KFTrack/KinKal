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
      double time() const override { return exing_->time();}
      bool active() const override { return  exing_->active(); }
      void process(FitState& kkdata,TimeDir tdir) override;
      void updateState(MetaIterConfig const& miconfig,bool first) override;
      void updateConfig(Config const& config) override {}
      void append(PTRAJ& fit,TimeDir tdir) override;
      void updateReference(PTRAJ const& ptraj) override;
      Chisq chisq(Parameters const& pdata) const override { return Chisq();}
      void print(std::ostream& ost=std::cout,int detail=0) const override;
      virtual ~Material(){}
      // create from the material and a trajectory
      Material(EXINGPTR const& dxing, PTRAJ const& ptraj);
      // accessors
      auto const& elementXing() const { return *exing_; }
      auto const& elementXingPtr() const { return exing_; }
      auto const& referenceTrajectory() const { return exing_->referenceTrajectory(); }
    private:
      EXINGPTR exing_; // element crossing for this effect
      Weights nextwt_; // cache of weight forwards of this effect.
  };

  template<class KTRAJ> Material<KTRAJ>::Material(EXINGPTR const& dxing, PTRAJ const& ptraj) : exing_(dxing) {
    updateReference(ptraj);
  }

  template<class KTRAJ> void Material<KTRAJ>::process(FitState& kkdata,TimeDir tdir) {
    if(exing_->active()){
      // cache for the forwards side of this effect
      // forwards
      if(tdir == TimeDir::forwards) {
        kkdata.append(exing_->params(),tdir);
        nextwt_ += kkdata.wData();
      } else {
        // backwards; note the append uses FORWARDS DIRECTION because the params funtcion already does the time ordering, using both would double-count
        nextwt_ += kkdata.wData();
        kkdata.append(exing_->params(),tdir);
      }
    }
  }

  template<class KTRAJ> void Material<KTRAJ>::updateState(MetaIterConfig const& miconfig,bool first) {
    // update the ElementXing
    exing_->updateState(miconfig,first);
    // reset the cached weights
    nextwt_ = Weights();
  }

  template<class KTRAJ> void Material<KTRAJ>::append(PTRAJ& ptraj,TimeDir tdir) {
    if(exing_->active()){
      // create a trajectory piece from the cached weight
      double etime = this->time();
      // make sure this effect is appendable
      if( (tdir == TimeDir::forwards && etime < ptraj.back().range().begin()) ||
          (tdir == TimeDir::backwards && etime > ptraj.front().range().end()) )
        throw std::invalid_argument("New piece overlaps existing");
      KTRAJ newpiece = (tdir == TimeDir::forwards) ? ptraj.back() : ptraj.front();
      // make sure the range includes the transit time
      newpiece.range() = (tdir == TimeDir::forwards) ? TimeRange(etime,std::max(ptraj.range().end(),etime+exing_->transitTime())) :
        TimeRange(std::min(ptraj.range().begin(),etime-exing_->transitTime()),etime);
      newpiece.params() = Parameters(nextwt_);
      if( tdir == TimeDir::forwards){
        ptraj.append(newpiece);
      } else {
        // Since the cache was forwards of this site, we have to apply the effect of this material xing to the parameters.
        auto temp = exing_->params();
        temp.parameters() *= -1; // reverse the sign of the parameter change
        newpiece.params() += temp;
        ptraj.prepend(newpiece);
      }
    }
    // update the xing
    if( tdir == TimeDir::forwards)
      exing_->updateReference(ptraj.backPtr());
    else
      exing_->updateReference(ptraj.frontPtr());
  }

  template<class KTRAJ> void Material<KTRAJ>::updateReference(PTRAJ const& ptraj) {
    exing_->updateReference(ptraj.nearestTraj(exing_->time()));
  }

  template<class KTRAJ> void Material<KTRAJ>::print(std::ostream& ost,int detail) const {
    ost << "Material " << static_cast<Effect<KTRAJ>const&>(*this);
    ost << " ElementXing ";
    exing_->print(ost,detail);
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
