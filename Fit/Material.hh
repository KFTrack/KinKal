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
      using KTRAJPTR = std::shared_ptr<KTRAJ>;
      using EXING = ElementXing<KTRAJ>;
      using EXINGPTR = std::shared_ptr<EXING>;
      double time() const override { return exing_->time();}
      bool active() const override { return  exing_->active(); }
      void process(FitState& kkdata,TimeDir tdir) override;
      void updateState(MetaIterConfig const& miconfig,bool first) override;
      void updateConfig(Config const& config) override {}
      void append(PTRAJ& fit,TimeDir tdir) override;
      void updateReference(KTRAJPTR const& ltrajptr) override;
      Chisq chisq(Parameters const& pdata) const override { return Chisq();}
      void print(std::ostream& ost=std::cout,int detail=0) const override;
      virtual ~Material(){}
      // create from the material and a trajectory
      Material(EXINGPTR const& dxing, PTRAJ const& ptraj);
      // accessors
      auto const& cache() const { return cache_; }
      auto const& elementXing() const { return *exing_; }
      auto const& elementXingPtr() const { return exing_; }
      auto const& referenceTrajectory() const { return exing_->referenceTrajectory(); }
    private:
      EXINGPTR exing_; // element crossing for this effect
      Weights cache_; // cache of weight processing in opposite directions, used to build the fit trajectory
  };

  template<class KTRAJ> Material<KTRAJ>::Material(EXINGPTR const& dxing, PTRAJ const& ptraj) : exing_(dxing) {}

  template<class KTRAJ> void Material<KTRAJ>::process(FitState& kkdata,TimeDir tdir) {
    if(exing_->active()){
      // forwards, set the cache AFTER processing this effect
      if(tdir == TimeDir::forwards) {
        kkdata.append(exing_->parameters(tdir));
        cache_ += kkdata.wData();
      } else {
        // backwards, set the cache BEFORE processing this effect, to avoid double-counting it
        cache_ += kkdata.wData();
        kkdata.append(exing_->parameters(tdir));
      }
    }
  }

  template<class KTRAJ> void Material<KTRAJ>::updateState(MetaIterConfig const& miconfig,bool first) {
    // update the ElementXing
    exing_->updateState(miconfig,first);
    // reset the cached weights
    cache_ = Weights();
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
      newpiece.params() = Parameters(cache_);
      // make sure the range includes the transit time
      newpiece.range() = (tdir == TimeDir::forwards) ? TimeRange(etime,std::max(ptraj.range().end(),etime+exing_->transitTime())) :
        TimeRange(std::min(ptraj.range().begin(),etime-exing_->transitTime()),etime);
      if( tdir == TimeDir::forwards)
        ptraj.append(newpiece);
      else
        ptraj.prepend(newpiece);
    }
    // update the xing
    if( tdir == TimeDir::forwards)
      exing_->updateReference(ptraj.backPtr());
    else
      exing_->updateReference(ptraj.frontPtr());
  }

  template<class KTRAJ> void Material<KTRAJ>::updateReference(KTRAJPTR const& ltrajptr) {
    exing_->updateReference(ltrajptr);
  }

  template<class KTRAJ> void Material<KTRAJ>::print(std::ostream& ost,int detail) const {
    ost << "Material " << static_cast<Effect<KTRAJ>const&>(*this);
    ost << " ElementXing ";
    exing_->print(ost,detail);
    if(detail >3){
      ost << " cache ";
      cache().print(ost,detail);
    }
  }

  template <class KTRAJ> std::ostream& operator <<(std::ostream& ost, Material<KTRAJ> const& kkmat) {
    kkmat.print(ost,0);
    return ost;
  }
}
#endif
