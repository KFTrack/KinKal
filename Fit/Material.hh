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
      Parameters const& effect() const { return mateff_; }
      Weights const& cache() const { return cache_; }
      EXING const& elementXing() const { return *exing_; }
      KTRAJ const& referenceTrajectory() const { return exing_->referenceTrajectory(); }
    private:
      // update the local cache representing the effect of this material on the reference parameters
      void updateCache();
      EXINGPTR exing_; // element crossing for this effect
      Parameters mateff_; // parameter space description of this effect
      Weights cache_; // cache of weight processing in opposite directions, used to build the fit trajectory
      double vscale_; // variance factor due to annealing 'temperature'
  };

  template<class KTRAJ> Material<KTRAJ>::Material(EXINGPTR const& dxing, PTRAJ const& ptraj) : exing_(dxing),
  vscale_(1.0) {
  }

  template<class KTRAJ> void Material<KTRAJ>::process(FitState& kkdata,TimeDir tdir) {
    if(exing_->active()){
      // forwards, set the cache AFTER processing this effect
      if(tdir == TimeDir::forwards) {
        kkdata.append(mateff_,tdir);
        cache_ += kkdata.wData();
      } else {
        // backwards, set the cache BEFORE processing this effect, to avoid double-counting it
        cache_ += kkdata.wData();
        kkdata.append(mateff_,tdir);
      }
    }
  }

  template<class KTRAJ> void Material<KTRAJ>::updateState(MetaIterConfig const& miconfig,bool first) {
    if(first)vscale_ = miconfig.varianceScale();
    exing_->updateState(miconfig,first);
    updateCache();
  }

  template<class KTRAJ> void Material<KTRAJ>::updateCache() {
    // reset the weight
    cache_ = Weights();
    // reset parameters before rebuilding from scratch
    mateff_ = Parameters();
    if(exing_->active()){
      // loop over the momentum change basis directions, adding up the effects on parameters from each
      std::array<double,3> dmom = {0.0,0.0,0.0}, momvar = {0.0,0.0,0.0};
      exing_->materialEffects(TimeDir::forwards, dmom, momvar);
      // get the parameter derivative WRT momentum
      DPDV dPdM = referenceTrajectory().dPardM(time());
      double mommag = referenceTrajectory().momentum(time());
      for(int idir=0;idir<MomBasis::ndir; idir++) {
        auto mdir = static_cast<MomBasis::Direction>(idir);
        auto dir = referenceTrajectory().direction(time(),mdir);
        // project the momentum derivatives onto this direction
        DVEC pder = mommag*(dPdM*SVEC3(dir.X(), dir.Y(), dir.Z()));
        // convert derivative vector to a Nx1 matrix
        ROOT::Math::SMatrix<double,NParams(),1> dPdm;
        dPdm.Place_in_col(pder,0,0);
        // update the transport for this effect; first the parameters.  Note these are for forwards time propagation (ie energy loss)
        mateff_.parameters() += pder*dmom[idir];
        // now the variance: this doesn't depend on time direction
        ROOT::Math::SMatrix<double, 1,1, ROOT::Math::MatRepSym<double,1>> MVar;
        MVar(0,0) = momvar[idir]*vscale_;
        mateff_.covariance() += ROOT::Math::Similarity(dPdm,MVar);
      }
    }
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
    ost << " effect ";
    effect().print(ost,detail-2);
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
