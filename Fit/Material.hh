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
      using PKTRAJ = ParticleTrajectory<KTRAJ>;
      using EXING = ElementXing<KTRAJ>;
      using EXINGPTR = std::shared_ptr<EXING>;
      double time() const override { return dxing_->crossingTime() + tbuff_ ;} 
      bool active() const override { return  dxing_->active(); }
      void update(PKTRAJ const& ref) override;
      void update(PKTRAJ const& ref, MetaIterConfig const& miconfig) override;
      void print(std::ostream& ost=std::cout,int detail=0) const override;
      void process(FitState& kkdata,TimeDir tdir) override;
      void append(PKTRAJ& fit) override;
      virtual ~Material(){}
      // create from the material and a trajectory 
      Material(EXINGPTR const& dxing, PKTRAJ const& pktraj);
      // accessors
      Parameters const& effect() const { return mateff_; }
      Weights const& cache() const { return cache_; }
      EXING const& detXing() const { return *dxing_; }
      KTRAJ const& refKTraj() const { return ref_; }
    private:
      // update the local cache
      void updateCache();
      EXINGPTR dxing_; // detector piece crossing for this effect
      KTRAJ ref_; // reference to local trajectory
      Parameters mateff_; // parameter space description of this effect
      Weights cache_; // cache of weight processing in opposite directions, used to build the fit trajectory
      double vscale_; // variance factor due to annealing 'temperature'
      static double tbuff_; // small time buffer to avoid ambiguity
  };

  template<class KTRAJ> double Material<KTRAJ>::tbuff_ = 1.0e-3;

  template<class KTRAJ> Material<KTRAJ>::Material(EXINGPTR const& dxing, PKTRAJ const& pktraj) : dxing_(dxing), 
  ref_(pktraj.nearestPiece(dxing->crossingTime())), vscale_(1.0) {
    update(pktraj);
  }

  template<class KTRAJ> void Material<KTRAJ>::process(FitState& kkdata,TimeDir tdir) {
    if(dxing_->active()){
      // forwards, set the cache AFTER processing this effect
      if(tdir == TimeDir::forwards) {
	kkdata.append(mateff_);
	cache_ += kkdata.wData();
      } else {
	// backwards, set the cache BEFORE processing this effect, to avoid double-counting it
	cache_ += kkdata.wData();
	// SUBTRACT the effect going backwards: covariance change is sign-independent
	Parameters reverse(mateff_);
	reverse.parameters() *= -1.0;
	kkdata.append(reverse);
      }
    }
    KKEFF::setState(tdir,KKEFF::processed);
  }

  template<class KTRAJ> void Material<KTRAJ>::update(PKTRAJ const& ref) {
    cache_ = Weights();
    ref_ = ref.nearestPiece(dxing_->crossingTime()); 
    updateCache();
    KKEFF::updateState();
  }

  template<class KTRAJ> void Material<KTRAJ>::update(PKTRAJ const& ref, MetaIterConfig const& miconfig) {
    vscale_ = miconfig.varianceScale();
      // update the detector Xings for this effect
    dxing_->update(ref,miconfig);
    update(ref);
  }

  template<class KTRAJ> void Material<KTRAJ>::updateCache() {
    mateff_ = Parameters();
    if(dxing_->active()){
      // loop over the momentum change basis directions, adding up the effects on parameters from each
      std::array<double,3> dmom = {0.0,0.0,0.0}, momvar = {0.0,0.0,0.0};
      dxing_->materialEffects(ref_,TimeDir::forwards, dmom, momvar);
      // get the parameter derivative WRT momentum
      DPDV dPdM = ref_.dPardM(time());
      double mommag = ref_.momentum(time());
      for(int idir=0;idir<MomBasis::ndir; idir++) {
	auto mdir = static_cast<MomBasis::Direction>(idir);
	auto dir = ref_.direction(time(),mdir);
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

  template<class KTRAJ> void Material<KTRAJ>::append(PKTRAJ& fit) {
    if(dxing_->active()){
      // create a trajectory piece from the cached weight
      double time = this->time();
      KTRAJ newpiece(ref_);
      newpiece.params() = Parameters(cache_);
      // extend as necessary: absolute time can shift during iterations
      newpiece.range() = TimeRange(time,std::max(time+tbuff_,fit.range().end()));
      // make sure the piece is appendable; if not, adjust
      if(time < fit.back().range().begin()){
      // if this is the first piece, simply extend it back
	if(fit.pieces().size() ==1){
	  fit.front().setRange(TimeRange(newpiece.range().begin()-tbuff_,fit.range().end()));
	} else {
	  newpiece.setRange(TimeRange(fit.back().range().begin()+tbuff_,fit.range().end()));
	}
      }
      fit.append(newpiece);
    }
  }

  template<class KTRAJ> void Material<KTRAJ>::print(std::ostream& ost,int detail) const {
    ost << "Material " << static_cast<Effect<KTRAJ>const&>(*this);
    ost << " effect ";
    effect().print(ost,detail-2);
    ost << " ElementXing ";
    dxing_->print(ost,detail);
    if(detail >3){
      ost << " cache ";
      cache().print(ost,detail);
      ost << "Reference " << ref_ << std::endl;
    }
  }

  template <class KTRAJ> std::ostream& operator <<(std::ostream& ost, Material<KTRAJ> const& kkmat) {
    kkmat.print(ost,0);
    return ost;
  }

}
#endif
