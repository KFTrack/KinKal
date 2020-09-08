#ifndef KinKal_KKMat_hh
#define KinKal_KKMat_hh
//
// Class to describe effect of a particle passing through discrete material on the fit (ie material transport)
// This effect adds no information content, just noise, and is processed in params space 
//
#include "KinKal/KKEff.hh"
#include "KinKal/DXing.hh"
#include "KinKal/TPoca.hh"
#include "KinKal/TDir.hh"
#include <iostream>
#include <stdexcept>
#include <array>
#include <ostream>

namespace KinKal {
  template<class KTRAJ> class KKMat : public KKEff<KTRAJ> {
    public:
      using KKEFF = KKEff<KTRAJ>;
      using PKTRAJ = PKTraj<KTRAJ>;
      using DXING = DXing<KTRAJ>;
      using DXINGPTR = std::shared_ptr<DXING>;
      virtual double time() const override { return dxing_->crossingTime() + 1.0e-3;} // small positive offset to disambiguate WRT hits should be a parameter FIXME!
      virtual bool isActive() const override { return active_ && dxing_->matXings().size() > 0; }
      virtual void update(PKTRAJ const& ref) override;
      virtual void update(PKTRAJ const& ref, MIConfig const& miconfig) override;
      virtual void print(std::ostream& ost=std::cout,int detail=0) const override;
      virtual void process(KKData& kkdata,TDir tdir) override;
      virtual void append(PKTRAJ& fit) override;
      void setTime(double time) { dxing_->crossingTime() = time; } // allow KKMHit to set the time (from the KKHit)
      virtual ~KKMat(){}
      // create from the material and a trajectory 
      KKMat(DXINGPTR const& dxing, PKTRAJ const& pktraj, bool active = true);
      // accessors
      PData const& effect() const { return mateff_; }
      WData const& cache() const { return cache_; }
      DXING const& detXing() const { return *dxing_; }
      KTRAJ const& refKTraj() const { return ref_; }
    private:
      // update the local cache
      void updateCache();
      DXINGPTR dxing_; // detector piece crossing for this effect
      KTRAJ ref_; // reference to local trajectory
      PData mateff_; // parameter space description of this effect
      WData cache_; // cache of weight processing in opposite directions, used to build the fit trajectory
      double vscale_; // variance factor due to annealing 'temperature'
      bool active_;
  };

   template<class KTRAJ> KKMat<KTRAJ>::KKMat(DXINGPTR const& dxing, PKTRAJ const& pktraj, bool active) : dxing_(dxing), 
   ref_(pktraj.nearestPiece(dxing->crossingTime())), vscale_(1.0), active_(active) {
     update(pktraj);
   }

  template<class KTRAJ> void KKMat<KTRAJ>::process(KKData& kkdata,TDir tdir) {
    if(active_){
      // forwards, set the cache AFTER processing this effect
      if(tdir == TDir::forwards) {
	kkdata.append(mateff_);
	cache_ += kkdata.wData();
      } else {
      // backwards, set the cache BEFORE processing this effect, to avoid double-counting it
	cache_ += kkdata.wData();
	// SUBTRACT the effect going backwards: covariance change is sign-independent
	PData reverse(mateff_);
	reverse.parameters() *= -1.0;
      	kkdata.append(reverse);
      }
    }
    KKEffBase::setStatus(tdir,KKEffBase::processed);
  }

  template<class KTRAJ> void KKMat<KTRAJ>::update(PKTRAJ const& ref) {
    cache_ = WData();
    ref_ = ref.nearestPiece(dxing_->crossingTime()); 
    updateCache();
    KKEffBase::updateStatus();
  }

  template<class KTRAJ> void KKMat<KTRAJ>::update(PKTRAJ const& ref, MIConfig const& miconfig) {
    vscale_ = miconfig.varianceScale();
    if(miconfig.updatemat_){
      // update the detector Xings for this effect
      dxing_->update(ref);
      // should check to see if this material is still active FIXME!
      update(ref);
    }
  }

  template<class KTRAJ> void KKMat<KTRAJ>::updateCache() {
    mateff_ = PData();
    if(dxing_->matXings().size() > 0){
      // loop over the momentum change basis directions, adding up the effects on parameters from each
      std::array<double,3> dmom = {0.0,0.0,0.0}, momvar = {0.0,0.0,0.0};
      dxing_->materialEffects(ref_,TDir::forwards, dmom, momvar);
      for(int idir=0;idir<MomBasis::ndir; idir++) {
	auto mdir = static_cast<MomBasis::Direction>(idir);
	// get the derivatives of the parameters WRT material effects
	// should call dPardM directly once and then project FIXME!
	DVEC pder = ref_.momDeriv(time(), mdir);
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

  template<class KTRAJ> void KKMat<KTRAJ>::append(PKTRAJ& fit) {
    if(active_){
      // create a trajectory piece from the cached weight
      double time = this->time();
      KTRAJ newpiece(ref_);
      newpiece.params() = PData(cache_);
      newpiece.range() = TRange(time,fit.range().end());
      // make sure the piece is appendable
      if(time > fit.back().range().begin()){
	fit.append(newpiece);
      } else {
	throw std::invalid_argument("KKMat: Can't append piece");
      }
    }
  }

  template<class KTRAJ> void KKMat<KTRAJ>::print(std::ostream& ost,int detail) const {
    ost << "KKMat " << static_cast<KKEff<KTRAJ>const&>(*this);
    ost << " effect ";
    effect().print(ost,detail-2);
    ost << " DXing ";
    dxing_->print(ost,detail);
    if(detail >3){
      ost << " cache ";
      cache().print(ost,detail);
      ost << "Reference " << ref_ << std::endl;
    }
  }

  template <class KTRAJ> std::ostream& operator <<(std::ostream& ost, KKMat<KTRAJ> const& kkmat) {
    kkmat.print(ost,0);
    return ost;
  }

}
#endif
