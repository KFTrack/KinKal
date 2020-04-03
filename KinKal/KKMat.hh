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
      typedef KKEff<KTRAJ> KKEFF;
      typedef PKTraj<KTRAJ> PKTRAJ;
      typedef DXing<KTRAJ> DXING;
      typedef typename KKEFF::PDATA PDATA; // forward the typedef
      typedef typename KKEFF::WDATA WDATA; // forward the typedef
      typedef KKData<PDATA::PDim()> KKDATA;
      typedef typename KTRAJ::PDER PDER; // forward the typedef
      virtual double time() const override { return dxing_.crossingTime() + 1.0e-3;} // small positive offset to disambiguate WRT hits should be a parameter FIXME!
      virtual bool isActive() const override { return active_; }
      virtual unsigned nDOF() const override { return 0; } 
      virtual double chisq(PDATA const& pars) const override { return 0.0; }
      virtual bool update(PKTRAJ const& ref) override;
      bool process(KKDATA& kkdata,TDir tdir) override;
      bool append(PKTRAJ& fit) override;
      // update for materials associated with a hit
      bool update(PKTRAJ const& ref, TPocaBase const& tpoca);
      PDATA const& effect() const { return pdata_; }
      WDATA const& cache() const { return wdata_; }
      virtual ~KKMat(){}
    // create from material and POCA
      KKMat(DXING& dxing, PKTRAJ const& pktraj, TPocaBase const& tdpoca, bool active = true);
      // create from just the material and a trajectory 
      KKMat(DXING& dxing, PKTRAJ const& pktraj, bool active = true); 
    private:
      // update the local cache
      void updateCache();
      // reset the cache. this must be done for each update cycle
      void resetCache() { wdata_ = WDATA(); pdata_ = PDATA();}
      DXING& dxing_; // detector piece crossing for this effect
      KTRAJ ref_; // reference to local trajectory
      PDATA pdata_; // parameter space description of this effect
      WDATA wdata_; // cache of weight processing in opposite directions, used to build the fit trajectory
      bool active_;
  };

   template<class KTRAJ> KKMat<KTRAJ>::KKMat(DXING& dxing, PKTRAJ const& pktraj, TPocaBase const& tdpoca, bool active) : dxing_(dxing),
   ref_(pktraj.nearestPiece(dxing.crossingTime())) , active_(active) {
     update(pktraj,tdpoca);
   }
   
   template<class KTRAJ> KKMat<KTRAJ>::KKMat(DXING& dxing, PKTRAJ const& pktraj, bool active) : dxing_(dxing), 
   ref_(pktraj.nearestPiece(dxing.crossingTime())), active_(active) {
     update(pktraj);
   }


  template<class KTRAJ> bool KKMat<KTRAJ>::process(KKDATA& kkdata,TDir tdir) {
    bool retval(false);
    if(this->isActive()){
      // forwards, set the cache AFTER processing this effect
      if(tdir == TDir::forwards) {
	kkdata.append(pdata_);
	wdata_ += kkdata.wData();
      } else {
      // backwards, set the cache BEFORE processing this effect, to avoid double-counting it
	wdata_ += kkdata.wData();
	// SUBTRACT the effect going backwards: covariance change is sign-independent
	PDATA reverse(pdata_);
	reverse.parameters() *= -1.0;
      	kkdata.append(reverse);
      }
      retval = kkdata.pData().matrixOK();
    }
    KKEffBase::setStatus(tdir,KKEffBase::processed);
    return retval;
  }

  template<class KTRAJ> bool KKMat<KTRAJ>::update(PKTRAJ const& ref) {
    dxing_.update(ref);
    ref_ = ref.nearestPiece(dxing_.crossingTime()); 
    updateCache();
    return true;
  }

  template<class KTRAJ> bool KKMat<KTRAJ>::update(PKTRAJ const& ref, TPocaBase const& tdpoca)  {
    dxing_.update(tdpoca);
    ref_ = ref.nearestPiece(dxing_.crossingTime()); 
    updateCache();
    return true;
  }

  template<class KTRAJ> void KKMat<KTRAJ>::updateCache() {
    // update the reference trajectory
    KKEffBase::updateStatus();
    resetCache();
    if(dxing_.matXings().size() > 0){
      // loop over the momentum change basis directions, adding up the effects on parameters from each
      std::array<double,3> dmom = {0.0,0.0,0.0}, momvar = {0.0,0.0,0.0};
      dxing_.momEffects(ref_,TDir::forwards, dmom, momvar);
      for(int idir=0;idir<=KInter::theta2; idir++) {
	auto mdir = static_cast<KInter::MDir>(idir);
	// get the derivatives of the parameters WRT material effects
	PDER pder;
	ref_.momDeriv(mdir, time(), pder);
	// convert derivative vector to a Nx1 matrix
	ROOT::Math::SMatrix<double,KTRAJ::NParams(),1> dPdm;
	dPdm.Place_in_col(pder,0,0);
	// update the transport for this effect; first the parameters.  Note these are for forwards time propagation (ie energy loss)
	this->pdata_.parameters() += pder*dmom[idir];
	// now the variance: this doesn't depend on time direction
	ROOT::Math::SMatrix<double, 1,1, ROOT::Math::MatRepSym<double,1> > MVar;
	MVar(0,0) = momvar[idir];
	this->pdata_.covariance() += ROOT::Math::Similarity(dPdm,MVar);
      }
    } else {
      // no material crossings: deactive this effect
      active_ = false;
    }
  }

  template<class KTRAJ> bool KKMat<KTRAJ>::append(PKTRAJ& fit) {
    // create a trajectory piece from the cached weight
    double time = this->time();
    KTRAJ newpiece(ref_);
    newpiece.params() = PDATA(wdata_,true);
    newpiece.range() = TRange(time,fit.range().high());
    // make sure the piece is appendable
    if(time > fit.back().range().low()){
      fit.append(newpiece);
    } else {
      return false;
    }
    return true;
  }

  template <class KTRAJ> std::ostream& operator <<(std::ostream& ost, KKMat<KTRAJ> const& kkmat) {
    ost << "KKMat " << static_cast<KKEff<KTRAJ>const&>(kkmat) << " effect " << kkmat.effect().parameters();
    return ost;
  }

}
#endif
