#ifndef KinKal_KKBFCorr_hh
#define KinKal_KKBFCorr_hh
//
// Class to correct for BField inhomogenity; adjust the parameters for the momentum change
// This effect adds no information content or noise (presently), just transports the parameters 
//
#include "KinKal/KKEff.hh"
#include "KinKal/TDir.hh"
#include <iostream>
#include <stdexcept>
#include <array>
#include <ostream>

namespace KinKal {
  template<class KTRAJ> class KKBFCorr : public KKEff<KTRAJ> {
    public:
      typedef KKEff<KTRAJ> KKEFF;
      typedef PKTraj<KTRAJ> PKTRAJ;
      typedef typename KKEFF::PDATA PDATA; // forward the typedef
      typedef typename KKEFF::WDATA WDATA; // forward the typedef
      typedef KKData<PDATA::PDim()> KKDATA;
      typedef typename KTRAJ::PDER PDER; // forward the typedef
      virtual float time() const override { return drange_.mid(); } // apply the correction at the middle of the range
      virtual bool isActive() const override { return active_;}
      virtual void update(PKTRAJ const& ref) override;
      virtual void update(PKTRAJ const& ref, MConfig const& mconfig) override;
      virtual void print(std::ostream& ost=std::cout,int detail=0) const override;
      virtual void process(KKDATA& kkdata,TDir tdir) override;
      virtual void append(PKTRAJ& fit) override;
      PDER const& effect() const { return bfeff_; }
      WDATA const& cache() const { return cache_; }
      virtual ~KKBFCorr(){}
      // create from the domain range, the effect, and the
      KKBFCorr(BField const& bfield, unsigned nsteps, PKTRAJ const& pktraj,TRange const& drange) : 
	bfield_(bfield), nsteps_(nsteps), ref_(pktraj.nearestPiece(drange.mid())), drange_(drange), active_(false) {} // not active until updated
    private:
      BField const& bfield_; // bfield
      unsigned nsteps_; // number of steps to integrate over this domain
      KTRAJ ref_; // reference to local trajectory
      TRange drange_; // extent of this domain
      PDATA bfeff_; // effect of the difference beween the actual BField and bnom integrated over this integral
      WDATA cache_; // cache of weight processing in opposite directions, used to build the fit trajectory
      bool active_; // activity state
  };

  template<class KTRAJ> void KKBFCorr<KTRAJ>::process(KKDATA& kkdata,TDir tdir) {
    if(this->isActive()){
      // forwards, set the cache AFTER processing this effect
      if(tdir == TDir::forwards) {
	kkdata.append(bfeff_);
	cache_ += kkdata.wData();
      } else {
      // backwards, set the cache BEFORE processing this effect, to avoid double-counting it
	cache_ += kkdata.wData();
	// SUBTRACT the effect going backwards: covariance change is sign-independent
	PDATA reverse(bfeff_);
	reverse.parameters() *= -1.0;
      	kkdata.append(reverse);
      }
    }
    KKEffBase::setStatus(tdir,KKEffBase::processed);
  }

  template<class KTRAJ> void KKBFCorr<KTRAJ>::update(PKTRAJ const& ref) {
    cache_ = WDATA();
    ref_ = ref.nearestPiece(drange_.mid()); 
    KKEffBase::updateStatus();
  }

  template<class KTRAJ> void KKBFCorr<KTRAJ>::update(PKTRAJ const& ref, MConfig const& mconfig) {
    if(mconfig.updatebfcorr_){
      active_ = true;
      // integrate the field difference over the new traj
      unsigned istep(0);
      float tstep = drange_.range()/(nsteps_-1);
      PDER dp; // change in parameters due to inhomogenous B integral
      float tsamp = drange_.low();
      while(istep < nsteps_){
	// find the local piece: this avoids search every call
	auto const& ktraj = ref.nearestPiece(tsamp);
	Vec3 tpos, bvec, bnom;
	bnom = ktraj.bnom();
	float bnommag = bnom.R();
	ktraj.position(tsamp,tpos);
	bfield_.fieldVect(bvec,tpos);
	// get the perp basis for the BField x-product at this point
	Vec3 t1hat, t2hat;
	ktraj.dirVector(KInter::theta1,tsamp,t1hat);
	ktraj.dirVector(KInter::theta2,tsamp,t2hat);
	// project the BField diff in these directions
	auto dbvec = bvec - bnom;
	float fdbt1 = dbvec.Dot(t1hat)/bnommag;
	float fdbt2 = dbvec.Dot(t2hat)/bnommag;
	// find the derivative on the local parameters
	PDER dpdt1, dpdt2;
	ktraj.momDeriv(KInter::theta1,tsamp,dpdt1);
	ktraj.momDeriv(KInter::theta2,tsamp,dpdt2);
	// add these to the drange_ calculation; this is the cross-product, for magnetic force
	dp += (tstep*CLHEP::c_light/ktraj.ebar())*(fdbt2*dpdt1 - fdbt1*dpdt2);
	istep++;
	tsamp += tstep;
      }
      bfeff_ = PDATA(dp); // no BF contribution to noise currently FIXME!
    }
    update(ref);
  }

  template<class KTRAJ> void KKBFCorr<KTRAJ>::append(PKTRAJ& fit) {
    if(isActive()){
      // create a trajectory piece from the cached weight
      float time = this->time();
      KTRAJ newpiece(ref_);
      newpiece.params() = PDATA(cache_);
      newpiece.range() = TRange(time,fit.range().high());
      // make sure the piece is appendable
      if(time > fit.back().range().low()){
	fit.append(newpiece);
      } else {
	throw std::invalid_argument("Can't append piece");
      }
    }
  }

  template<class KTRAJ> void KKBFCorr<KTRAJ>::print(std::ostream& ost,int detail) const {
    ost << "KKBFCorr " << static_cast<KKEff<KTRAJ>const&>(*this);
    ost << " effect " << bfeff_
    << " domain range " << drange_ << std::endl;
    if(detail >0){
      ost << " cache ";
      cache().print(ost,detail);
      ost << "Reference " << ref_ << std::endl;
    }
  }

  template <class KTRAJ> std::ostream& operator <<(std::ostream& ost, KKBFCorr<KTRAJ> const& kkmat) {
    kkmat.print(ost,0);
    return ost;
  }

}
#endif
