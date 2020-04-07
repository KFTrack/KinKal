#ifndef KinKal_KKEnd_hh
#define KinKal_KKEnd_hh
//
// End of the KK fit effect set, used transiently to start the fit and trajectory building
//
#include "KinKal/KKWEff.hh"
#include <stdexcept>
#include <limits>
#include <ostream>

namespace KinKal {
  template<class KTRAJ> class KKEnd : public KKWEff<KTRAJ> {
    public:
      typedef KKEff<KTRAJ> KKEFF;
      typedef KKWEff<KTRAJ> KKWEFF;
      typedef PKTraj<KTRAJ> PKTRAJ;
      typedef typename KTRAJ::PDATA PDATA; // forward derivative type
      typedef typename KKEFF::WDATA WDATA; // forward the typedef
      typedef typename KKEFF::KKDATA KKDATA;

      // provide interface
      virtual bool update(PKTRAJ const& ref) override;
      virtual float time() const override { return (tdir_ == TDir::forwards) ? -std::numeric_limits<float>::max() : std::numeric_limits<float>::max(); } // make sure this is always at the end
      virtual unsigned nDOF() const override { return 0; }
      virtual bool isActive() const override { return true; }
      virtual float chisq(PDATA const& pars) const override { return 0.0; }
      virtual bool process(KKDATA& kkdata,TDir tdir) override;
      virtual bool append(PKTRAJ& fit) override;
      virtual void print(std::ostream& ost=std::cout,int detail=0) const override;
      virtual ~KKEnd(){}
      // accessors
      TDir const& tDir() const { return tdir_; }
      double deWeighting() const { return dwt_; }
      KTRAJ const& end() const { return end_; }

      // construct from trajectory and direction.  Deweighting must be tuned to balance stability vs bias
      KKEnd(PKTRAJ const& pktraj,TDir tdir, double dweight=1e6); 
    private:
      TDir tdir_; // direction for this effect; note the early end points forwards, the late backwards
      double dwt_; // deweighting factor
      KTRAJ end_; // cache of parameters at the end of processing this direction, used in traj creation
 };

  template <class KTRAJ> KKEnd<KTRAJ>::KKEnd(PKTRAJ const& pktraj, TDir tdir, double dweight) :
    tdir_(tdir) , dwt_(dweight), end_(tdir == TDir::forwards ? pktraj.front() : pktraj.back()){
      update(pktraj);
    }

  template<class KTRAJ> bool KKEnd<KTRAJ>::process(KKDATA& kkdata,TDir tdir) {
    bool retval(true);
    if(tdir == tdir_) 
    // start the fit with the de-weighted info from the previous iteration or seed
      retval = KKWEFF::process(kkdata,tdir);
    else {
    // at the opposite end, cache the final parameters
      end_.params() = kkdata.pData();
      KKEffBase::setStatus(tdir,KKEffBase::processed);
    }
    return retval;
  }

  template<class KTRAJ> bool KKEnd<KTRAJ>::update(PKTRAJ const& ref) {
    auto refend = ref.nearestPiece(time()).params();
    refend.covariance() *= dwt_;
    // convert this to a weight (inversion)
    KKWEFF::wdata_ = WData<PKTRAJ::NParams()>(refend);
    KKEffBase::updateStatus();
    return KKWEFF::wData().matrixOK();
  }

  template<class KTRAJ> bool KKEnd<KTRAJ>::append(PKTRAJ& fit) {
    // if the fit is empty and we're going in the right direction, take the end cache and
    // seed the fit with it
    if(tdir_ == TDir::forwards) {
      if(fit.pieces().size() == 0){
	// start with a very large range 
	end_.range() = TRange(-std::numeric_limits<float>::max(),std::numeric_limits<float>::max());
	// append this to the (empty) fit
	fit.append(end_);
      } else
	return false;	
    }
    return true;
  }

  template<class KTRAJ> void KKEnd<KTRAJ>::print(std::ostream& ost,int detail) const {
    ost << "KKEnd " << static_cast<KKEff<KTRAJ>const&>(*this) << " direction " << tDir() << " deweight " << deWeighting() << std::endl;
    end().print(ost,detail);
    if(detail > 0){
      ost << "End ";
      end().print(ost,detail);
      ost << "Weight " << this->wData() << std::endl;
    }


  }
  template <class KTRAJ> std::ostream& operator <<(std::ostream& ost, KKEnd<KTRAJ> const& kkend) {
    kkend.print(ost,0);
    return ost;
  }


}
#endif
