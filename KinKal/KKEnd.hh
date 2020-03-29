#ifndef KinKal_KKEnd_hh
#define KinKal_KKEnd_hh
//
// End of the KK fit effect set, used transiently to start the fit and trajectory building
//
#include "KinKal/KKWEff.hh"
#include <stdexcept>
#include <limits>

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
      virtual double time() const override { return (tdir_ == TDir::forwards) ? -std::numeric_limits<double>::max() : std::numeric_limits<double>::max(); } // make sure this is always at the end
      virtual unsigned nDOF() const override { return 0; }
      virtual bool isActive() const override { return true; }
      virtual double chisq(PDATA const& pars) const override { return 0.0; }
      virtual bool process(KKDATA& kkdata,TDir tdir) override;
      virtual bool append(PKTRAJ& fit) override;
      // accessors
      TDir const& tDir() const { return tdir_; }
      double dWeighting() const { return dwt_; }
      KTRAJ const& end() const { return end_; }

      // construct from trajectory and direction.  Deweighting must be tuned to balance stability vs bias
      KKEnd(PKTRAJ const& pktraj,TDir tdir, double dweight=1e6); 
      virtual ~KKEnd(){}
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
    else
    // at the opposite end, cache the final parameters
      end_.params() = kkdata.pData();
    return retval;
  }

  template<class KTRAJ> bool KKEnd<KTRAJ>::update(PKTRAJ const& ref) {
    auto refend = ref.nearestPiece(time()).params();
    refend.covariance() *= dwt_;
    // convert this to a weight (inversion)
    KKWEFF::wdata_ = WData<PKTRAJ::NParams()>(refend,true);
    KKEffBase::updateStatus();
    return KKWEFF::wData().matrixOK();
  }

  template<class KTRAJ> bool KKEnd<KTRAJ>::append(PKTRAJ& fit) {
    // if the fit is empty and we're going in the right direction, take the end cache and
    // seed the fit with it
    if(tdir_ == TDir::forwards && fit.pieces().size() == 0){
    // start with the reference traj, and override the range and parameters
      end_.range() = fit.range();
      // append this to the (empty) fit
      bool ok = fit.append(end_);
      if(!ok)throw std::invalid_argument("append failed");
    }
    return true;
  }

}
#endif
