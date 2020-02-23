#ifndef KinKal_KKEnd_hh
#define KinKal_KKEnd_hh
//
// End of the KK fit site set, used transiently to start the fit 
//
#include "KinKal/KKWeight.hh"
#include "KinKal/PKTraj.hh"
#include "KinKal/PData.hh"
#include "KinKal/WData.hh"
#include <stdexcept>
#include <limits>

namespace KinKal {
  template<class KTRAJ> class KKEnd : public KKWeight<KTRAJ> {
    public:
      typedef KKEff<KTRAJ> KKEFF;
      typedef KKWeight<KTRAJ> KKWT;
      typedef typename KKEFF::WDATA WDATA; // forward the typedef
      typedef PKTraj<KTRAJ> PKTRAJ;

      // provide interface
      virtual bool process(KKEFF const& other,TDir tdir) override;
      virtual bool update(PKTRAJ const& ref) override;
      virtual double time() const override { return (tdir_ == TDir::forwards) ? -std::numeric_limits<double>::max() : std::numeric_limits<double>::max(); } // make sure this is always at the end
      virtual unsigned nDOF() const override { return 0; }

      // construct from trajectory and direction.  Deweighting must be tuned to balance stability vs bias
      KKEnd(PKTRAJ const& pktraj,TDir tdir, double dweight=1e6); 
      virtual ~KKEnd(){}
    private:
      TDir tdir_; // direction for this site
      double dwt_; // deweighting factor
  };

  template <class KTRAJ> bool KKEnd<KTRAJ>::process(KKEFF const& other,TDir tdir) {
    bool retval(false);
    if(tdir == tdir_){
      throw std::runtime_error("Cannot process KKEnd in its own direction");
    } else {
      // copy over previous state
      auto idir = static_cast<std::underlying_type<TDir>::type>(tdir);
      KKEFF::weights_[idir] = other.weights(tdir);
      KKEFF::params_[idir] = other.params(tdir);
      retval = true;
    }
    return retval;
  }

  template <class KTRAJ> KKEnd<KTRAJ>::KKEnd(PKTRAJ const& pktraj, TDir tdir, double dweight) : 
    KKWT(pktraj.front()), tdir_(tdir) , dwt_(dweight) {
      auto idir = static_cast<std::underlying_type<TDir>::type>(tdir);
      update(pktraj);
      // auto-process this site by copying the state
      KKEFF::weights_[idir] = KKWT::weight_;
      KKEFF::setStatus(tdir, KKEFF::processed); 
    }

  template<class KTRAJ> bool KKEnd<KTRAJ>::update(PKTRAJ const& ref) {
    KKEFF::resetRefTraj(ref);
    // extract the parameter vector and deweight the covariance
    auto endpars = KKEFF::referenceTraj().params();
    endpars.covariance() *= dwt_;
    // convert this to a weight (inversion)
    KKWT::weight_ = WData<PKTRAJ::NParams()>(endpars,true);
    return KKWT::weight_.matrixOK();
  }

}
#endif
