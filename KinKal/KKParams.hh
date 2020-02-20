#ifndef KinKal_KKParams_hh
#define KinKal_KKParams_hh
//
// Class to describe parameter transport/noise i
// This effect provides information content and is processed in params space 
//
#include "KinKal/KKEffect.hh"
namespace KinKal {
  template<class KTRAJ> class KKParams : public KKEffect<KTRAJ> {
    public:
      typedef typename KKEffect::WDATA WDATA; // forward the typedef
      // process this site given the adjacent site
      virtual bool process(KKEffect const& other,TDir tdir) override;
      // construct from a params
      KKParams(WDATA const& params) : params_(params) {}
    protected:
      KKParams() {}
      WDATA params_; // params of this site
  };

  template<> bool KKParams<KTRAJ>::process(KKEffect const& other,TDir tdir) {
    bool retval(false);
    if( other.status() == KKEffect<KTRAJ>::processed) {
      if(isActive()){
      // copy in the other sites params to my cache
	params_[tdir] = other.params(tdir);
	if(params_[tdir].matrixOK()){
	// add this site's information
	  params_[tdir] += params_;
	  setStatus(tdir,TData::valid);
	  retval = true;
	}
      } else {
	// simply copy over the statefrom the previous site
	params_[tdir] = other.params(tdir);
	params_[tdir] = other.params(tdir);
	retval = true;
      }
    }
    return retval;
  }

}
#endif
