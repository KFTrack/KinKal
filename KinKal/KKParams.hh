#ifndef KinKal_KKParams_hh
#define KinKal_KKParams_hh
//
// Class to describe parameter transport/noise i
// This effect provides information content and is processed in params space 
//
#include "KinKal/KKEff.hh"
namespace KinKal {
  template<class KTRAJ> class KKParams : public KKEff<KTRAJ> {
    public:
      typedef typename KKEff::PDATA PDATA; // forward the typedef
      // process this effect given the adjacent effect
      virtual bool process(KKEff const& other,TDir tdir) override;
      virtual unsigned nDOF() const { return 0; } // no information content, therefore no NDOF
      virtual double chisq(PDATA const& pars) const { return 0.0; }
      // construct from a params
      KKParams(PDATA const& params) : params_(params) {}
      virtual ~KKParams(){}
    protected:
      KKParams() {}
      PDATA params_; // params of this effect
  };

  template<> bool KKParams<KTRAJ>::process(KKEff const& other,TDir tdir) {
    bool retval(false);
    if( other.status() == KKEff<KTRAJ>::processed) {
      if(isActive()){
      // copy in the other effects params to my cache
	params_[tdir] = other.params(tdir);
	if(params_[tdir].matrixOK()){
	// add this effect's information
	  params_[tdir] += params_;
	  setStatus(tdir,TData::valid);
	  retval = true;
	}
      } else {
	// simply copy over the statefrom the previous effect
	params_[tdir] = other.params(tdir);
	params_[tdir] = other.params(tdir);
	retval = true;
      }
    }
    return retval;
  }

}
#endif
