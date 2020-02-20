#ifndef KinKal_KKWeight_hh
#define KinKal_KKWeight_hh
//
// Class to add information; a measurement or constraint
// This effect provides information content and is processed in weight space 
//
#include "BTrk/KinKal/KKEffect.hh"
namespace KinKal {
  template<class PTRAJ> class KKWeight : public KKEffect<PTRAJ> {
    public:
      typedef typename KKEffect::WDATA WDATA; // forward the typedef
      // process this site given the adjacent site
      virtual bool process(KKEffect const& other,TDir tdir) override;
      // construct from a weight
      KKWeight(PTRAJ const& reftraj, WDATA const& weight) : KKEffect<PTRAJ>(reftraj), weight_(weight) {}
      virtual unsigned nDOF() const = 0;
    protected:
      KKWeight(PTRAJ const& reftraj) : KKEffect<PTRAJ>(reftraj) {}
      WDATA weight_; // weight representation of this site's constraint/measurement
  };

  template<> bool KKWeight<PTRAJ>::process(KKEffect const& other,TDir tdir) {
    bool retval(false);
    if( other.status() == KKEffect<PTRAJ>::processed) {
      if(isActive()){
      // copy in the other sites weight to my cache
	weight_[tdir] = other.weight(tdir);
	if(weight_[tdir].matrixOK()){
	// add this site's information
	  weight_[tdir] += weight_;
	  setStatus(tdir,TData::valid);
	  retval = true;
	}
      } else {
	// copy over the raw state from the previous site
	weight_[tdir] = other.weight_[tdir];
	params_[tdir] = other.params_[tdir];
	retval = true;
      }
    }
    return retval;
  }

}
#endif
