#ifndef KinKal_KKWeight_hh
#define KinKal_KKWeight_hh
//
// Class to add information; a measurement or constraint
// This effect provides information content and is processed in weight space 
//
#include "KinKal/KKEffect.hh"
namespace KinKal {
  template<class KTRAJ> class KKWeight : public KKEffect<KTRAJ> {
    public:
    typedef typename KKEffect<KTRAJ>::WDATA WDATA; // forward the typedef
      // process this site given the adjacent site
      virtual bool process(KKEffect<KTRAJ> const& other,TDir tdir) override;
      // construct from a weight
      KKWeight(KTRAJ const& reftraj, WDATA const& weight) : KKEffect<KTRAJ>(reftraj), weight_(weight) {}
      virtual unsigned nDOF() const = 0;
      virtual ~KKWeight(){}
      WDATA const& weight() const { return weight_; }
    protected:
    // allow subclasses to construct without a weight
      KKWeight(KTRAJ const& reftraj) : KKEffect<KTRAJ>(reftraj) {}
      WDATA weight_; // weight representation of this site's constraint/measurement
  };

  template<class KTRAJ> bool KKWeight<KTRAJ>::process(KKEffect<KTRAJ> const& other,TDir tdir) {
    bool retval(false);
    if( other.status(tdir) == KKEffect<KTRAJ>::processed) {
      if(KKEffect<KTRAJ>::isActive()){
      // copy in the other sites weight to my cache
	KKEffect<KTRAJ>::weights_[tdir] = other.weights(tdir);
	if(KKEffect<KTRAJ>::weights_[tdir].matrixOK()){
	// add this site's information
	  KKEffect<KTRAJ>::weights_[tdir] += weight_;
	  KKEffect<KTRAJ>::status_[tdir] = KKEffect<KTRAJ>::processed;
	  retval = true;
	}
      } else {
	// copy over the raw state from the previous site
	KKEffect<KTRAJ>::weights_[tdir] = other.weights(tdir);
	KKEffect<KTRAJ>::params_[tdir] = other.params(tdir);
	retval = true;
      }
    }
    return retval;
  }

}
#endif
