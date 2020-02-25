#ifndef KinKal_KKWeight_hh
#define KinKal_KKWeight_hh
//
// Class to add information; a measurement or constraint
// This effect provides information content and is processed in weight space 
//
#include "KinKal/KKEff.hh"
namespace KinKal {
  template<class KTRAJ> class KKWeight : public KKEff<KTRAJ> {
    public:
      typedef KKEff<KTRAJ> KKEFF;
      typedef typename KKEFF::WDATA WDATA; // forward the typedef
      // process this effect given the adjacent effect
      virtual bool process(KKEFF const& other,TDir tdir) override;
      // construct from a weight
      KKWeight(KTRAJ const& reftraj, WDATA const& weight) : KKEFF(reftraj), weight_(weight) {}
      virtual ~KKWeight(){}
      WDATA const& weight() const { return weight_; }
    protected:
    // allow subclasses to construct without a weight
      KKWeight(KTRAJ const& reftraj) : KKEFF(reftraj) {}
      KKWeight() {} // no local reftraj; must be set
      WDATA weight_; // weight representation of this effect's constraint/measurement
  };

  template<class KTRAJ> bool KKWeight<KTRAJ>::process(KKEFF const& other,TDir tdir) {
    bool retval(false);
    auto idir = static_cast<std::underlying_type<TDir>::type>(tdir);
    if( other.status(tdir) == KKEFF::processed) {
      if(this->isActive()){
      // copy in the other effects weight to my cache
	KKEFF::weights_[idir] = other.weights(tdir);
	if(KKEFF::weights_[idir].matrixOK()){
	// add this effect's information
	  KKEFF::weights_[idir] += weight_;
	  retval = true;
	}
      } else {
	// copy over the raw state from the previous effect
	KKEFF::weights_[idir] = other.weights(tdir);
	KKEFF::params_[idir] = other.params(tdir);
	retval = true;
      }
      KKEFF::setStatus(tdir,KKEFF::processed);
      return retval;
    }
    if(retval)KKEffBase::setStatus(tdir,KKEffBase::processed);
    return retval;
  }

}
#endif
