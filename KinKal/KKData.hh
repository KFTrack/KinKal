#ifndef KinKal_KKData_hh
#define KinKal_KKData_hh
//
//  Data payload for processing the fit.  This object exists in both
// parameter and weight space, with lazy evaluation to go between the
// two with the minimum of matrix inversions
//
#include "KinKal/WData.hh"
#include "KinKal/PData.hh"
namespace KinKal {
  class KKData {
    public:
      KKData() : hasPData_(false), hasWData_(false) {}
      KKData(PData const& pdata) : pdata_(pdata), hasPData_(true), hasWData_(false) {}
      KKData(WData const& wdata) : wdata_(wdata), hasPData_(false), hasWData_(true) {}
      // accessors
      bool hasPData() const { return hasPData_; }
      bool hasWData() const { return hasWData_; }
      // add to either parameters or weights
      void append(PData const& pdata) {
	pData() += pdata;
	// this invalidates the weight information
	hasPData_ = true;
	hasWData_ = false;
      }
      void append(WData const& wdata) {
	wData() += wdata;
	// this invalidates the parameter information
	hasWData_ = true;
	hasPData_ = false;
      }
      PData& pData() { 
	if(!hasPData_ && hasWData_ ){
	  // invert the weight
	  pdata_ = PData(wdata_);
	  hasPData_ = true;
	}
	return pdata_;
      }
      WData& wData() { 
	if(!hasWData_ && hasPData_ ){
	  // invert the parameters
	  wdata_ = WData(pdata_);
	  hasWData_ = true;
	}
	return wdata_;
      }
    private:
      PData pdata_; // parameters space representation of (intermediate) fit data
      WData wdata_; // weight space representation of fit data
      bool hasPData_, hasWData_;  // keep track of validity for lazy evaluation (cache coherence)
  };
}
#endif
