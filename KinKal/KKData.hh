#ifndef KinKal_KKData_hh
#define KinKal_KKData_hh
//
//  Data payload for processing the fit.  This object exists in both
// parameter and weight space, with lazy evaluation to go between the
// two with the minimum of matrix inversions
//
#include "KinKal/WData.hh"
#include "KinKal/PData.hh"
#include <array>
namespace KinKal {
  template <size_t DDIM> class KKData {
    public:
      typedef PData<DDIM> PDATA; // forward the type declarations
      typedef WData<DDIM> WDATA; 
      KKData() : hasPData_(false), hasWData_(false) {}
      KKData(PDATA const& pdata) : pdata_(pdata), hasPData_(true), hasWData_(false) {}
      KKData(WDATA const& wdata) : wdata_(wdata), hasPData_(false), hasWData_(true) {}
      // accessors
      bool hasPData() const { return hasPData_; }
      bool hasWData() const { return hasWData_; }
      // add to either parameters or weights
      void append(PDATA const& pdata) {
	pData() += pdata;
	// this invalidates the weight information
	hasPData_ = true;
	hasWData_ = false;
      }
      void append(WDATA const& wdata) {
	wData() += wdata;
	// this invalidates the parameter information
	hasWData_ = true;
	hasPData_ = false;
      }
      PDATA& pData() { 
	if(!hasPData_ && hasWData_ ){
	  // invert the weight
	  pdata_ = PDATA(wdata_,true);
	  hasPData_ = pdata_.matrixOK();
	}
	return pdata_;
      }
      WDATA& wData() { 
	if(!hasWData_ && hasPData_ ){
	  // invert the parameters
	  wdata_ = WDATA(pdata_,true);
	  hasWData_ = wdata_.matrixOK();
	}
	return wdata_;
      }
    private:
      PDATA pdata_; // parameters space representation of (intermediate) fit data
      WDATA wdata_; // weight space representation of fit data
      bool hasPData_, hasWData_;  // keep track of validity for lazy evaluation
  };
}
#endif
