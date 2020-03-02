#ifndef KinKal_KKData_hh
#define KinKal_KKData_hh
//
//  Data payload for processing the fit.  This object exists in both
// parameter and weight space, with lazy evaluation to go between the
// two with the minimum of matrix inversions
//
#include "KinKal/PData.hh"
#include "KinKal/WData.hh"
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
      PDATA const& pData() const { return pdata_; } // note these just return the local cache, which can be invalid!
      WDATA const& wData() const { return wdata_; }
      // add to either parameters or weights
      void addPData(PDATA const& pdata) {
	pData() += pdata;
	// this invalidates the weight information
	hasWData_ = false;
      }
      void addWData(WDATA const& wdata) {
	wData() += wdata;
	// this invalidates the parameter information
	hasPData_ = false;
      }
      // non-const accessors; these perform auto-translation if necessary
      PDATA& pData()  { 
	if(!hasPData_){
	  // invert the weight
	  pdata_.invert(wdata_);
	  hasPData_ = pdata_.matrixOK();
	}
	return pdata_;
      }
      WDATA& wData()  { 
	if(!hasWData_){
	  // invert the weight
	  wdata_.invert(pdata_);
	  hasWData_ = wdata_.matrixOK();
	}
	return wdata_;
      }
    private:
      PDATA pdata_; // parameters space representation of (intermediate) fit data
      WDATA wdata_; // weight space representation of fit data
      bool hasPData_, hasWData_; 
  };
}
#endif
