#ifndef KinKal_FitState_hh
#define KinKal_FitState_hh
//
//  Data payload for processing the fit.  This object exists in both
// parameter and weight space, with lazy evaluation to go between the
// two with the minimum of matrix inversions
//
#include "KinKal/General/Weights.hh"
#include "KinKal/General/Parameters.hh"
namespace KinKal {
  class FitState {
    public:
      FitState() : hasParameters_(false), hasWeights_(false) {}
      FitState(Parameters const& pdata) : pdata_(pdata), hasParameters_(true), hasWeights_(false) {}
      FitState(Weights const& wdata) : wdata_(wdata), hasParameters_(false), hasWeights_(true) {}
      // accessors
      bool hasParameters() const { return hasParameters_; }
      bool hasWeights() const { return hasWeights_; }
      // add to either parameters or weights
      void append(Parameters const& pdata) {
	pData() += pdata;
	// this invalidates the weight information
	hasParameters_ = true;
	hasWeights_ = false;
      }
      void append(Weights const& wdata) {
	wData() += wdata;
	// this invalidates the parameter information
	hasWeights_ = true;
	hasParameters_ = false;
      }
      Parameters& pData() { 
	if(!hasParameters_ && hasWeights_ ){
	  // invert the weight
	  pdata_ = Parameters(wdata_);
	  hasParameters_ = true;
	}
	return pdata_;
      }
      Weights& wData() { 
	if(!hasWeights_ && hasParameters_ ){
	  // invert the parameters
	  wdata_ = Weights(pdata_);
	  hasWeights_ = true;
	}
	return wdata_;
      }
    private:
      Parameters pdata_; // parameters space representation of (intermediate) fit data
      Weights wdata_; // weight space representation of fit data
      bool hasParameters_, hasWeights_;  // keep track of validity for lazy evaluation (cache coherence)
  };
}
#endif
