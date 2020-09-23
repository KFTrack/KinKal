#ifndef KinKal_Weights_hh
#define KinKal_Weights_hh
//
//  Data object describing weight-space information
//  used as part of the kinematic kalman fit
//
#include "KinKal/Fit/FitData.hh"
#include <ostream>
namespace KinKal {
  class Parameters;
  class Weights {
    public:
      // construct from vector and matrix
      Weights(DVEC const& wvec, DMAT const& wmat) : fitdata_(wvec,wmat) {}
      Weights(DVEC const& wvec) : fitdata_(wvec) {}
      Weights(Parameters const& pdata);
      Weights() {}
      // accessors; just re-interpret the base class accessors
      DVEC const& weightVec() const { return fitdata_.vec(); }
      DMAT const& weightMat() const { return fitdata_.mat(); }
      DVEC& weightVec() { return fitdata_.vec(); }
      DMAT& weightMat() { return fitdata_.mat(); }
      FitData const& fitData() const { return fitdata_; }
      FitData& fitData() { return fitdata_; }
      // addition: only works for other weights
      Weights & operator +=(Weights const& other) {
	fitdata_ += other.fitdata_;
	return *this;
      }
      Weights & operator -=(Weights const& other) {
	fitdata_ -= other.fitdata_;
	return *this;
      }
      Weights & operator *=(double scale){
	fitdata_.vec() *= scale;
	fitdata_.mat() *= scale;
	return *this;
      }
      void print(std::ostream& ost=std::cout,int detail=0) const {
	ost << "Weights wVec " << weightVec() << std::endl;
	if(detail > 1)
	  ost << "weight " << weightMat() << std::endl;
      }
    private:
      FitData fitdata_; // data payload
  };
  std::ostream& operator << (std::ostream& ost, Weights const& wdata);
}
#endif
