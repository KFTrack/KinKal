#ifndef KinKal_Weights_hh
#define KinKal_Weights_hh
//
//  Data object describing weight-space information
//  used as part of the kinematic kalman fit
//
#include "KinKal/Data.hh"
#include <ostream>
namespace KinKal {
  class Parameters;
  class Weights {
    public:
      // construct from vector and matrix
      Weights(DVEC const& wvec, DMAT const& wmat) : tdata_(wvec,wmat) {}
      Weights(DVEC const& wvec) : tdata_(wvec) {}
      Weights(Parameters const& pdata);
      Weights() {}
      // accessors; just re-interpret the base class accessors
      DVEC const& weightVec() const { return tdata_.vec(); }
      DMAT const& weightMat() const { return tdata_.mat(); }
      DVEC& weightVec() { return tdata_.vec(); }
      DMAT& weightMat() { return tdata_.mat(); }
      Data const& tData() const { return tdata_; }
      Data& tData() { return tdata_; }
      // addition: only works for other weights
      Weights & operator +=(Weights const& other) {
	tdata_ += other.tdata_;
	return *this;
      }
      Weights & operator -=(Weights const& other) {
	tdata_ -= other.tdata_;
	return *this;
      }
      void print(std::ostream& ost=std::cout,int detail=0) const {
	ost << "Weights wVec " << weightVec() << std::endl;
	if(detail > 1)
	  ost << "weight " << weightMat() << std::endl;
      }
    private:
      Data tdata_; // data payload
  };
  std::ostream& operator << (std::ostream& ost, Weights const& wdata);
}
#endif
