#ifndef KinKal_Parameters_hh
#define KinKal_Parameters_hh
//
//  Data object describing fit parameters
//  used as part of the kinematic kalman fit
//
#include "KinKal/Data.hh"
#include <ostream>
namespace KinKal {
  class Weights;
  class Parameters {
    public:
      // construct from vector and matrix
      Parameters(DVEC const& pars, DMAT const& pcov) : tdata_(pars,pcov) {}
      Parameters(DVEC const& pars) : tdata_(pars) {}
      Parameters(Weights const& wdata);
      Parameters() {}
      // accessors; just re-interpret the base class accessors
      DVEC const& parameters() const { return tdata_.vec(); }
      DMAT const& covariance() const { return tdata_.mat(); }
      DVEC& parameters() { return tdata_.vec(); }
      DMAT& covariance() { return tdata_.mat(); }
      Data const& tData() const { return tdata_; }
      Data& tData() { return tdata_; }
      // scale the matrix
      void scale(double sfac) { tdata_.scale(sfac); }
// addition: only works for other parameters
      Parameters & operator +=(Parameters const& other) {
	tdata_ += other.tdata_;
	return *this;
      }
      void print(std::ostream& ost=std::cout,int detail=0) const {
	ost << "Parameters params " << parameters() << std::endl;
	if(detail > 1)
	  ost << "covariance " << covariance() << std::endl;
      }
    private:
      Data tdata_; // data payload
  };
  std::ostream& operator << (std::ostream& ost, Parameters const& pdata);

}
#endif
