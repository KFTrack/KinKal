#ifndef KinKal_Parameters_hh
#define KinKal_Parameters_hh
//
//  Data object describing fit parameters
//  used as part of the kinematic kalman fit
//
#include "KinKal/FitData.hh"
#include <ostream>
namespace KinKal {
  class Weights;
  class Parameters {
    public:
      // construct from vector and matrix
      Parameters(DVEC const& pars, DMAT const& pcov) : fitdata_(pars,pcov) {}
      Parameters(DVEC const& pars) : fitdata_(pars) {}
      Parameters(Weights const& wdata);
      Parameters() {}
      // accessors; just re-interpret the base class accessors
      DVEC const& parameters() const { return fitdata_.vec(); }
      DMAT const& covariance() const { return fitdata_.mat(); }
      DVEC& parameters() { return fitdata_.vec(); }
      DMAT& covariance() { return fitdata_.mat(); }
      FitData const& fitData() const { return fitdata_; }
      FitData& fitData() { return fitdata_; }
      // scale the matrix
      void scale(double sfac) { fitdata_.scale(sfac); }
// addition: only works for other parameters
      Parameters & operator +=(Parameters const& other) {
	fitdata_ += other.fitdata_;
	return *this;
      }
      void print(std::ostream& ost=std::cout,int detail=0) const {
	ost << "Parameters params " << parameters() << std::endl;
	if(detail > 1)
	  ost << "covariance " << covariance() << std::endl;
      }
    private:
      FitData fitdata_; // data payload
  };
  std::ostream& operator << (std::ostream& ost, Parameters const& pdata);

}
#endif
