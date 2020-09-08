#ifndef KinKal_PData_hh
#define KinKal_PData_hh
//
//  Data object describing fit parameters
//  used as part of the kinematic kalman fit
//
#include "KinKal/TData.hh"
#include <ostream>
namespace KinKal {
  class WData;
  class PData {
    public:
      // construct from vector and matrix
      PData(DVEC const& pars, DMAT const& pcov) : tdata_(pars,pcov) {}
      PData(DVEC const& pars) : tdata_(pars) {}
      PData(WData const& wdata);
      PData() {}
      // accessors; just re-interpret the base class accessors
      DVEC const& parameters() const { return tdata_.vec(); }
      DMAT const& covariance() const { return tdata_.mat(); }
      DVEC& parameters() { return tdata_.vec(); }
      DMAT& covariance() { return tdata_.mat(); }
      TData const& tData() const { return tdata_; }
      TData& tData() { return tdata_; }
      // scale the matrix
      void scale(double sfac) { tdata_.scale(sfac); }
// addition: only works for other parameters
      PData & operator +=(PData const& other) {
	tdata_ += other.tdata_;
	return *this;
      }
      void print(std::ostream& ost=std::cout,int detail=0) const {
	ost << "PData params " << parameters() << std::endl;
	if(detail > 1)
	  ost << "covariance " << covariance() << std::endl;
      }
    private:
      TData tdata_; // data payload
  };
  std::ostream& operator << (std::ostream& ost, PData const& pdata);

}
#endif
