#ifndef KinKal_PData_hh
#define KinKal_PData_hh
//
//  Data object describing fit parameters
//  used as part of the kinematic kalman fit
//
#include "KinKal/TData.hh"
#include "KinKal/WData.hh"
#include <ostream>
namespace KinKal {
  template <size_t DDIM> class WData;
  template <size_t DDIM> class PData {
    public:
      constexpr static size_t PDim() { return DDIM; }
    // forward the typedefs
      typedef TData<DDIM> TDATA;
      typedef WData<DDIM> WDATA;
      typedef typename TDATA::DVEC DVEC;
      typedef typename TDATA::DMAT DMAT;
      // construct from vector and matrix
      PData(DVEC const& pars, DMAT const& pcov) : tdata_(pars,pcov) {}
      PData(DVEC const& pars) : tdata_(pars) {}
      PData(WDATA const& wdata) : tdata_(wdata.tData(),true) {}
      PData() {}
      // accessors; just re-interpret the base class accessors
      DVEC const& parameters() const { return tdata_.vec(); }
      DMAT const& covariance() const { return tdata_.mat(); }
      DVEC& parameters() { return tdata_.vec(); }
      DMAT& covariance() { return tdata_.mat(); }
      TDATA const& tData() const { return tdata_; }
      TDATA& tData() { return tdata_; }
      // scale the matrix
      void scale(double sfac) { tdata_.scale(sfac); }
      // diagnostic access to diagonal vector
      DVEC diagonal() const { 
	DVEC retval;
	for(size_t idim=0;idim < DDIM; idim++){
	  retval(idim) = sqrt(tdata_.mat()(idim,idim));
	}
	return retval;
      }
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
      TDATA tdata_; // data payload
  };

  template<size_t DDIM> std::ostream& operator << (std::ostream& ost, PData<DDIM> const& pdata) {
    pdata.print(ost,0);
    return ost;
  }
}
#endif
