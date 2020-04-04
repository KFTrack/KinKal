#ifndef KinKal_WData_hh
#define KinKal_WData_hh
//
//  Data object describing weight-space information
//  used as part of the kinematic kalman fit
//
#include "KinKal/TData.hh"
#include <ostream>
namespace KinKal {
  template <size_t DDIM> class WData : public TData<DDIM> {
    public:
    // forward the typedefs
      typedef typename TData<DDIM>::DVEC DVEC;
      typedef typename TData<DDIM>::DMAT DMAT;
      // construct from vector and matrix
      WData(DVEC const& wvec, DMAT const& wmat) : TData<DDIM>(wvec,wmat) {}
      WData(DVEC const& wvec) : TData<DDIM>(wvec) {}
      WData(TData<DDIM> const& tdata,bool inv=false) : TData<DDIM>(tdata,inv) {}
      WData() : TData<DDIM>() {}
      // accessors; just re-interpret the base class accessors
      DVEC const& weightVec() const { return TData<DDIM>::vec(); }
      DMAT const& weightMat() const { return TData<DDIM>::mat(); }
      DVEC& weightVec() { return TData<DDIM>::vec(); }
      DMAT& weightMat() { return TData<DDIM>::mat(); }
      // addition: only works for other weights
      WData & operator +=(WData const& other) {
	TData<DDIM>::operator +=(other);
	return *this;
      }
      void print(std::ostream& ost=std::cout,int detail=0) const {
	ost << "WData wVec " << weightVec() << std::endl;
	if(detail > 0)
	  ost << "weight " << weightMat() << std::endl;
      }
  };
  template<size_t DDIM> std::ostream& operator << (std::ostream& ost, WData<DDIM> const& wdata) {
    wdata.print(ost,0);
    return ost;
  }
}
#endif
