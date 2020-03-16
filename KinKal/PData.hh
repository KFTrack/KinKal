#ifndef KinKal_PData_hh
#define KinKal_PData_hh
//
//  Data object describing fit parameters
//  used as part of the kinematic kalman fit
//
#include "KinKal/TData.hh"
#include <iostream>
namespace KinKal {
  template <size_t DDIM> class PData : public TData<DDIM> {
    public:
    // forward the typedefs
      typedef typename TData<DDIM>::DVEC DVEC;
      typedef typename TData<DDIM>::DMAT DMAT;
      // construct from vector and matrix
      PData(DVEC const& pars, DMAT const& pcov) : TData<DDIM>(pars,pcov) {}
      PData(DVEC const& pars) : TData<DDIM>(pars) {}
      PData(TData<DDIM> const& tdata,bool inv=false) : TData<DDIM>(tdata,inv) {}
      PData() : TData<DDIM>() {}
      // accessors; just re-interpret the base class accessors
      DVEC const& parameters() const { return TData<DDIM>::vec(); }
      DMAT const& covariance() const { return TData<DDIM>::mat(); }
      DVEC& parameters() { return TData<DDIM>::vec(); }
      DMAT& covariance() { return TData<DDIM>::mat(); }
      // addition: only works for other parameters
      PData & operator +=(PData const& other) {
	TData<DDIM>::operator +=(other);
	return *this;
      }
      PData & operator -=(PData const& other) {
	TData<DDIM>::operator -=(other);
	return *this;
      }
  };
}
#endif
