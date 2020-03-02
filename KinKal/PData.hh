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
      typedef typename TData<DDIM>::DVec DVec;
      typedef typename TData<DDIM>::DMat DMat;
      // construct from vector and matrix
      PData(DVec const& pars, DMat const& pcov) : TData<DDIM>(pars,pcov) {}
      PData(DVec const& pars) : TData<DDIM>(pars) {}
      PData(TData<DDIM> const& tdata,bool inv=false) : TData<DDIM>(tdata,inv) {}
      PData() : TData<DDIM>() {}
      // accessors; just re-interpret the base class accessors
      DVec const& parameters() const { return TData<DDIM>::vec(); }
      DMat const& covariance() const { return TData<DDIM>::mat(); }
      DVec& parameters() { return TData<DDIM>::vec(); }
      DMat& covariance() { return TData<DDIM>::mat(); }
      // addition: only works for other parameters
      PData & operator +=(PData const& other) {
	TData<DDIM>::operator +=(other);
	return *this;
      }
  };
}
#endif
