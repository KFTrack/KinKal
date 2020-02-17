#ifndef KinKal_PData_hh
#define KinKal_PData_hh
//
//  Data object describing fit parameters
//  used as part of the kinematic kalman fit
//
#include "BTrk/KinKal/TData.hh"
//#include "BTrk/KinKal/WData.hh"
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
      PData() : TData<DDIM>() {}
      // construct from a WData object: this requires inversion and may result in an unusable object
//      PData(WData const& wdata) : PData(wdata,true) {}
      PData(TData<DDIM> const& tdata,bool invert=false) : TData<DDIM>(tdata,invert) {}
      // accessors; just re-interpret the base class accessors
      DVec const& parameters() const { return TData<DDIM>::vec(); }
      DMat const& covariance() const { return TData<DDIM>::mat(); }
      DVec& parameters() { return TData<DDIM>::vec(); }
      DMat& covariance() { return TData<DDIM>::mat(); }
      // addition: only works for other parameters
      PData & operator +=(PData const& other) {
	parameters() += other.parameters();
	covariance() += other.covariance();
	return *this;
      }
  };
}
#endif
