#ifndef KinKal_WData_hh
#define KinKal_WData_hh
//
//  Data object describing weight-space information
//  used as part of the kinematic kalman fit
//
#include "KinKal/TData.hh"
#include <iostream>
namespace KinKal {
  template <size_t DDIM> class WData : public TData<DDIM> {
    public:
    // forward the typedefs
      typedef typename TData<DDIM>::DVec DVec;
      typedef typename TData<DDIM>::DMat DMat;
      // construct from vector and matrix
      WData(DVec const& wvec, DMat const& wmat) : TData<DDIM>(wvec,wmat) {}
      WData(DVec const& wvec) : TData<DDIM>(wvec) {}
      WData(TData<DDIM> const& tdata) : TData<DDIM>(tdata) {}
      WData() : TData<DDIM>() {}
      // construct from a PData object: this requires inversion and may result in an unusable object
//      WData(PData const& pdata) : WData(pdata,true) {}
      WData(TData<DDIM> const& tdata,bool invert=false) : TData<DDIM>(tdata,invert) {}
      // accessors; just re-interpret the base class accessors
      DVec const& weightVec() const { return TData<DDIM>::vec(); }
      DMat const& weightMat() const { return TData<DDIM>::mat(); }
      DVec& weightVec() { return TData<DDIM>::vec(); }
      DMat& weightMat() { return TData<DDIM>::mat(); }
      // addition: only works for other weights
      WData & operator +=(WData const& other) {
	weightVec() += other.weightVec();
	weightMat() += other.weightMat();
	return *this;
      }
      void invert(TData<DDIM> const& other) { TData<DDIM>::invert(other); }
  };
}
#endif
