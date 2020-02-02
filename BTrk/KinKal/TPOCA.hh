#ifndef KinKal_TPOCA_hh
#define KinKal_TPOCA_hh
///
//  This functor class finds the (spacetime) points of closest approach between two TTrajs.
//  Concrete instances are specializations and must be implemented explicity for
//  each trajectory type pair.
//  Used as part of the kinematic Kalman fit
//
#include "BTrk/KinKal/Types.hh"
#include "BTrk/KinKal/TTraj.hh"
#include "Math/SMatrix.h"
#include <typeinfo>
#include <iostream>
#include <array>

namespace KinKal {
  // base class for POCA calculations
  class TPOCABase {
    public:
      enum pstat{converged=0,unconverged,outsiderange,pocafailed,derivfailed,unknown};
      //accessors
      Vec4 const& poca(size_t itraj) const { return poca_[itraj]; }
      Vec4 const& poca0() const { return poca_[0]; }
      Vec4 const& poca1() const { return poca_[1]; }
      pstat status() const { return status_; }
      double doca() const { return doca_; }
      double precision() const { return precision_; }
      // utility functions
      void delta(Vec4& ds) const { ds = poca_[1]-poca_[0]; }
      double dt() const { return poca_[1].T() - poca_[0].T(); }
      void doca(Vec3& ds) const { Vec3 p0(poca_[0]); Vec3 p1(poca_[1]); ds = p1-p0; }
      bool usable() const { return status_ != pocafailed && status_ != unknown; }
      virtual TTraj const& ttraj0() const { return *ttraj_[0]; }
      virtual TTraj const& ttraj1() const { return *ttraj_[1]; }
      // actions
      virtual void findPOCA() = 0;
      void testRange() {
	if(usable() && ((!ttraj_[0]->inRange(poca_[0].T())) || (!ttraj_[1]->inRange(poca_[1].T()))) )
	  status_ = outsiderange;
      }
    protected:
      pstat status_; // status of computation
      double doca_; // geometric distance of closest approach (cached)
      double precision_; // precision used to define convergence
      std::array<Vec4,2> poca_; // spacetime points at GEOMETRIC closest approach
      std::array<const TTraj*,2> ttraj_; // base class pointers to trajs
      // override default constructor:
      // default precision = 10 um on DOCA
      TPOCABase(TTraj const& traj0,TTraj const& traj1, double precision=0.01) : status_(unknown), doca_(-1.0), precision_(precision),
      ttraj_{&traj0,&traj1} {}
      void reset() {status_ = unknown;}
  };

  // Class to calculate POCA using time parameterized trajectories.
  // Templated on the types of trajectories. The actual implementations must be specializations for particular trajectory classes.
  template<class T0, class T1> class TPOCA : public TPOCABase {
    public:
      // construct from a pair of trajs
      TPOCA(T0 const& t0, T1 const& t1, double precision=0.01) : TPOCABase(t0,t1,precision)  { }
      // accessors
      T0 const& ttraj0() const override { return static_cast<T0 const&>(TPOCABase::ttraj0()); }
      T1 const& ttraj1() const override { return static_cast<T1 const&>(TPOCABase::ttraj1()); }
      // actions
      // Find POCA (or fail trying);
      virtual void findPOCA() override {
	// no default implementation so throw: this must be specialized to be useful
	std::cout << "No findPOCA implementation for types " << typeid(T0).name()
	  << " and " << typeid(T1).name() << std::endl; // should throw FIXME!
	status_ = TPOCABase::pocafailed;
      }
  };

  // Add the derivatives of POCA WRT the parameters of the 1st traj to the above.
  template<class T0, class T1> class TDPOCA : public TPOCA<T0,T1> {
    public:
      typedef ROOT::Math::SMatrix<double,T0::NParams(),1> DMat; // derivative type, dimensioned on the 0th traj parameters
      TDPOCA(T0 const& t0, T1 const& t1, double precision=0.01) : TPOCA<T0,T1>(t0,t1,precision){}
      TDPOCA(TPOCA<T0,T1> const& tpoca) : TPOCA<T0,T1>(tpoca) {} // 'upgrade' constructor
      DMat const& derivs() const { return dDdP_; }
      // find POCA and it's derivatives WRT the 1st trajectory parameters.
      virtual void findPOCA() override {
	// Start by calling the base class version
	TPOCABase::pstat status = TPOCA<T0,T1>::FindPOCA();
	if(TPOCABase::usable()){
	  fillDerivatives();
	}
      }
      void fillDerivatives() { 
	std::cout << "No fillDerivatives implementation for types " << typeid(T0).name()
	  << " and " << typeid(T1).name() << std::endl; // should throw FIXME!
      }
    private:
      DMat dDdP_; // derivative of DOCA WRT Parameters
  };

}
#endif
