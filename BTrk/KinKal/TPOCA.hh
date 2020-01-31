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
#include <typeinfo>
#include <typeindex>
#include <iostream>
#include <array>
#include <tuple>

namespace KinKal {

// base class for POCA calculations
  class POCABase {
    public:
      enum pstat{converged=0,unconverged,outsidelimits,failed,unknown};
      typedef std::array<double,2> trange;
      //accessors
      Vec4 const& poca(size_t itraj) const { return poca_[itraj]; }
      pstat status() const { return status_; }
      double doca() const { return doca_; }
      double precision() const { return precision_; }
      // utility functions
      void delta(Vec4& ds) const { ds = poca_[1]-poca_[0]; }
      double dt() const { return poca_[1].T() - poca_[0].T(); }
      void doca(Vec3& ds) const { Vec3 p0(poca_[0]); Vec3 p1(poca_[1]); ds = p1-p0; }
      bool usable() const { return status_ != failed && status_ != unknown; }
      // actions
      virtual void findPOCA() = 0;
    protected:
      pstat status_; // status of computation
      double doca_; // geometric distance of closest approach (cached)
      double precision_; // precision used to define convergence
      std::array<Vec4,2> poca_; // spacetime points at GEOMETRIC closest approach
      std::array<trange,2> trange_; // time range to find POCA; Negative range = all allowed
      // override default constructor:
      // default precision = 10 um on DOCA
      POCABase(trange tr1={1.0,-1.0},trange tr2={1.0,-1.0}, double precision=0.01) : status_(unknown), doca_(-1.0), precision_(precision), trange_{tr1,tr2} {}
      void reset() {status_ = unknown;}
  };

  template<class T1, class T2> class TPOCA : public POCABase {
    public:
      // construct from a pair of trajs
      TPOCA(T1 const& t1, T2 const& t2,trange tr1={1.0,-1.0},trange tr2={1.0,-1.0}, double precision=0.01) : POCABase(tr1,tr2,precision), t1_(t1), t2_(t2) { }
      // accessors
      T1 const& traj1() const { return t1_; }
      T2 const& traj2() const { return t2_; }
      // actions
      virtual void findPOCA() override {
	// no default implementation so throw: this must be specialized to be useful
	std::cout << "No findPOCA implementation for types " << typeid(T1).name()
	  << " and " << typeid(T2).name() << std::endl; // should throw FIXME!
      }
      
    // Derivatives at the specified time along the specified direction WRT the given trajectory derviatives
      // this must be templated on the type of the derivative for this 
//    template<class PDer> void findDPOCA(PDer& dervis) {
//	std::cout << "No findDPOCA implementation for types " << typeid(T1).name()
//	  << " and " << typeid(T2).name() << std::endl; // should throw FIXME!
//    }

    private:
      T1 const& t1_;
      T2 const& t2_;
  };

  
}
#endif
