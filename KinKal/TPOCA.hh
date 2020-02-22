#ifndef KinKal_TPOCA_hh
#define KinKal_TPOCA_hh
///
//  This functor class finds the (spacetime) points of closest approach between two TTrajs.
//  Concrete instances are specializations and must be implemented explicity for
//  each trajectory type pair.
//  Used as part of the kinematic Kalman fit
//
#include "KinKal/TTraj.hh"
#include "KinKal/TPOCABase.hh"
#include "Math/SMatrix.h"
#include <typeinfo>
#include <stdexcept>
#include <array>

namespace KinKal {

  // Class to calculate POCA using time parameterized trajectories.
  // Templated on the types of trajectories. The actual implementations must be specializations for particular trajectory classes.
  template<class T0, class T1> class TPOCA : public TPOCABase {
    public:
      // construct from a pair of trajs; POCA is computed on construction
      TPOCA(T0 const& t0, T1 const& t1, double precision=0.01);
      // accessors
      T0 const& ttraj0() const override { return static_cast<T0 const&>(TPOCABase::ttraj0()); }
      T1 const& ttraj1() const override { return static_cast<T1 const&>(TPOCABase::ttraj1()); }
  };

  // Compute POCA and the derivatives of DOCA WRT the parameters of T0
  template<class T0, class T1> class TDPOCA : public TPOCA<T0,T1> {
    public:
      typedef typename T0::PDer PDer; // forward derivative type from the 0th traj parameters
      TDPOCA(T0 const& t0, T1 const& t1, double precision=0.01);
      TDPOCA(TPOCA<T0,T1> const& tpoca); // 'upgrade' a regular POCA 
      PDer const& dDdP() const { return dDdP_; }
      PDer const& dTdP() const { return dTdP_; } 
    private:
      PDer dDdP_; // derivative of DOCA WRT Parameters
      PDer dTdP_; // derivative of Dt WRT Parameters
  };

}
#endif
