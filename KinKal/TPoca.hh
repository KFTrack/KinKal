#ifndef KinKal_TPoca_hh
#define KinKal_TPoca_hh
///
//  This functor class finds the (spacetime) points of closest approach between two TTrajs.
//  Concrete instances are specializations and must be implemented explicity for
//  each trajectory type pair.
//  Used as part of the kinematic Kalman fit
//
#include "KinKal/TTraj.hh"
#include "KinKal/TPocaBase.hh"
#include "Math/SMatrix.h"
#include <typeinfo>
#include <stdexcept>
#include <array>

namespace KinKal {

  // Class to calculate POCA using time parameterized trajectories.
  // Templated on the types of trajectories. The actual implementations must be specializations for particular trajectory classes.
  template<class T0, class T1> class TPoca : public TPocaBase {
    public:
      // construct from a pair of trajs; POCA is computed on construction
      TPoca(T0 const& t0, T1 const& t1, double precision=0.01);
      TPoca() {} 
      // accessors
      T0 const& ttraj0() const override { return static_cast<T0 const&>(TPocaBase::ttraj0()); }
      T1 const& ttraj1() const override { return static_cast<T1 const&>(TPocaBase::ttraj1()); }
  };

  // Compute POCA and the derivatives of DOCA WRT the parameters of T0
  template<class T0, class T1> class TDPoca : public TPoca<T0,T1> {
    public:
      typedef typename T0::PDER PDER; // forward derivative type from the 0th traj parameters
      TDPoca(T0 const& t0, T1 const& t1, double precision=0.01);
      TDPoca(TPoca<T0,T1> const& tpoca); // 'upgrade' a regular POCA 
      TDPoca() {} 
      PDER const& dDdP() const { return dDdP_; }
      PDER const& dTdP() const { return dTdP_; } 
    private:
      PDER dDdP_; // derivative of DOCA WRT Parameters
      PDER dTdP_; // derivative of Dt WRT Parameters
  };

}
#endif
