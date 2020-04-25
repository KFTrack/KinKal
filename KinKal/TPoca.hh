#ifndef KinKal_TPoca_hh
#define KinKal_TPoca_hh
///
//  This functor class finds the (spacetime) points of closest approach between a particle and sensor trajectory
//  Both trajectories must satisfy the 'TTraj' interface
//  Concrete instances are specializations and must be implemented explicity for each trajectory pair
//  Used as part of the kinematic Kalman fit
//
#include "KinKal/TPocaBase.hh"
#include <ostream>

namespace KinKal {
  // Class to calculate DOCA and TOCA using time parameterized trajectories.
  // Templated on the types of trajectories. The actual implementations must be specializations for particular trajectory classes.
  template<class KTRAJ, class STRAJ> class TPoca : public TPocaBase {
    public:
    // forward the base interface.
      typedef typename KTRAJ::DVEC DVEC; // forward derivative type from the particle trajectory
      // derviatives of TOCA and DOCA WRT particle trajectory parameters
      DVEC const& dDdP() const { return dDdP_; }
      DVEC const& dTdP() const { return dTdP_; }
      // construct from the particle and sensor trajectories; POCA is computed on construction, using possible hints
      // default precision = 1 Ps (~300 um) along the trajectories
      TPoca(KTRAJ const& ktraj, STRAJ const& straj, TPocaHint const& hint=TPocaHint(), float precision=0.001);
      // accessors
      KTRAJ const& particleTraj() const { return *ktraj_; }
      STRAJ const& sensorTraj() const { return *straj_; }
      bool inRange() const { return particleTraj().inRange(particleToca()) && sensorTraj().inRange(sensorToca()); }
      void print(std::ostream& ost=std::cout,int detail=0) const;
    private:
      const KTRAJ* ktraj_; // kinematic particle trajectory
      const STRAJ* straj_; // sensor trajectory
      DVEC dDdP_; // derivative of DOCA WRT Parameters
      DVEC dTdP_; // derivative of Dt WRT Parameters
  };

  template<class KTRAJ, class STRAJ> void TPoca<KTRAJ,STRAJ>::print(std::ostream& ost,int detail) const {
    ost << "TPoca " << TPocaBase::statusName(status()) << " Doca " << doca() << " +- " << sqrt(docaVar())
      << " dToca " << deltaT() << " +- " << sqrt(tocaVar()) << " cos(theta) " << dirDot() << " Precision " << precision() << std::endl;
    if(detail > 0)
      ost << "Particle Poca " << particlePoca() << " Sensor Poca " << sensorPoca() << std::endl;
    if(detail > 1)
      ost << "dDdP " << dDdP() << " dTdP " << dTdP() << std::endl;
    if(detail > 2){
      ost << "Particle ";
      particleTraj().print(ost,0);
      ost << "Sensor ";
      sensorTraj().print(ost,0);
    }
  }

}
#endif
