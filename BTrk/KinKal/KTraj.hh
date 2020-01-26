#ifndef KinKal_KTraj_hh
#define KinKal_KTraj_hh
//
//  Base class for a trajectory used in the kinematic Kalman fit
//  Templated on the dimension of the parameter vector
//  The geometric and kinematic interpretation of the parameters is defined in the subclasses
//
#include "Math/SVector.h"
#include "Math/SMatrix.h"
#include "KinKal/Types.hh"
#include "KinKal/Constants.hh"
namespace KinKal {
  class KTraj<size_t PDIM> {
    public:
      // define the parameter types
      typedef SVector<double,PDIM> PVec; // vector of parameters
      typedef SMatrix<double,PDIM,PDIM,MatRepSym<double,PDIM> > PMat;  // covariance matrix of parameters
      typedef SMatrix<double,PDIM,1> PDer; // derivative of parameters

      PVec const& params() const { return pars_;}
      PMat const& covar() const { return pcov_; }
      double mass() const { return mass_;}
      int charge() const { return charge_;}

      void position(FourV& pos) const =0; // position as a function of time
      void momentum(double t,FourV& mom) const =0; // momentum as a function of time

    protected:
      KTraj() : mass_(-1.0), charge_(0) {} //  
      // construct from particle properties and parameters 
      KTraj(PVec const& pars, PMat const& pcov,double mass,int charge) : pars_(pars), pcov_(pcov), mass_(mass), charge_(charge) {}
      PVec pars_; // parameters for this trajectory
      PMat pcov_; // covariance for this trajectory
      // kinematic parameters
      double mass_;  // in units of MeV
      int charge_; // charge in units of proton charge
  };
}
#endif
