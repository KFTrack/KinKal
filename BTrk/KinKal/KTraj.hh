#ifndef KinKal_KTraj_hh
#define KinKal_KTraj_hh
//
//  Base class for a trajectory used in the kinematic Kalman fit
//  Templated on the dimension of the parameter vector
//  The geometric and kinematic interpretation of the parameters is defined in the subclasses
//
#include "Math/SVector.h"
#include "Math/SMatrix.h"
#include "BTrk/KinKal/Types.hh"
#include "BTrk/KinKal/Constants.hh"
namespace KinKal {
  template <size_t PDIM> class KTraj {
    public:
      // define local basis vector indices; along and perp to the local momentum.  theta2 is also perpendicular to z
      enum trajdir {momdir=0,theta1,theta2};
      // unit vectors in the different directions
      virtual void dirVector(trajdir dir,double time,Pos3& unit) const = 0;
      
      // define the parameter types
      typedef ROOT::Math::SVector<double,PDIM> PVec; // vector of parameters
      typedef ROOT::Math::SMatrix<double,PDIM,PDIM,ROOT::Math::MatRepSym<double,PDIM> > PMat;  // covariance matrix of parameters
      typedef ROOT::Math::SMatrix<double,PDIM,1> PDer; // derivative of parameters

    // direct accessors
      PVec const& params() const { return pars_;}
      PMat const& covar() const { return pcov_; }
      double mass() const { return mass_;}
      int charge() const { return charge_;}

    // geometric accessors
      virtual void position(Pos4& pos) const =0; // position as a function of time
      virtual void momentum(double t,Mom4& mom) const =0; // momentum as a function of time
      void momentum(Pos4 const& pos, Mom4& mom) const { return momentum(pos.T(),mom); }

    // Parameter derivatives for a change in momentum, along the different directions
      virtual void paramDeriv(trajdir dir, double time, PDer& der) const = 0;

    protected:
      KTraj(double mass, int charge) : mass_(mass), charge_(charge) {}
      // construct from particle properties and parameters 
      KTraj(PVec const& pars, PMat const& pcov,double mass,int charge) : pars_(pars), pcov_(pcov), mass_(mass), charge_(charge) {}
      PVec pars_; // parameters for this trajectory
      PMat pcov_; // covariance for this trajectory
      // kinematic parameters
      double mass_;  // in units of MeV
      int charge_; // charge in units of proton charge
    private:
      KTraj() : mass_(-1.0), charge_(0) {} //  invalid
  };
}
#endif
