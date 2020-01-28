#ifndef KinKal_TTraj_hh
#define KinKal_TTraj_hh
//
//  Base class for a trajectory (1-dimensional path in space) with time as parametric variable. 
//  Used as part of the kinematic Kalman fit
//  Templated on the dimension of the parameter vector
//  The geometric interpretation of the parameters is defined in the subclasses
//
#include "Math/SVector.h"
#include "Math/SMatrix.h"
#include "BTrk/KinKal/Types.hh"
namespace KinKal {
  template <size_t PDIM> class TTraj {
    public:
      // define the parameter types
      typedef ROOT::Math::SVector<double,PDIM> PVec; // vector of parameters
      typedef ROOT::Math::SMatrix<double,PDIM,PDIM,ROOT::Math::MatRepSym<double,PDIM> > PMat;  // covariance matrix of parameters
      typedef ROOT::Math::SMatrix<double,PDIM,1> PDer; // derivative of parameters

      // construct from parameters and covariance
      TTraj(PVec const& pars, PMat const& pcov) : pars_(pars), pcov_(pcov) {}
      // direct accessors
      PVec const& params() const { return pars_;}
      PMat const& covar() const { return pcov_; }

      // geometric accessors
      virtual void position(Vec4& pos) const =0; // position as a function of time
      virtual void position(double time, Vec3& pos) const=0;
      virtual void velocity(double time, Vec3& vel) const =0; // velocity vector in mm/ns
      virtual void direction(double time, Vec3& dir) const =0; // unit vector in the direction of positive time

    protected:
      TTraj(){}
      PVec pars_; // parameters for this trajectory
      PMat pcov_; // covariance for this trajectory
  };
}
#endif

