#ifndef kinkal_tvec_hh
#define kinkal_tvec_hh
//
//  Parameter types to describe trajectories 
//  used as part of the kinematic kalman fit
//  templated on the dimension of the parameter vector
//
#include "Math/SVector.h"
#include "Math/SMatrix.h"
namespace KinKal {
  template <size_t PDIM> struct TPars {
    // define the parameter types
    typedef ROOT::Math::SVector<double,PDIM> PVec; // vector
    typedef ROOT::Math::SMatrix<double,PDIM,PDIM,ROOT::Math::MatRepSym<double,PDIM> > PMat;  // matrix

    // construct from vector and matrix
    TPars(PVec const& pars, PMat const& pcov) : vec_(pars), mat_(pcov) {}
    TPars() {} // not sure why the compiler can't write this
    PVec vec_; // parameter or weight vector
    PMat mat_; // covariance or weight matrix
  };
}
#endif

