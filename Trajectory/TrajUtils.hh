#ifndef KinKal_TrajUtils_hh
#define KinKal_TrajUtils_hh
//
// Utility functions related to trajectories
//
#include "KinKal/General/Vectors.hh"
#include "KinKal/General/MomBasis.hh"
namespace KinKal {
  // Create the transform matrix from the (trajectory-specific) momentum basis to global cartesian
  template <class KTRAJ> SSMAT momToGlobal(KTRAJ const& ktraj,double time) {
    // loop over the momentum change basis directions and create the transform matrix between that and global cartesian basis
    SSMAT dmdxyz; // momentum -> cartesian conversion matrix
    for(int idir=0;idir<MomBasis::ndir; idir++) {
      auto mdir = ktraj.direction(time,static_cast<MomBasis::Direction>(idir));
      SVEC3 vmdir(mdir.X(), mdir.Y(), mdir.Z());
      dmdxyz.Place_in_col(vmdir,0,idir);
    }
    return dmdxyz;
  }

}
#endif

