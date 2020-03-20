#ifndef KinKal_MIsect_hh
#define KinKal_MIsect_hh
//
//  Struct to describe the intersection of a kinematic trajectory with a piece of material
//
#include "MatEnv/DetMaterial.hh"
namespace KinKal {
  struct MIsect {
    MatEnv::DetMaterial const& dmat_; // material
    double plen_; // path length through this material
    MIsect(MatEnv::DetMaterial const& dmat,double plen) : dmat_(dmat), plen_(plen) {}
  };
}
#endif

