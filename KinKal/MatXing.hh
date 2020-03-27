#ifndef KinKal_MatXing_hh
#define KinKal_MatXing_hh
//
//  Struct to describe a path crossing a piece of material
//
#include "MatEnv/DetMaterial.hh"
namespace KinKal {
  struct MatXing {
    MatEnv::DetMaterial const& dmat_; // material
    double plen_; // path length through this material
    MatXing(MatEnv::DetMaterial const& dmat,double plen) : dmat_(dmat), plen_(plen) {}
  };
}
#endif

