#ifndef KinKal_MaterialXing_hh
#define KinKal_MaterialXing_hh
//
//  Struct to describe a path crossing a piece of material
//
#include "KinKal/MatEnv/DetMaterial.hh"
#include <vector>
namespace KinKal {
  struct MaterialXing {
    MatEnv::DetMaterial const& dmat_; // material
    double plen_; // path length through this material
    MaterialXing(MatEnv::DetMaterial const& dmat,double plen) : dmat_(dmat), plen_(plen) {}
  };
  typedef std::vector<MaterialXing> MaterialXingCol;
}
#endif

