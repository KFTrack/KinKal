#ifndef KinKal_MatXing_hh
#define KinKal_MatXing_hh
//
//  Struct to describe a path crossing a piece of material
//
#include "MatEnv/DetMaterial.hh"
#include <vector>
namespace KinKal {
  struct MatXing {
    MatEnv::DetMaterial const& dmat_; // material
    float plen_; // path length through this material
    MatXing(MatEnv::DetMaterial const& dmat,float plen) : dmat_(dmat), plen_(plen) {}
  };
  typedef std::vector<MatXing> MatXingCol;
}
#endif

