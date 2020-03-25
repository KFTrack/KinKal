#ifndef KinKal_DMat_hh
#define KinKal_DMat_hh
//
//  Interface describing a piece of detector material
//
#include "KinKal/MatXing.hh"
#include "KinKal/TPocaBase.hh"
namespace KinKal {
  class DMat {
    public:
      virtual void intersect(TPocaBase const& poca,std::vector<MatXing>& mxings) const =0;
  };
}
#endif
