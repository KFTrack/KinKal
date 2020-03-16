#ifndef KinKal_DMat_hh
#define KinKal_DMat_hh
//
//  Interface describing a piece of detector material
//
#include "KinKal/MIsect.hh"
#include "KinKal/TPocaBase.hh"
namespace KinKal {
  class DMat {
    public:
      virtual void intersect(TPocaBase const& poca,std::vector<MIsect>& misects) const =0;
  };
}
#endif
