#ifndef KinKal_DMat_hh
#define KinKal_DMat_hh
//
//  Interface describing the material of a detector piece
//
#include "KinKal/MatXing.hh"
#include "KinKal/TPocaBase.hh"
namespace KinKal {
  class DMat {
    public:
    // find the material crossings give TPOCA from a hit
      virtual void findXings(TPocaBase const& poca,std::vector<MatXing>& mxings) const =0;
  };
}
#endif
