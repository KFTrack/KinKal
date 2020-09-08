#ifndef KinKal_KKCons_hh
#define KinKal_KKCons_hh
//
// Class to add a constraint to the fit
// This effect provides information content and is processed in weight space 
//
#include "KinKal/KKWEff.hh"
#include <array>

namespace KinKal {
  template<class KTRAJ> class KKCons : public KKWEff<KTRAJ> {
    public:
      using PMASK = std::array<bool,KTRAJ::PData::PDim()>; // parameter mask
      // process this effect given the adjacent effect
      virtual bool process(KKEff const& other,TDir tdir) override;
      // construct from (masked) parameters
      KKCons(double time, PData const& params, PMASK CONST& pmask);
      virtual ~KKCons(){}
  };

}
#endif
