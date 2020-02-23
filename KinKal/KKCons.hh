#ifndef KinKal_KKCons_hh
#define KinKal_KKCons_hh
//
// Class to add a constraint to the fit
// This effect provides information content and is processed in weight space 
//
#include "KinKal/KKWeight.hh"
#include <array>

namespace KinKal {
  template<class KTRAJ> class KKCons : public KKWeight<KTRAJ> {
    public:
      typedef std::array<bool,KTRAJ::PDATA::PDim()> PMASK; // parameter mask
      typedef typename KKEff::WDATA WDATA; // forward the typedef
      // process this site given the adjacent site
      virtual bool process(KKEff const& other,TDir tdir) override;

      // construct from (masked) parameters
      KKCons(double time, PDATA const& params, PMASK CONST& pmask);
      virtual ~KKCons(){}
  };

}
#endif
