#ifndef KinKal_KKConstraint_hh
#define KinKal_KKConstraint_hh
//
// Class to add a constraint to the fit
// This effect provides information content and is processed in weight space 
//
#include "BTrk/KinKal/KKWeight.hh"
#include <array>

namespace KinKal {
  template<class KTRAJ> class KKConstraint : public KKWeight<KTRAJ> {
    public:
      typedef std::array<bool,KTRAJ::PDATA::PDim()> PMASK; // parameter mask
      typedef typename KKEffect::WDATA WDATA; // forward the typedef
      // process this site given the adjacent site
      virtual bool process(KKEffect const& other,TDir tdir) override;
    
  // construct from (masked) parameters
   KKConstraint(double time, PDATA const& params, PMASK CONST& pmask);
  };

}
#endif
