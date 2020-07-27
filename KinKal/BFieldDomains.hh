#ifndef KinKal_BFieldDomains_hh
#define KinKal_BFieldDomains_hh
//
// Define domains of the BField for a specific trajectory over which the impact of the
// changing BField on the particle compared to the nominal average in the domain is
// below the specified threshold.
// Used to model BField inhomogenity in the KinKal track fit
//
#include "KinKal/TRange.hh"
#include "KinKal/Vectors.hh"
#include <vector>
namespace KinKal {
  struct BFDomain {
    TRange trange_; // time range associated with this domain
      Vec3 bavg_; // average BField in this domain
      Vec3 dmom_; // integrated momentum difference over this domain of the particle trajectory for the true BField Map compared to the nominal over this range
  };
  template <class KTRAJ> class BFieldDomains {
    public:
      typedef BDoms std::vector<BFDomain>;
      BFieldDomains(KTRAJ const& ktraj);
    private:
      KTRAJ const& ktraj_; // kinematic trajectory used to calculate these domains
      BDoms bdomains_; // domains for this trajectory
  };

}
#endif
