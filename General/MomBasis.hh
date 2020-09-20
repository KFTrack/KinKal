#ifndef KinKal_MomBasis_hh
#define KinKal_MomBasis_hh
//
// logical defintion for a local basis defined WRT the local momentum direction (momdir_).  phidir_ is perpendicular to
// the momentum and Z, perpdir_ is perpendicular to momdir_ and perpdir_.  This basis is used to define directions and
// derivatives in the kinematic kalman fit
//
#include <string>
namespace KinKal {
  struct MomBasis {
    enum Direction {momdir_=0,perpdir_,phidir_,ndir};
    static std::string directionName(Direction tdir) {
      switch (tdir) {
	case momdir_:
	  return std::string("MomentumDirection");
	case perpdir_:
	  return std::string("PerpendicularDirection");
	case phidir_:
	  return std::string("PhiDirection");
	default:
	  return std::string("Unknown");
      }
    }
  };
}
#endif
