#ifndef KinKal_LocalBasis_hh
#define KinKal_LocalBasis_hh
//
// logical defintion for a local basis defined WRT the local momentum direction (momdir).  phidir is perpendicular to
// the momentum and Z, perpdir is perpendicular to momhat and perpdir.  This basis is used to define directions and
// derivatives in the kinematic kalman fit
//
#include <string>
namespace KinKal {
  struct LocalBasis {
    enum LocDir {momdir=0,perpdir,phidir,ndir};
    static std::string directionName(LocDir tdir) {
      switch (tdir) {
	case momdir:
	  return std::string("momdir");
	case perpdir:
	  return std::string("perpdir");
	case phidir:
	  return std::string("phidir");
	default:
	  return std::string("unknown");
      }
    }
  };
}
#endif
