#ifndef BTrk_KinKal_TrkDir_hh 
#define BTrk_KinKal_TrkDir_hh 
// -----------------------------------------------------------------------
// enum defining propagation directions: positive is increasing time, negative is decreasing
//------------------------------------------------------------------------

namespace KinKal {
  enum TrkDir { ttneg=0, ttpos=1  };

  std::ostream& operator<<(std::ostream& os, const TrkDir& tdir) {
    switch(tdir) {
      case ttneg:
	os << " Decreasing time";
	break;
      case ttpos:
  	os << " Increasing time";
	break;
      default:
	os << " Error";
	break;
    }
    return os;
  }

}
#endif
