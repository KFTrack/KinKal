#ifndef BTrk_KinKal_TrkDirection_hh 
#define BTrk_KinKal_TrkDirection_hh 
// -----------------------------------------------------------------------
// enum defining directions: positive is increasing time, negative is decreasing
//------------------------------------------------------------------------

namespace KinKal {
  enum trkDirection { ttneg=-1, pos=1  };

  std::ostream& operator<<(std::ostream& os, const trkDirection& tdir) {
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
