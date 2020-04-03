
#include "KinKal/LRAmbig.hh"
namespace KinKal {
  std::ostream& operator <<(std::ostream& ost, LRAmbig const& lrambig) {
    if(lrambig == LRAmbig::left)
      ost << "left";
    else if (lrambig == LRAmbig::right)
      ost << "right";
    else if (lrambig == LRAmbig::null)
      ost << "null";
    else
      ost <<"Unknown";
    return ost;
  }
}
