#include "KinKal/TDir.hh"
namespace KinKal {
  TDir& operator ++ (TDir& tdir) {
    if (tdir == TDir::end) throw std::out_of_range("TDir& operator ++");
    tdir = TDir(static_cast<std::underlying_type<TDir>::type>(tdir) + 1);
    return tdir;
  }
  std::ostream& operator <<(std::ostream& ost, TDir const& tdir) {
    if(tdir == TDir::forwards)
      ost << "forwards";
    else if (tdir == TDir::backwards)
      ost << "backwards";
    else
      ost <<"Unknown";
    return ost;
  }
}
