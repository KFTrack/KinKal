#include "KinKal/General/TimeDir.hh"
namespace KinKal {
  TimeDir& operator ++ (TimeDir& tdir) {
    if (tdir == TimeDir::end) throw std::out_of_range("TimeDir& operator ++");
    tdir = TimeDir(static_cast<std::underlying_type<TimeDir>::type>(tdir) + 1);
    return tdir;
  }
  std::ostream& operator <<(std::ostream& ost, TimeDir const& tdir) {
    if(tdir == TimeDir::forwards)
      ost << "forwards";
    else if (tdir == TimeDir::backwards)
      ost << "backwards";
    else
      ost <<"Unknown";
    return ost;
  }
}
