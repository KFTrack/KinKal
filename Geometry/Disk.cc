#include "KinKal/Geometry/Disk.hh"

std::ostream& operator <<(std::ostream& ost, KinKal::Disk const& disk) {
  KinKal::Plane const& plane = static_cast<KinKal::Plane const&>(disk);
  ost << "Disk in " << plane << " radius " << disk.outerRadius();
  return ost;
}
