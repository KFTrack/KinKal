#include "KinKal/Geometry/Ray.hh"
std::ostream& operator <<(std::ostream& ost, KinKal::Ray const& ray){
  ost << "Ray starting " << ray.start() << " direction " << ray.direction();
  return ost;
}
