#include "KinKal/Geometry/Ray.hh"
std::ostream& operator <<(std::ostream& ost, KinKal::Ray const& ray){
  ost << "Ray starting " << ray.start_ << " direction " << ray.dir_;
  return ost;
}
