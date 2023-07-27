//
//  Description  a geometric ray as a point and direction
//  original author: David Brown (LBN) 2023
//
#ifndef KinKal_Ray_hh
#define KinKal_Ray_hh
#include "KinKal/General/Vectors.hh"
#include <ostream>
namespace KinKal {
  struct Ray {
    // construct, making sure direciton is unit
    Ray(VEC3 const& dir, VEC3 const& start) : dir_(dir.Unit()), start_(start){}
    VEC3 dir_; // direction
    VEC3 start_; // starting position
    VEC3 position(double distance) const { return start_ + distance*dir_; }
  };
}
std::ostream& operator <<(std::ostream& ost, KinKal::Ray const& ray);
#endif
