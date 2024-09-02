//
//  Description  a geometric ray as a point and direction
//  original author: David Brown (LBN) 2023
//
#ifndef KinKal_Ray_hh
#define KinKal_Ray_hh
#include "KinKal/General/Vectors.hh"
#include <ostream>
namespace KinKal {
  class Ray {
    public:
      // construct, making sure direciton is unit
      Ray(VEC3 const& dir, VEC3 const& start) : dir_(dir.Unit()), start_(start){}
      VEC3 position(double distance) const { return start_ + distance*dir_; }
      VEC3 const& direction() const { return dir_; }
      VEC3 const& start() const { return start_; }
      // Distance along the the ray to the point of Closest Approach to a point
      double distanceTo(VEC3 const& point) const { return (point - start_).Dot(dir_); }
      void reverse() { dir_ *= -1.0; } // reverse direction
    private:
      VEC3 dir_; // direction
      VEC3 start_; // starting position
  };
}
std::ostream& operator <<(std::ostream& ost, KinKal::Ray const& ray);
#endif
