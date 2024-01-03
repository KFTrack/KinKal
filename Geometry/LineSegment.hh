#ifndef KinKal_LineSegment_hh
#define KinKal_LineSegment_hh
//
// directed line segment
//
#include "KinKal/General/Vectors.hh"
#include "KinKal/Geometry/Ray.hh"

namespace KinKal {
  class LineSegment : public Ray {
    public:
      // construct from a point and direction
      LineSegment(VEC3 const& dir, VEC3 const& start, double length) : Ray(dir,start), length_(length) {}
      // construct from 2 points
      LineSegment(VEC3 const& start, VEC3 const& end) : Ray(end-start,start), length_((end-start).R()) {}
      // accessors
      VEC3 end() const { return position(length_); }
      double length() const { return length_; }
      void print(std::ostream& ost, int detail) const;
    private:
      double length_; // segment length
  };
  std::ostream& operator <<(std::ostream& ost, LineSegment const& tLineSegment);
}
#endif
