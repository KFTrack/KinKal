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
      LineSegment(VEC3 const& dir, VEC3 const& start, double length) : Ray(dir,start), length_(length) {
        end_ = position(length_);
      }
      // construct from 2 points
      LineSegment(VEC3 const& start, VEC3 const& end) : Ray(end-start,start), length_((end-start).R()), end_(end) {}
      // accessors
      VEC3 const& end() const { return end_; }
      double length() const { return length_; }
      void print(std::ostream& ost, int detail) const;
      // position from the end
      VEC3 endPosition(double dist) const { return end_ - dist*direction(); }
    private:
      double length_; // segment length
      VEC3 end_; // end position
  };
  std::ostream& operator <<(std::ostream& ost, LineSegment const& tLineSegment);
}
#endif
