//
//  Abstract base class for a solid
//  These are used to describe material objects
//  original author: David Brown (LBN) 2023
//
#ifndef KinKal_Solid_hh
#define KinKal_Solid_hh
namespace KinKal {
  class Solid {
    public:
      virtual ~Solid() {}
      // determin if a point is inside the solid or not
      virtual bool isInside(VEC3 const& point) const = 0;
  };
}
#endif
