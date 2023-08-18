//
//  Implementation of a 'thin' solid as a surface with a given thickness
//
#ifndef KinKal_ThinSolid_hh
#define KinKal_ThinSolid_hh
#include "KinKal/Geometry/Surface.hh"
#include "KinKal/Geometry/Solid.hh"
#include <memory>
namespace KinKal {
  class ThinSolid : public Solid {
    public:
      virtual ~ThinSolid() {}
      // construct from a pointer to a surface.  Note this class TAKES OWNERSHIP of that object
      ThinSolid(std::unique_ptr<Surface> surf, double halfthick) : surf_(std::move(surf)), halfthick_(halfthick) {}
      // Solid interface
      bool isInside(VEC3 const& point) const override;
      // thin solid specific interface
      auto const& surface() const { return surf_; }
    private:
      std::unique_ptr<Surface> surf_; // surface pointer
      double halfthick_; // half-thickness, perpendicular to the surface
  };
}
#endif
