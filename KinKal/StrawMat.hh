#ifndef KinKal_StrawMat_hh
#define KinKal_StrawMat_hh
//
//  description of a local segment of a straw, including a
//  mixture for the straw, the gas, and the wire, allowing for
//  offset between the wire and the straw
//
#include "KinKal/DMat.hh"
#include "KinKal/MIsect.hh"
#include "MatEnv/DetMaterial.hh"
#include "KinKal/Vectors.hh"
namespace KinKal {
  class StrawMat : public DMat {
    public:
      StrawMat(float rad, float thick, float rwire,
	  DetMaterial const& wallmat, DetMaterial const& gasmat, DetMaterial const& wiremat) :
	rad_(rad), thick_(thick), rwire_(rwire), wallmat_(wallmat), gasmat_(gasmat), wiremat_(wiremat) { 
	  rad2_ = rad_*rad_;
	  rdmax_ = (rad_ - thick_)/rad_;
	}
      // pathlength through gas, give DOCA to the axis, uncertainty on that,
      // and the dot product of the path direction WRT the axis.
      float gasPath(float doca, float ddoca, float adot) const;
      // same for wall material
      float wallPath(float doca, float ddoca, float adot) const; 
      // add up the material effects for a given doca.
      void intersect(float doca, float ddoca, float adot, std::vector<MIsect>& misects) const;
      float strawRadius() const { return rad_; }
      float wallThickness() const { return thick_; }
      float wireRadius() const { return rwire_; }
    private:
      float rad_; // outer transverse radius of the straw
      float rad2_; // outer transverse radius of the straw squared
      float rdmax_; // maximum relative DOCA
      float thick_; // wall thickness
      float rwire_; // transverse radius of the wire
      Pol2 axis_; // cylinder axis offset WRT wire.  Not used FIXME!
      DetMaterial const& wallmat_; // material of the straw wall
      DetMaterial const& gasmat_; // material of the straw gas
      DetMaterial const& wiremat_; // material of the wire
  };
}
#endif
