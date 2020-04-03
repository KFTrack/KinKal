#ifndef KinKal_StrawMat_hh
#define KinKal_StrawMat_hh
//
//  description of a local segment of a straw, including a
//  mixture for the straw, the gas, and the wire,
//
#include "MatEnv/DetMaterial.hh"
#include "KinKal/MatXing.hh"
namespace KinKal {
  class StrawMat {
    public:
    // explicit constructor from geometry and materials
      StrawMat(float srad, float thick, float wrad,
	  MatEnv::DetMaterial const& wallmat, MatEnv::DetMaterial const& gasmat, MatEnv::DetMaterial const& wiremat) :
	srad_(srad), thick_(thick), wrad_(wrad), wallmat_(wallmat), gasmat_(gasmat), wiremat_(wiremat) { 
	  srad2_ = srad_*srad_;
	  rdmax_ = (srad_ - thick_)/srad_;
	  wpmax_ = sqrt(8.0*srad_*thick_);
	  ddmax_ = 0.05*srad_;
	}
      // pathlength through gas, give DOCA to the axis, uncertainty on that,
      // and the dot product of the path direction WRT the axis.
      float gasPath(float doca, float ddoca, float adot) const;
      // same for wall material
      float wallPath(float doca, float ddoca, float adot) const; 
      // should add function to compute wire effect (probabilstically) FIXME!
      // find the material crossings given doca and error on doca.  Should allow for straw and wire to have different axes FIXME!
      void findXings(float doca, float ddoca, float adot, std::vector<MatXing>& mxings) const;
      float strawRadius() const { return srad_; }
      float wallThickness() const { return thick_; }
      float wireRadius() const { return wrad_; }
    private:
      float srad_; // outer transverse radius of the straw
      float srad2_; // outer transverse radius of the straw squared
      float rdmax_; // maximum relative DOCA
      float wpmax_; // maximum wall path
      float ddmax_; // max ddoca to integrate
      float thick_; // straw wall thickness
      float wrad_; // transverse radius of the wire
      MatEnv::DetMaterial const& wallmat_; // material of the straw wall
      MatEnv::DetMaterial const& gasmat_; // material of the straw gas
      MatEnv::DetMaterial const& wiremat_; // material of the wire
  };
}
#endif
