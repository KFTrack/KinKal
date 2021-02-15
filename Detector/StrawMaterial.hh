#ifndef KinKal_StrawMaterial_hh
#define KinKal_StrawMaterial_hh
//
//  description of a local segment of a straw, including a
//  mixture for the straw, the gas, and the wire,
//
#include "KinKal/MatEnv/DetMaterial.hh"
#include "KinKal/Detector/MaterialXing.hh"
#include "KinKal/MatEnv/MatDBInfo.hh"
#include "KinKal/Trajectory/ClosestApproachData.hh"

namespace KinKal {
// simple struct to hold crossing calculation configuration parameters
  struct StrawXingConfig {
    double minsigdoca_; // minimum doca sigma to integrate
    double nsig_; // number of sigma past wall to consider 'inside' the straw
    StrawXingConfig(double ddmax, double nsig) : minsigdoca_(ddmax), nsig_(nsig) {}
  };

  class StrawMaterial {
    public:
      // explicit constructor from geometry and materials
      StrawMaterial(double srad, double thick, double wrad,
	  const MatEnv::DetMaterial *wallmat, const MatEnv::DetMaterial *gasmat, const MatEnv::DetMaterial *wiremat) :
	srad_(srad), thick_(thick), wrad_(wrad), wallmat_(wallmat), gasmat_(gasmat), wiremat_(wiremat) { 
	  srad2_ = srad_*srad_;
	  wpmax_ = sqrt(8.0*srad_*thick_);
	}
      // construct using default materials
      StrawMaterial(MatEnv::MatDBInfo const& matdbinfo,double srad, double thick, double wrad) :
	StrawMaterial(srad,thick,wrad, matdbinfo.findDetMaterial("straw-wall"),
	    matdbinfo.findDetMaterial("straw-gas"),
	    matdbinfo.findDetMaterial("straw-wire")) {}
      // pathlength through gas, give DOCA to the axis, uncertainty on that,
      // and the dot product of the path direction WRT the axis.
      double gasPath(ClosestApproachData const& cadata,StrawXingConfig const& caconfig) const;
      // same for wall material
      double wallPath(ClosestApproachData const& cadata,StrawXingConfig const& caconfig) const; 
      double wirePath(ClosestApproachData const& tpdata,StrawXingConfig const& config) const;
      // find the material crossings given doca and error on doca.  Should allow for straw and wire to have different axes TODO
      void findXings(ClosestApproachData const& cadata,StrawXingConfig const& caconfig, std::vector<MaterialXing>& mxings) const;
      double strawRadius() const { return srad_; }
      double wallThickness() const { return thick_; }
      double wireRadius() const { return wrad_; }
      MatEnv::DetMaterial const& wallMaterial() const { return *wallmat_; }
      MatEnv::DetMaterial const& gasMaterial() const { return *gasmat_; }
      MatEnv::DetMaterial const& wireMaterial() const { return *wiremat_; }
    private:
      double srad_; // outer transverse radius of the straw
      double srad2_; // outer transverse radius of the straw squared
      double wpmax_; // maximum wall path
      double thick_; // straw wall thickness
      double wrad_; // transverse radius of the wire
      const MatEnv::DetMaterial* wallmat_; // material of the straw wall
      const MatEnv::DetMaterial* gasmat_; // material of the straw gas
      const MatEnv::DetMaterial* wiremat_; // material of the wire
  };
}
#endif
