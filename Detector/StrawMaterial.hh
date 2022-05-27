#ifndef KinKal_StrawMaterial_hh
#define KinKal_StrawMaterial_hh
//
//  description of a local segment of a straw, including a
//  mixture for the straw, the gas, and the wire,
//
#include "KinKal/MatEnv/DetMaterial.hh"
#include "KinKal/Detector/MaterialXing.hh"
#include "KinKal/Detector/StrawXingConfig.hh"
#include "KinKal/MatEnv/MatDBInfo.hh"
#include "KinKal/Trajectory/ClosestApproachData.hh"

namespace KinKal {
  class StrawMaterial {
    public:
      // explicit constructor from geometry and materials
      StrawMaterial(double srad, double thick, double wrad,
          const MatEnv::DetMaterial *wallmat, const MatEnv::DetMaterial *gasmat, const MatEnv::DetMaterial *wiremat) :
        srad_(srad), thick_(thick), wrad_(wrad), wallmat_(wallmat), gasmat_(gasmat), wiremat_(wiremat) {
          srad2_ = srad_*srad_;
        }
      // construct using default materials
      StrawMaterial(MatEnv::MatDBInfo const& matdbinfo,double srad, double thick, double wrad,
      const char* wallmat="straw-wall", const char* gasmat="straw-gas", const char* wiremat="straw-wire") :
        StrawMaterial(srad,thick,wrad, matdbinfo.findDetMaterial(wallmat),
            matdbinfo.findDetMaterial(gasmat),
            matdbinfo.findDetMaterial(wiremat)) {}
      // pathlength through straw components, given closest approach
      void pathLengths(ClosestApproachData const& cadata,StrawXingConfig const& caconfig, double& wallpath, double& gaspath, double& wirepath) const;
      // transit length given closest approach
      double transitLength(ClosestApproachData const& cadata) const;
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
      double thick_; // straw wall thickness
      double wrad_; // transverse radius of the wire
      const MatEnv::DetMaterial* wallmat_; // material of the straw wall
      const MatEnv::DetMaterial* gasmat_; // material of the straw gas
      const MatEnv::DetMaterial* wiremat_; // material of the wire
  };
}
#endif
