#ifndef KinKal_DStraw_hh
#define KinKal_DStraw_hh
//
//  description of a local segment of a straw, including a
//  mixture for the straw, the gas, and the wire, allowing for
//  offset between the wire and the straw
//
#include "KinKal/StrawMat.hh"
#include "KinKal/DPiece.hh"
#include "KinKal/TLine.hh"
namespace KinKal {
  template <class TT> class DStraw : public DPiece<TT> {
    public:
    // explicit constructor from geometry and materials
      DStraw(StrawMat const& smat, TLine const& tline) : smat_(smat), tline_(tline) {}
      // DPiece interface; first, for pieces associated with a hit
      virtual void findXings(TPocaBase const& poca,std::vector<MatXing>& mxings) const override;
      // also for materials unassociated with a hit
      virtual void findXings(TT const& ttraj, MatXingCol& mxings) const =0;

      // pathlength through gas, give DOCA to the axis, uncertainty on that,
      // and the dot product of the path direction WRT the axis.
      float gasPath(float doca, float ddoca, float adot) const;
      // same for wall material
      float wallPath(float doca, float ddoca, float adot) const;  // should add doca to the wire for wire material effects FIXME!
      // find the material crossings given doca.
      void findXings(float doca, float ddoca, float adot, std::vector<MatXing>& mxings) const;
      float strawRadius() const { return rad_; }
      float wallThickness() const { return thick_; }
      float wireRadius() const { return rwire_; }
    private:
      float rad_; // outer transverse radius of the straw
      float rad2_; // outer transverse radius of the straw squared
      float rdmax_; // maximum relative DOCA
      float wpmax_; // maximum wall path
      float ddmax_; // max ddoca to integrate
      float thick_; // wall thickness
      float rwire_; // transverse radius of the wire
      MatEnv::DetMaterial const& wallmat_; // material of the straw wall
      MatEnv::DetMaterial const& gasmat_; // material of the straw gas
      MatEnv::DetMaterial const& wiremat_; // material of the wire
  };
}
#endif
