#ifndef KinKal_Context_hh
#define KinKal_Context_hh
//
// References to quantities supplied externally to KinKal
//
#include "BTrk/KinKal/BField.hh"
namespace KinKal {
  class Context {
    public:
      BField const& bField() const { return bfield_; }
      double bNom() const { return bz_; }
      Context(BField const& bfield) : bfield_(bfield) {
	Vec3 bnom; bfield.fieldVect(bnom); bz_ = bnom.Z(); }
    private:
      BField const& bfield_;
      double bz_; // nominal value of the longitudinal magentic field
  };

}
#endif
