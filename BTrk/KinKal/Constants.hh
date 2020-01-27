#ifndef KinKal_Constants_hh
#define KinKal_Constants_hh
//
// values that don't change.  These should come from an authoritative source eventually
//
namespace KinKal {
  static constexpr double c_ = 299.792; // speed of light in mm/nsec
  // signed mass in mm from mass in MeV and Bz in tesla
  constexpr double reducedMass(double mass,double Bz,int charge) {
    return mass*1000.0/(charge*Bz*c_); }
}
#endif
