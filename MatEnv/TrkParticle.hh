// Class to define particle types which can be reconstructed as tracks.
// the type values are copied from the PDG codes.
#ifndef TrkParticle_hh
#define TrkParticle_hh
// for now just the 5 stable 'fundamental' particles
//  plus deuterons, tritons, alphas.  It should include nuclear fragments, Cascade, someday
// types are a subset of PDGCodes
#include <string>
#include <map>
class TrkParticle {
  public:
    enum Type {
      e_minus = 11 ,
      e_plus = -11 ,
      mu_minus = 13 ,
      mu_plus = -13 ,
      pi_plus = 211 ,
      pi_minus = -211 ,
      K_plus = 321 ,
      K_minus = -321 ,
      p_plus = 2212 ,
      anti_p_minus = -2212
	//	deuterium = 1000010020,
	//	tritium = 1000010030,
	//	He3 = 1000020030,
	//	He4 = 100002004
    };
    // particle masses: good enough for tracking
    static constexpr double e_mass_ = 5.10998910E-01; // electron mass in MeVC^2
    static constexpr double mu_mass_ = 1.05658367E+02; 
    static constexpr double pi_mass_ = 1.3957018E+02; 
    static constexpr double K_mass_ = 4.93677E+02; 
    static constexpr double p_mass_ = 9.3827203E+02;
    // construct from a type
    TrkParticle(Type ptype=e_minus);
    ~TrkParticle(){}
    // basic accessor
    Type particleType() const { return type_; }
    bool operator == (TrkParticle const& other) const { return type_ == other.type_; }
    bool operator != (TrkParticle const& other) const { return ! this->operator==(other); }
    // return basic information
    double mass() const { return mass_; }
    int charge() const;
    std::string const& name() const;
    // basic kinematics; provide momentum in CLHEP units
    double energy(double momentum) const; // return value in CLHEP units
    double beta(double momentum) const;
    double betagamma(double momentum) const;
    double gamma(double momentum) const;
    static  double mass(TrkParticle::Type type);
  private:
    Type type_;
    double mass_, m2_; // particle mass cache
};

#endif
