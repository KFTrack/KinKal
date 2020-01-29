#ifndef KinKal_PWire_hh
#define KinKal_PWire_hh
//
//  Linear trajectory representing a wire sensor (approximately) perpendicular to the z axis
//  Used as part of the kinematic Kalman fit
//
#include "BTrk/KinKal/TTraj.hh"
namespace KinKal {
  class PWire : public TTraj<5> {
    public:
      enum paramIndex {d0_=0,phi0_=1,z0_=2,cost_=3,t0_=4,npars_=5};
      static size_t nParams() { return npars_; }
      static std::vector<std::string> const& paramNames(); 
      static std::vector<std::string> const& paramTitles();
      static std::string const& paramName(paramIndex index);
      static std::string const& paramTitle(paramIndex index);
 
      // construct from parameters and covariance
      PWire(PVec const& pars, PMat const& pcov, double vel) : TTraj(pars,pcov),vel_(vel) {}
      // construct from a wire position, signal propagation velocity (mm/ns), and measurement time at the position
      PWire(Vec3 const& p0, Vec3 const& svel, double tmeas);
      // direct accessors
      double velocity() const { return vel_; }
      double cosTheta() const { return pars_[cost_]; }
      double sinTheta() const { return sqrt(1.0-pars_[cost_]*pars_[cost_]); }

      // named parameter accessors
      double d0() const { return pars_[d0_]; }
      double phi0() const { return pars_[phi0_]; }
      double z0() const { return pars_[z0_]; }
      double cost() const { return pars_[cost_]; }
      double t0() const { return pars_[t0_]; }
      // position at t=t0
      void pos0(Vec3& pos) const;

      // geometric accessors
      void position(Vec4& pos) const override;
      void position(double time, Vec3& pos) const override;
      void velocity(double time, Vec3& vel) const override;
      void direction(double time, Vec3& dir) const override;

    private:
      double vel_; // signed linear velocity, translates time to distance along the trajectory (mm/nsec)

      static std::vector<std::string> paramTitles_;
      static std::vector<std::string> paramNames_;
  };
}
#endif

