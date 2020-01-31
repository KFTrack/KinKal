#ifndef KinKal_PLine_hh
#define KinKal_PLine_hh
//
//  Linear time-based trajectory (approximately) perpendicular to the z axis
//  Models a sensor with constant signal propagation velocity
//  Used as part of the kinematic Kalman fit
//
#include "BTrk/KinKal/TTraj.hh"
#include "BTrk/KinKal/TPars.hh"
namespace KinKal {
  class PLine : public TTraj {
    public:
      enum paramIndex {d0_=0,phi0_=1,z0_=2,cost_=3,t0_=4,npars_=5};
      static size_t nParams() { return npars_; }
      static std::vector<std::string> const& paramNames(); 
      static std::vector<std::string> const& paramTitles();
      static std::string const& paramName(paramIndex index);
      static std::string const& paramTitle(paramIndex index);
 
      // construct from a wire position, signal propagation velocity (mm/ns), and measurement time at the position
      PLine(Vec3 const& p0, Vec3 const& svel, double tmeas);

    // named parameter accessors
      double param(size_t index) const { return pars_.vec_[index]; }
      double d0() const { return pars_.vec_[d0_]; }
      double phi0() const { return pars_.vec_[phi0_]; }
      double z0() const { return pars_.vec_[z0_]; }
      double cost() const { return pars_.vec_[cost_]; }
      double t0() const { return pars_.vec_[t0_]; }
    
      // simple functions 
      double velocity() const { return vel_; }
      double cosTheta() const { return cost(); }
      double sinTheta() const { return sqrt(1.0-cost()*cost()); }

      // position at t=t0
      void pos0(Vec3& pos) const;

      // geometric accessors
      void position(Vec4& pos) const override;
      void position(double time, Vec3& pos) const override;
      void velocity(double time, Vec3& vel) const override;
      void direction(double time, Vec3& dir) const override;

    private:
      TPars<npars_> pars_; // parameters
      double vel_; // signed linear velocity, translates time to distance along the trajectory (mm/nsec)

      static std::vector<std::string> paramTitles_;
      static std::vector<std::string> paramNames_;
  };
}
#endif

