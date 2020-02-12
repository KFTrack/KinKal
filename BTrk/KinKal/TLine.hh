#ifndef KinKal_TLine_hh
#define KinKal_TLine_hh
//
//  Linear time-based trajectory with a constant velocity 
//  Used as part of the kinematic Kalman fit
//
#include "BTrk/KinKal/TTraj.hh"
#include "BTrk/KinKal/TData.hh"
namespace KinKal {
  class TLine : public TTraj {
    public:
      enum paramIndex {d0_=0,phi0_=1,z0_=2,cost_=3,t0_=4,npars_=5};
      constexpr static size_t NParams() { return npars_; }
      static std::vector<std::string> const& paramNames(); 
      static std::vector<std::string> const& paramTitles();
      static std::string const& paramName(paramIndex index);
      static std::string const& paramTitle(paramIndex index);
 
      // construct from a wire position, signal propagation velocity (mm/ns), and measurement time at the position
      TLine(Vec3 const& p0, Vec3 const& svel, double tmeas, TRange const& range=TRange());

    // named parameter accessors
      double param(size_t index) const { return pars_.vec()[index]; }
      double d0() const { return pars_.vec()[d0_]; }
      double phi0() const { return pars_.vec()[phi0_]; }
      double z0() const { return pars_.vec()[z0_]; }
      double cost() const { return pars_.vec()[cost_]; }
      double t0() const { return pars_.vec()[t0_]; }
    
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
      TData<npars_> pars_; // parameters
      double vel_; // signed linear velocity, translates time to distance along the trajectory (mm/nsec)

      static std::vector<std::string> paramTitles_;
      static std::vector<std::string> paramNames_;
  };
}
#endif

