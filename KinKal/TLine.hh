#ifndef KinKal_TLine_hh
#define KinKal_TLine_hh
//
//  Linear time-based trajectory with a constant velocity.
//  Used as part of the kinematic Kalman fit
//
#include "KinKal/Vectors.hh"
#include "KinKal/TRange.hh"
#include "KinKal/PData.hh"
namespace KinKal {
  class TLine {
    public:
      enum ParamIndex {d0_=0,phi0_=1,z0_=2,cost_=3,t0_=4,npars_=5};
      constexpr static size_t NParams() { return npars_; }
      typedef PData<npars_> PDATA; // Data payload for this class
      typedef typename PDATA::DVEC DVEC; // derivative of parameters type
      static std::vector<std::string> const& paramNames(); 
      static std::vector<std::string> const& paramUnits(); 
      static std::vector<std::string> const& paramTitles();
      static std::string const& paramName(ParamIndex index);
      static std::string const& paramTitle(ParamIndex index);
      static std::string const& paramUnit(ParamIndex index);

      // construct from a spacepoint and propagation velocity (mm/ns)
      // by default, the line has infinite unforced range
      TLine(Vec4 const& p0, Vec3 const& svel, TRange const& range=TRange(),bool forcerange=false);
      TLine(Vec3 const& p0, Vec3 const& svel, double tmeas, TRange const& range=TRange(),bool forcerange=false);
    
//      TLine(PDATA const& pdata) : pars_(pdata){std::cout<<" T Constructor 3 "<<pdata<<std::endl;  };
      TLine(PDATA::DVEC const &pvec, PDATA::DMAT const &pcov) : pars_(pvec, pcov){std::cout<<" T Constructor 4 "<<std::endl;  };
      // named parameter accessors
      double paramVal(size_t index) const { return pars_.parameters()[index]; }
      PDATA const &params() const { return pars_; }
      PDATA &params() { return pars_; }
      double d0() const { return paramVal(d0_); }
      double phi0() const { return paramVal(phi0_); }
      double z0() const { return paramVal(z0_); }
      double cost() const { return paramVal(cost_); }
      double t0() const { return paramVal(t0_); }
    
      // simple functions 
      double cosTheta() const { return cost(); }
      double sinTheta() const { return sqrt(1.0-cost()*cost()); }
      double cosPhi0() const { return cos(phi0()); }
      double sinPhi0() const { return sin(phi0()); }
      double theta() const { return acos(cost());} 
      double tanTheta() const { return sqrt(1.0-cost()*cost())/cost(); }

      // access position and direction
      Vec3 const& pos0() const { return pos0_; }
      double speed() const { 
         return speed_; 
      }

      // are we forcing the range?
      bool forceRange() const { return forcerange_; }

      // TOCA for a given point
      double TOCA(Vec3 point) const;

      // geometric accessors
      void position(Vec4& pos) const ;
      Vec3 position(double time) const ;
      Vec3 velocity(double time) const ;

      Vec3 const& direction(double time) const  {
        return dir_;
      };
      Vec3 const& dir() const { return dir_;}

      double speed(double time) const ;
      void setspeed(double speed) { 
        speed_ = speed;
      }

      void print(std::ostream& ost, int detail) const ;

      TRange const& range() const { return trange_; }
      TRange& range() { std::cout<<"Setting TRange "<<trange_<<std::endl; return trange_;  }
      virtual void setRange(TRange const& trange) { trange_ = trange; std::cout<<"Setting TRange from LH "<<trange_<<std::endl; }

      bool inRange(double time) const { return trange_.inRange(time); }

    protected:
      TRange trange_;
      PDATA pars_; // parameters
      double speed_; // signed linear velocity, translates time to distance along the trajectory (mm/nsec)
      Vec3 pos0_, dir_; // caches
      bool forcerange_; // if set, strictly enforce the range

      static std::vector<std::string> paramTitles_;
      static std::vector<std::string> paramNames_;
      static std::vector<std::string> paramUnits_;

     // nonconst accessors
      double& param(size_t index)  { return pars_.parameters()[index]; }

  };
   std::ostream& operator <<(std::ostream& ost, TLine const& line);
}
#endif
