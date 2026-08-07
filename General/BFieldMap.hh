#ifndef KinKal_BFieldMap_hh
#define KinKal_BFieldMap_hh
// class defining a BFieldMap Map interface for use in KinKal.
#include "KinKal/General/Vectors.hh"
#include "KinKal/General/TimeRange.hh"
#include "KinKal/General/TimeDir.hh"
#include "KinKal/General/PhysicalConstants.h"
#include "Math/SMatrix.h"
#include <vector>
#include <limits>
#include <algorithm>
#include <cstdarg>
#include <cmath>
#include <functional>
#include <ostream>

namespace KinKal {
  class BFieldMap {
    public:
      using Grad = ROOT::Math::SMatrix<double,3>; // field gradient: ie dBi/d(x,y,z)
      // return value of the field at a point
      virtual VEC3 fieldVect(VEC3 const& position) const = 0;
      // return BFieldMap gradient = dB_i/dx_j, at a given point
      virtual Grad fieldGrad(VEC3 const& position) const = 0;
      // return the BFieldMap derivative at a given point along a given velocity, WRT time
      virtual VEC3 fieldDeriv(VEC3 const& position, VEC3 const& velocity) const = 0;
      // is the point inside the range of this map?
      virtual bool inRange(VEC3 const& position) const = 0;
      virtual ~BFieldMap(){}
      virtual void print(std::ostream& os ) const = 0;
      // smallest |B| (T) usable for field-corrected transport; 0 (default) disables low-field
      // protection, leaving an unusable sample a hard failure for the caller
      BFieldMap(double minfield=0.0) : minfield_(minfield) {}
      double minField() const { return minfield_; }
      bool protecting() const { return minfield_ > 0.0; } // is low-field protection enabled?
      // Can this position support field-corrected transport? Answered from the region computed once at
      // construction, so this costs no field evaluation; maps registering no region use the direct test.
      bool usable(VEC3 const& position) const {
        if(!region_.built_) return directlyUsable(position);
        return region_.usable(position);
      }
      // pre-registered value returned for positions outside the usable region; no lookup needed
      VEC3 const& invalidField() const { return invalidfield_; }
      void printUsableRegion(std::ostream& os) const {
        if(!region_.built_){
          os << "BField usable region: none precomputed (direct test, minField " << minfield_ << " T)" << std::endl;
        } else {
          size_t ncell = region_.usable_.size();
          os << "BField usable region: minField " << minfield_ << " T, cylindrical r<" << region_.rmax_
            << " z[" << region_.zlow_ << "," << region_.zlow_+region_.zdim_ << "] mm, cell " << region_.cell_
            << " mm x " << region_.nphi_ << " phi bins, " << region_.nr_ << "x" << region_.nphi_ << "x"
            << region_.nz_ << " = " << ncell << " cells, " << region_.nsampled_ << " field samples, "
            << region_.nusable_ << " usable (" << (ncell>0 ? 100.0*region_.nusable_/ncell : 0.0) << "%)"
            << std::endl;
        }
      }
      // disallow copy and equivalence
      BFieldMap(BFieldMap const& ) = delete;
      BFieldMap& operator =(BFieldMap const& ) = delete;
      // speed of light in units to convert Tesla to mm (bending radius)
      static double constexpr cbar() { return CLHEP::c_light/1000.0; }
      // |B| below this (T) is treated as physically zero when deciding to hand off extrapolation;
      // fit paths must not sample here
      static double constexpr zeroField() { return 1.0e-6; }
      static bool isZeroField(VEC3 const& bvec) { return bvec.R() < zeroField(); }
      // templated interface for interacting with kinematic trajectory classes
      // how far can you go along the given kinematic trajectory till BField inhomogeneity makes the momentum accuracy out of (fractional) tolerance
      template<class KTRAJ> double rangeInTolerance(KTRAJ const& ktraj, double tstart, double tol) const;
      // the domain step at tstart: rangeInTolerance with a floor applied. Callers must have checked
      // usable(position) first -- rangeInTolerance samples fieldDeriv, which is only defined in range.
      template<class KTRAJ> double domainStep(KTRAJ const& ktraj, double tstart, double tol, double mindtstep) const {
        return std::max(rangeInTolerance(ktraj,tstart,tol),mindtstep);
      }
      // integrate the residual magentic force over the given kinematic trajectory and range due to the difference between the true field and the nominal field in the
      template<class KTRAJ> VEC3 integrate(KTRAJ const& ktraj, TimeRange const& trange) const;
    protected:
      // Map the usable region once, at construction, on a cylindrical lattice about the z axis; called by
      // subclasses once their extents are known. nphi need only follow the departure from axial symmetry.
      void buildUsableRegion(double rmax, double zlow, double zhigh, double cellsize, unsigned nphi);
      // True only when the precomputed region positively excludes this point, so a map can return its
      // pre-registered invalid value without a lookup. Must not call usable(): that would recurse.
      bool preRegisteredInvalid(VEC3 const& position) const {
        return region_.built_ && !region_.usable(position);
      }
      // the direct, per-sample test. Used to build the region, and as the fallback for maps that
      // register none.
      bool directlyUsable(VEC3 const& position) const {
        if(!inRange(position)) return false;             // outside: fieldDeriv is undefined here
        if(!protecting()) return true;                   // unprotected: in-range was the only requirement
        return fieldVect(position).R() >= std::max(minfield_,zeroField());
      }
    private:
      double minfield_; // smallest usable |B| (T); 0 disables low-field protection
      VEC3 invalidfield_ = VEC3(0.0,0.0,0.0); // pre-registered return outside the usable region
      // Precomputed occupancy of the usable region, in cylindrical coordinates about the z axis: for a
      // near-axially-symmetric field the boundary is close to r = rmax(phi,z), so azimuth needs few bins.
      struct UsableRegion {
        bool built_ = false;
        double rmax_ = 0.0, zlow_ = 0.0, zdim_ = 0.0;
        double cell_ = 0.0;      // radial and axial bin size
        double dphi_ = 0.0;      // azimuthal bin size
        unsigned nr_ = 0, nphi_ = 0, nz_ = 0;
        std::vector<bool> usable_;
        size_t nsampled_ = 0;    // field samples taken while building (the run-start cost)
        size_t nusable_ = 0;
        bool usable(VEC3 const& p) const {
          double z = p.Z()-zlow_;
          if(z < 0.0 || z >= zdim_) return false;
          double r = std::sqrt(p.X()*p.X()+p.Y()*p.Y());
          if(r >= rmax_) return false;
          double phi = std::atan2(p.Y(),p.X());
          if(phi < 0.0) phi += 2.0*M_PI;
          unsigned ir = static_cast<unsigned>(r/cell_);
          unsigned iz = static_cast<unsigned>(z/cell_);
          unsigned ip = static_cast<unsigned>(phi/dphi_);
          if(ip >= nphi_) ip = nphi_-1; // guard the wrap at exactly 2pi
          return usable_[(static_cast<size_t>(ir)*nphi_ + ip)*nz_ + iz];
        }
      };
      UsableRegion region_;
  };

  // Corners are sampled once each and a cell is marked usable only if all eight are, so the recorded
  // region is conservative: it never claims space the field cannot support.
  inline void BFieldMap::buildUsableRegion(double rmax, double zlow, double zhigh, double cellsize, unsigned nphi) {
    auto& r = region_;
    // a non-positive cell (or radius) cannot describe a region: leave the direct test in place rather
    // than building a degenerate one that would report everything unusable
    if(cellsize <= 0.0 || rmax <= 0.0 || zhigh <= zlow) return;
    r.rmax_ = rmax; r.zlow_ = zlow; r.zdim_ = zhigh-zlow; r.cell_ = cellsize;
    r.nr_ = std::max(1u,static_cast<unsigned>(std::ceil(rmax/cellsize)));
    r.nz_ = std::max(1u,static_cast<unsigned>(std::ceil(r.zdim_/cellsize)));
    r.nphi_ = std::max(1u,nphi);
    r.dphi_ = 2.0*M_PI/r.nphi_;
    r.rmax_ = r.nr_*cellsize; r.zdim_ = r.nz_*cellsize;
    // corner lattice: (nr+1) x nphi x (nz+1), azimuth periodic so no extra plane
    unsigned mr = r.nr_+1, mz = r.nz_+1;
    std::vector<bool> node(static_cast<size_t>(mr)*r.nphi_*mz,false);
    for(unsigned ir=0; ir<mr; ++ir){
      double rad = ir*cellsize;
      for(unsigned ip=0; ip<r.nphi_; ++ip){
        double phi = ip*r.dphi_;
        double cphi = std::cos(phi), sphi = std::sin(phi);
        for(unsigned iz=0; iz<mz; ++iz){
          VEC3 pos(rad*cphi,rad*sphi,zlow+iz*cellsize);
          node[(static_cast<size_t>(ir)*r.nphi_ + ip)*mz + iz] = directlyUsable(pos);
        }
      }
    }
    r.nsampled_ = node.size();
    r.usable_.assign(static_cast<size_t>(r.nr_)*r.nphi_*r.nz_,false);
    for(unsigned ir=0; ir<r.nr_; ++ir){
      for(unsigned ip=0; ip<r.nphi_; ++ip){
        unsigned ipn = (ip+1)%r.nphi_; // azimuth wraps
        for(unsigned iz=0; iz<r.nz_; ++iz){
          bool ok = true;
          for(unsigned dr=0; dr<2 && ok; ++dr)
            for(unsigned dp=0; dp<2 && ok; ++dp)
              for(unsigned dz=0; dz<2 && ok; ++dz)
                ok = node[(static_cast<size_t>(ir+dr)*r.nphi_ + (dp ? ipn : ip))*mz + (iz+dz)];
          r.usable_[(static_cast<size_t>(ir)*r.nphi_ + ip)*r.nz_ + iz] = ok;
        }
      }
    }
    r.nusable_ = std::count(r.usable_.begin(),r.usable_.end(),true);
    r.built_ = true;
  }

  template<class KTRAJ> VEC3 BFieldMap::integrate(KTRAJ const& ktraj, TimeRange const& trange) const {
    // take a fixed number of steps.  This may fail for long ranges FIXME!
    unsigned nsteps(10);
    double dt = trange.range()/nsteps;
    // now integrate
    VEC3 dmom;
    for(unsigned istep=0; istep< nsteps; istep++){
      double tstep = trange.begin() + (0.5+istep)*dt;
      VEC3 vel = ktraj.velocity(tstep);
      VEC3 db = fieldVect(ktraj.position3(tstep)) - ktraj.bnom(tstep);
      dmom += cbar()*ktraj.charge()*dt*vel.Cross(db);
    }
    return dmom;
  }

  // estimate how long the momentum vector from the given trajectory will stay within the given (fractional) tolerance given the field spatial variation
  // ie mag(P_true(tstart+dt) - P_traj(tstart+dt)) < tol.  This is a 2nd order local calculation
  template<class KTRAJ> double BFieldMap::rangeInTolerance(KTRAJ const& ktraj, double tstart, double tol) const {
    auto tpos = ktraj.position3(tstart); // starting position
    double dp = ktraj.momentum(tstart)*tol; // fractional tolerance on momentum
    auto vel = ktraj.velocity(tstart); // starting velocity
    auto dBdt = fieldDeriv(tpos,vel); // change in field WRT time along this velocity
    double d2pdt2 = (dBdt.Cross(vel)).R()*cbar()*fabs(ktraj.charge()); // 2nd derivative of momentum due to B change along the path
    if(d2pdt2 > 1e-10)
      return sqrt(dp/d2pdt2);
    else
      return ktraj.range().range();
  }

  // trivial instance of the above, used for testing
  class UniformBFieldMap : public BFieldMap {
    public:
      VEC3 fieldVect(VEC3 const& position) const override { return fvec_; }
      Grad fieldGrad(VEC3 const& position) const override { return Grad(); }
      VEC3 fieldDeriv(VEC3 const& position, VEC3 const& velocity) const override { return VEC3(); }
      void print(std::ostream& os =std::cout) const override { os << "Uniform BField, B = " << fvec_ << std::endl; }
      UniformBFieldMap(VEC3 const& bnom) : fvec_(bnom) {}
      UniformBFieldMap(double BZ) : UniformBFieldMap(VEC3(0.0,0.0,BZ)) {}
      bool inRange(VEC3 const& position) const override { return true; };
      virtual ~UniformBFieldMap(){}
      // disallow copy and equivalence
      UniformBFieldMap(UniformBFieldMap const& ) = delete;
      UniformBFieldMap& operator =(UniformBFieldMap const& ) = delete;
    private:
      VEC3 fvec_; // constant field
  };

  // use superposition to create a composite field
  class CompositeBFieldMap : public BFieldMap {
    public:
      using FCOL = std::vector<const BFieldMap*>;
      VEC3 fieldVect(VEC3 const& position) const override;
      Grad fieldGrad(VEC3 const& position) const override;
      VEC3 fieldDeriv(VEC3 const& position, VEC3 const& velocity) const override;
      void print(std::ostream& os =std::cout) const override;
      bool inRange(VEC3 const& position) const override;
      CompositeBFieldMap () {}
      CompositeBFieldMap(FCOL const& fields) : fields_(fields) {}
      void addField(BFieldMap const& field) { fields_.push_back(&field); }
      virtual ~CompositeBFieldMap() {}
      // disallow copy and equivalence
      CompositeBFieldMap(CompositeBFieldMap const& ) = delete;
      CompositeBFieldMap& operator =(CompositeBFieldMap const& ) = delete;

    private:
      FCOL fields_; // fields
  };

  // simple Z gradient field, used to test Field corrections
  class GradientBFieldMap : public BFieldMap {
    public:
      GradientBFieldMap(double b0, double b1, double zg0, double zg1);
      VEC3 fieldVect(VEC3 const& position) const override;
      Grad fieldGrad(VEC3 const& position) const override { return fgrad_; }
      VEC3 fieldDeriv(VEC3 const& position, VEC3 const& velocity) const override;
      bool inRange(VEC3 const& position) const override { return true; }
      void print(std::ostream& os =std::cout) const override { os << "BField with  constant gradient of " << grad_ << " Tesla/mm" << std::endl; }
      double gradient() const { return grad_; }
      virtual ~GradientBFieldMap(){}
      // disallow copy and equivalence
      GradientBFieldMap(GradientBFieldMap const& ) = delete;
      GradientBFieldMap& operator =(GradientBFieldMap const& ) = delete;
    private:
      double b0_, b1_;
      double z0_;
      double grad_; // gradient in tesla/mm, computed from the fvec values
      Grad fgrad_;
  };

}
#endif
