#ifndef KinKal_CylBFieldMap_hh
#define KinKal_CylBFieldMap_hh

#include "KinKal/General/CylBMap.hh"
#include "KinKal/General/InterpBilinear.hh"
#include "KinKal/General/BFieldMap.hh"

namespace KinKal {
  class CylBFieldMap : public BFieldMap {
    public:
      using VEC = std::vector<double>;
      using MAT = std::vector<std::vector<double>>;
      CylBFieldMap(std::string const& file); // supply data file
      VEC3 fieldVect(VEC3 const& position) const override;
      Grad fieldGrad(VEC3 const& position) const override;
      VEC3 fieldDeriv(VEC3 const& position, VEC3 const& velocity) const override;
      bool inRange(VEC3 const& position) const override;
      // TO DO!
      void print(std::ostream& os=std::cout) const override {
        os << "Cylindrically symmetric Bfield with boundaries z=[" << zMin() << ", " << zMax() << "] and r=["
          << rMin() << ", " << rMax() << "] with " << B_dat.ntot_
          << " field values from (lower left, upper left): Br=(" << B_dat.Br_.front().front() << ", "
          << B_dat.Br_.back().back() << ") Tesla and Bz=(" << B_dat.Bz_.front().front() << ", "
          << B_dat.Bz_.back().back() << ") Tesla" << std::endl;
      }
      virtual ~CylBFieldMap(){}
      // disallow copy and equivalence
      CylBFieldMap(CylBFieldMap const& ) = delete;
      CylBFieldMap& operator =(CylBFieldMap const& ) = delete;
      // getters
      double zMax() const { return B_dat.z_.back(); }
      double zMin() const { return B_dat.z_.front(); }
      double rMax() const { return B_dat.r_.back(); }
      double rMin() const { return B_dat.r_.front(); }
    private:
      const CylBMap B_dat;
      const InterpBilinear Br_interp, Bz_interp;
  };
}
#endif
