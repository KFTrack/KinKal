#ifndef KinKal_AxialBFieldMap_hh
#define KinKal_AxialBFieldMap_hh
//
// BFieldMap implementation taking external model of the axial field and making simple extrapoloations to compute off-axis elements
//
#include "KinKal/General/BFieldMap.hh"
#include <vector>
namespace KinKal {
  class AxialBFieldMap : public BFieldMap {
    public:
      AxialBFieldMap(std::string const& file); // read from a file of Z positions and field values; must be evenly spaced
      AxialBFieldMap(double zmin, double zmax, std::vector<double> const& field); // field values are assumed to be evenly spaced
      VEC3 fieldVect(VEC3 const& position) const override;
      Grad fieldGrad(VEC3 const& position) const override;
      VEC3 fieldDeriv(VEC3 const& position, VEC3 const& velocity) const override;
      bool inRange(VEC3 const& position) const override { return position.Z() > zmin_ && position.Z() < zmax_; }
      double zMin() const { return zmin_; }
      double zMax() const { return zmax_; }
      void print(std::ostream& os=std::cout) const override {
        os << "Axial Bfield between " << zMin() << " and " << zMax() << " with " << axial_.size()
          << " field values from "  << axial_.front() << " to "  << axial_.back() << std::endl;
      }
      virtual ~AxialBFieldMap(){}
      // disallow copy and equivalence
      AxialBFieldMap(AxialBFieldMap const& ) = delete;
      AxialBFieldMap& operator =(AxialBFieldMap const& ) = delete;
      auto const& field() const { return axial_; }
    private:
      double zmin_, zmax_, zstep_; // z limits of the field
      std::vector<double> axial_; // axial field
      size_t lowBound(double zval) const;
      double zval(size_t index) const { return zmin_ + index*zstep_; }
  };
}
#endif
