#ifndef KinKal_WireCell_hh
#define KinKal_WireCell_hh
//
// base class describing a wire cell.  Eventually this could have more structure
//
#include "KinKal/General/Vectors.hh"
namespace KinKal {
  class WireCell {
    public:
      virtual ~WireCell() {}
      WireCell() {}
      // disallow copy and equivalence
      WireCell(WireCell const& ) = delete; 
      WireCell& operator =(WireCell const& ) = delete;
      virtual double size() const = 0; // approximate size perpendicular to the wire
      // given a drift DOCA and direction in the cell, compute most probable drift time, expected RMS of drift time, and the local drift speed
      virtual void distanceToTime(POL2 const& drift, double& tdrift, double& tdriftvar, double& dspeed) const = 0;
      virtual double averageDriftSpeed() const = 0; // average drift speed
      virtual double maximumDriftTime() const = 0; // maximum drift
      // add functions to describe material TODO
  };
}
#endif
