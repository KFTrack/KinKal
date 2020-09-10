#ifndef KinKal_StrawHit_hh
#define KinKal_StrawHit_hh
//
//  class representing a straw sensor measurement.  It assumes a (possibly displaced)
//  circular outer cathode locally parallel to the wire.  All the work is done in the WireHit parent.
//  Used as part of the kinematic Kalman fit
//
#include "KinKal/WireHit.hh"
#include "KinKal/StrawXing.hh"
#include <memory>
namespace KinKal {
  template <class KTRAJ> class StrawHit : public WireHit<KTRAJ> {
    public:
      using WHIT = WireHit<KTRAJ>;
      using STRAWXING = StrawXing<KTRAJ>;
      using STRAWXINGPTR = std::shared_ptr<STRAWXING>;
      StrawHit(BFieldMap const& bfield, Line const& straj, WireCell const& cell, STRAWXINGPTR const& sxing,LRAmbig ambig=LRAmbig::null) :
	WHIT(sxing, bfield,straj,cell,ambig) {}
      void print(std::ostream& ost=std::cout,int detail=0) const override;
      virtual ~StrawHit(){}
    private:
      // add state for longitudinal resolution, transverse resolution; could be a 2ndary measurement TODO
  };

  template<class KTRAJ> void StrawHit<KTRAJ>::print(std::ostream& ost, int detail) const {
    if(this->isActive())
      ost<<"Active ";
    else
      ost<<"Inactive ";
    ost << " StrawHit LRAmbig " << this-> ambig() << " " << std::endl;
  }
}
#endif
