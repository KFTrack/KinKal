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
      typedef PKTraj<KTRAJ> PKTRAJ;
      typedef WireHit<KTRAJ> WHIT;
      typedef THit<KTRAJ> THIT;
      typedef StrawXing<KTRAJ> STRAWXING;
      typedef std::shared_ptr<STRAWXING> STRAWXINGPTR;
      StrawHit(BField const& bfield, TLine const& straj, D2T const& d2t, STRAWXINGPTR const& sxing,LRAmbig ambig=LRAmbig::null) :
	WHIT(sxing, bfield,straj,d2t,sxing->strawMat().strawRadius(),ambig) {}
      virtual double tension() const override { return 0.0; } // check against straw diameter, length, any other measurement content FIXME!
      virtual void print(std::ostream& ost=std::cout,int detail=0) const override;
      virtual ~StrawHit(){}
    private:
      // add state for longitudinal resolution, transverse resolution, to use in tension measurement FIXME!
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
