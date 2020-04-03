#ifndef KinKal_StrawHit_hh
#define KinKal_StrawHit_hh
//
//  class representing a straw sensor measurement.  It assumes a (possibly displaced)
//  circular outer cathode locally parallel to the wire.
//  Used as part of the kinematic Kalman fit
//
#include "KinKal/WireHit.hh"
#include "KinKal/StrawXing.hh"
namespace KinKal {
  template <class KTRAJ> class StrawHit : public WireHit<KTRAJ> {
    public:
      typedef PKTraj<KTRAJ> PKTRAJ;
      typedef WireHit<KTRAJ> WHIT;
      typedef StrawXing<KTRAJ> STRAWXING;
      StrawHit(BField const& bfield, PKTRAJ const& pktraj, TLine const& straj, D2T const& d2t, STRAWXING const& sxing,float nulldoca, LRAmbig ambig=LRAmbig::null,bool active=true) :
	WireHit<KTRAJ>(bfield,pktraj,straj,d2t,std::min(nulldoca,sxing.strawMat().strawRadius())*std::min(nulldoca,sxing.strawMat().strawRadius())/3.0,ambig,active), sxing_(sxing) {}
      virtual float tension() const override { return 0.0; } // check against straw diameter, length, any other measurement content FIXME!
      virtual STRAWXING* detCrossing() override { return &sxing_; }
      virtual ~StrawHit(){}
    private:
      STRAWXING sxing_; // straw material crossing information
      // add state for longitudinal resolution, transverse resolution FIXME!
  };
}
#endif
