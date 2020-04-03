#ifndef KinKal_StrawXing_hh
#define KinKal_StrawXing_hh
//
//  Describe the material effects of a kinematic trajectory crossing a straw
//  Used in the kinematic Kalman fit
//
#include "KinKal/DXing.hh"
#include "KinKal/StrawMat.hh"
#include "KinKal/TPoca.hh"

namespace KinKal {
  template <class KTRAJ> class StrawXing : public DXing<KTRAJ> {
    public:
      typedef PKTraj<KTRAJ> PKTRAJ;
      typedef DXing<KTRAJ> DXING;
      typedef TPoca<PKTRAJ,TLine> TPOCA;

      // construct from a trajectory and a time:
      StrawXing(PKTRAJ const& pktraj,double xtime, StrawMat const& smat, TLine const& axis) : DXING(pktraj,xtime), smat_(smat), axis_(axis) {
	update(pktraj); } 
      // construct from TPOCA (for use with hits)
      StrawXing(TPOCA const& tpoca, StrawMat const& smat) : DXING(tpoca) , smat_(smat), axis_(tpoca.sensorTraj()) {
	update(tpoca); }
      virtual ~StrawXing() {}
      // DXing interface
      virtual void update(PKTRAJ const& pktraj) override;
      virtual void update(TPocaBase const& tpoca) override;
      // accessors
      StrawMat const& strawMat() const { return smat_; }

    private:
      StrawMat const& smat_;
      TLine axis_; // straw axis, expressed as a timeline
  };

  template <class KTRAJ> void StrawXing<KTRAJ>::update(TPocaBase const& tpoca) {
    DXING::mxings_.clear();
    smat_.findXings(tpoca.doca(),tpoca.dDoca(),tpoca.dirDot(),DXING::mxings_);
    DXING::xtime_ = tpoca.particleToca();
  }

  template <class KTRAJ> void StrawXing<KTRAJ>::update(PKTRAJ const& pktraj) {
    TPOCA tpoca(pktraj,axis_);
    update(tpoca);
  }
}
#endif
