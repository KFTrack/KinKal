#ifndef KinKal_StrawHit_hh
#define KinKal_StrawHit_hh
//
//  class representing a straw sensor measurement.  It assumes a (possibly displaced)
//  circular outer cathode locally parallel to the wire.
//  Used as part of the kinematic Kalman fit
//
#include "KinKal/WireHit.hh"
#include "KinKal/StrawMat.hh"
namespace KinKal {
  class StrawHit : public WireHit {
    public:
      StrawHit(TLine const& straj, BField const& bfield, D2T const& d2t,StrawMat const& smat,float nulldoca, LRAmbig ambig=null,bool active=true) : 
	WireHit(straj,bfield,d2t,std::min(nulldoca,smat.strawRadius())*std::min(nulldoca,smat.strawRadius())/3.0,ambig,active), smat_(smat) {}
      virtual float inRange(TPocaBase const& tpoca) const override;
      virtual void update(TPocaBase const& tpoca) const override;
      virtual ~StrawHit(){}
      virtual const StrawMat* material() const override { return &smat_; }
    private:
      StrawMat const& smat_; // straw material information
      // add state for longitudinal resolution, transverse resolution FIXME!
  };
}
#endif
