#ifndef KinKal_StrawXing_hh
#define KinKal_StrawXing_hh
//
//  Describe the material effects of a kinematic trajectory crossing a straw
//  Used in the kinematic Kalman fit
//
#include "KinKal/DetectorXing.hh"
#include "KinKal/StrawMat.hh"
#include "KinKal/Line.hh"
#include "KinKal/PieceClosestApproach.hh"

namespace KinKal {
  template <class KTRAJ> class StrawXing : public DetectorXing<KTRAJ> {
    public:
      using PKTRAJ = ParticleTrajectory<KTRAJ>;
      using DXING = DetectorXing<KTRAJ>;
      using PTCA = PieceClosestApproach<KTRAJ,Line>;
      // construct from PTCA (for use with hits)
      StrawXing(PTCA const& tpoca, StrawMat const& smat) : DXING(tpoca.particleToca()) , smat_(smat), axis_(tpoca.sensorTraj()) {
	update(tpoca); }
      virtual ~StrawXing() {}
      // DetectorXing interface
      virtual void update(PKTRAJ const& pktraj,double precision) override;
      // specific interface: this xing is based on PTCA
      void update(PTCA const& tpoca);
      virtual void print(std::ostream& ost=std::cout,int detail=0) const override;
      // accessors
      StrawMat const& strawMat() const { return smat_; }
    private:
      StrawMat const& smat_;
      Line axis_; // straw axis, expressed as a timeline
  };

  template <class KTRAJ> void StrawXing<KTRAJ>::update(PTCA const& tpoca) {
    if(tpoca.usable()){
      DXING::mxings_.clear();
      smat_.findXings(tpoca.doca(),sqrt(tpoca.docaVar()),tpoca.dirDot(),DXING::mxings_);
      DXING::xtime_ = tpoca.particleToca();
    } else
      throw std::runtime_error("CA failure");
  }

  template <class KTRAJ> void StrawXing<KTRAJ>::update(PKTRAJ const& pktraj,double precision) {
    // use current xing time create a hint to the CA calculation: this speeds it up
    CAHint tphint(DXING::xtime_, DXING::xtime_);
    PTCA tpoca(pktraj,axis_,tphint,precision);
    update(tpoca);
  }

  template <class KTRAJ> void StrawXing<KTRAJ>::print(std::ostream& ost,int detail) const {
    ost <<"Straw Xing time " << this->crossingTime();
    if(detail > 0){
      for(auto const& mxing : this->matXings()){
	ost << " " << mxing.dmat_.name() << " pathLen " << mxing.plen_;
      }
    }
    if(detail > 1){
      ost << " Axis ";
      axis_.print(ost,0);
    }
    ost << std::endl;
  }

}
#endif
