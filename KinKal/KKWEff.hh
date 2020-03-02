#ifndef KinKal_KKWEff_hh
#define KinKal_KKWEff_hh
//
// Class to add information; a measurement or constraint
// This effect provides information content and is processed in weight space 
//
#include "KinKal/KKEff.hh"
namespace KinKal {
  template<class KTRAJ> class KKWEff : public KKEff<KTRAJ> {
    public:
      typedef KKEff<KTRAJ> KKEFF;
      typedef PKTraj<KTRAJ> PKTRAJ;
      typedef typename KKEFF::WDATA WDATA; // forward the typedef
      typedef typename KKEFF::KKDATA KKDATA; // forward the typedef
      // process this effect given the adjacent effect
      virtual bool process(KKDATA& kkdata,TDir tdir) override;
      virtual bool append(PKTRAJ& fit) override {return true;} // wdatas don't change the output fit
      // construct from a wdata
      KKWEff(KTRAJ const& reftraj, WDATA const& wdata) : KKEFF(reftraj), wdata_(wdata) {}
      virtual ~KKWEff(){}
      WDATA const& wData() const { return wdata_; }
    protected:
      // allow subclasses to construct without a wdata
      KKWEff(KTRAJ const& reftraj) : KKEFF(reftraj) {}
      KKWEff() {} // no local reftraj; must be set
      WDATA wdata_; // wdata representation of this effect's constraint/measurement
  };

  template<class KTRAJ> bool KKWEff<KTRAJ>::process(KKDATA& kkdata,TDir tdir) {
    // direction is irrelevant for adding information
    bool retval(false);
    if(KKEffBase::status(tdir) == KKEffBase::processed){
      retval = true;
    } else if(this->isActive()){
      // add this effect's information
      kkdata.append(wdata_);
      retval = kkdata.wData().matrixOK();
    }
    KKEffBase::setStatus(tdir,KKEffBase::processed);
    return retval;
  }

}
#endif
