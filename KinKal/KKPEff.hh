#ifndef KinKal_KKPEff_hh
#define KinKal_KKPEff_hh
//
// Class to describe parameter transport and noise due to material effects or model descrepencies
// This effect provides information content and is processed in params space 
//
#include "KinKal/KKEff.hh"
namespace KinKal {
  template<class KTRAJ> class KKPEff : public KKEff<KTRAJ> {
    public:
      typedef KKEff<KTRAJ> KKEFF;
      typedef PKTraj<KTRAJ> PKTRAJ;
      typedef typename KKEFF::PDATA PDATA; // forward the typedefs
      typedef typename KKEFF::WDATA WDATA;
      typedef typename KKEFF::KKDATA KKDATA;
      // process this effect given the adjacent effect
      virtual bool process(KKDATA& kkdata,TDir tdir) override;
      virtual unsigned nDOF() const { return 0; } // no information content, therefore no NDOF
      virtual double chisq(PDATA const& pars) const { return 0.0; }
      virtual bool append(PKTRAJ& fit) override;
      // construct from a pdata
      KKPEff(PDATA const& pdata) : pdata_(pdata) {}
      // accessors
      PDATA const& pData() const { return pdata_; }
      virtual ~KKPEff(){}
    protected:
      KKPEff() {}
      // reset the cached weight; this must be done for each update cycle
      void resetCache() { wdata_ = WDATA(); }
      PDATA pdata_; // parameter space description of this effect
      WDATA wdata_; // cache of weight processing in opposite directions, used to build the fit trajectory
  };

  template<> bool KKPEff<KTRAJ>::process(KKDATA& kkdata,TDir tdir) {
    bool retval(false);
    if(this->isActive()){
      // backwards, set the cache BEFORE processing this effect, to avoid double-counting it
      if(tdir == TDir::backwards)wdata_ += kkdata.wData();
      kkdata.append(pdata_);
      // forwards, set the cache to AFTER processing this effect
      if(tdir == TDir::forwards)wdata_ += kkdata.wData();
      retval = kkdata.pData().matrixOK();
    }
    KKEffBase::setStatus(tdir,KKEffBase::processed);
    return retval;
  }

  template<> bool KKPEff<KTRAJ>::append(PKTRAJ& fit) override {
    // create a trajectory piece from the cached weight
    KTRAJ endpiece(KKEFF::refTraj());
    endpiece.params() = PDATA(wdata_);
    // adjust the range appropriately
    endpiece.range() = TRange(time(),fit.range().high());
    // append this to the fit
    fit.append(endpiece);
  }

}
#endif
