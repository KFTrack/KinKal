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
      virtual unsigned nDOF() const override { return 0; } // no information content, therefore no NDOF
      virtual double chisq(PDATA const& pars) const override { return 0.0; }
      virtual bool append(PKTRAJ& fit) override;
      // construct from a pdata
      KKPEff(PDATA const& pdata) : pdata_(pdata) {}
      // accessors
      PDATA const& pData() const { return pdata_; }
      virtual ~KKPEff(){}
    protected:
      KKPEff() {}
      // reset the cache. this must be done for each update cycle
      void resetCache() { wdata_ = WDATA(); pdata_ = PDATA();}
      PDATA pdata_; // parameter space description of this effect
      WDATA wdata_; // cache of weight processing in opposite directions, used to build the fit trajectory
  };

  template<class KTRAJ> bool KKPEff<KTRAJ>::process(KKDATA& kkdata,TDir tdir) {
    bool retval(false);
    if(this->isActive()){
      // forwards, set the cache AFTER processing this effect
      if(tdir == TDir::forwards) {
	kkdata.append(pdata_);
	wdata_ += kkdata.wData();
      } else {
      // backwards, set the cache BEFORE processing this effect, to avoid double-counting it
	wdata_ += kkdata.wData();
	// reverse the sign; this is clumsy FIXME
	PData revpar(pdata_);
	revpar.parameters() *= -1.0;
	kkdata.append(revpar);
      }
      retval = kkdata.pData().matrixOK();
    }
    KKEffBase::setStatus(tdir,KKEffBase::processed);
    return retval;
  }

  template<class KTRAJ> bool KKPEff<KTRAJ>::append(PKTRAJ& fit) {
    // create a trajectory piece from the cached weight
    double time = this->time();
    KTRAJ newpiece(KKEFF::refTraj());
    newpiece.params() = PDATA(wdata_,true);
    // make sure there's enough range to append as a physical piece
    newpiece.range() = TRange(time,std::max(fit.range().high(),time+1.0));
    fit.append(newpiece);
    return true;
  }
}
#endif
