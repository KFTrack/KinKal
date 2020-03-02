#ifndef KinKal_KKPEff_hh
#define KinKal_KKPEff_hh
//
// Class to describe parameter transport/noise i
// This effect provides information content and is processed in params space 
//
#include "KinKal/KKEff.hh"
namespace KinKal {
  template<class KTRAJ> class KKPEff : public KKEff<KTRAJ> {
    public:
      typedef KKEff<KTRAJ> KKEFF;
      typedef typename KKEFF::PDATA PDATA; // forward the typedef
      typedef typename KKEFF::KKDATA KKDATA; // forward the typedef
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
      PDATA pdata_; // parameter data for  this effect
      std::array<WDATA,2> weights_; // cache of weight processing in opposite directions, used to build the traj piece
  };

  template<> bool KKPEff<KTRAJ>::process(KKDATA& kkdata,TDir tdir) {
    bool retval(false);
    if(this->isActive()){
      kkdata.addPData(pdata_);
      retval = kkdata.pData().matrixOK();
    }
    KKEffBase::setStatus(tdir,KKEffBase::processed);
    // set the cache for this direction
    auto idir = static_cast<std::underlying_type<TDir>::type>(tdir);
    if(kkdata.hasWData()){

    }
    return retval;
  }

  template<> bool KKPEff<KTRAJ>::append(PKTRAJ& fit) override {
  // 

  }

}
#endif
