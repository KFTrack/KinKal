#ifndef KinKal_KKMHit_hh
#define KinKal_KKMHit_hh
//
//  Class defining a hit with associated material
//
#include "KinKal/KKHit.hh"
#include "KinKal/KKMat.hh"
#include "KinKal/KKEff.hh"
#include <stdexcept>
#include <ostream>
#include <memory>

namespace KinKal {
  template <class KTRAJ> class KKMHit : public KKEff<KTRAJ> {
    public:
      typedef KKHit<KTRAJ> KKHIT;
      typedef KKMat<KTRAJ> KKMAT;
      typedef PKTraj<KTRAJ> PKTRAJ;
      typedef THit<KTRAJ> THIT;
      typedef std::shared_ptr<THIT> THITPTR;
      typedef typename KTRAJ::PDATA PDATA;
      typedef KKData<PDATA::PDim()> KKDATA;
      KKMHit(KKHIT& kkhit, KKMAT& kkmat) : kkhit_(kkhit), kkmat_(kkmat) {}
      KKMHit(THITPTR const& thit, PKTRAJ const& reftraj);
      // override the interface
      virtual double time() const override { return kkhit_.time(); }
      virtual unsigned nDOF() const override { return kkhit_.nDOF(); }
      virtual bool isActive() const override { return kkhit_.isActive(); }
      virtual bool process(KKDATA& kkdata,TDir tdir) override;
      virtual double chisq(PDATA const& pars) const override { return kkhit_.chisq(pars); }
      virtual bool update(PKTRAJ const& ref) override;
      virtual bool append(PKTRAJ& fit) override { return kkmat_.append(fit); }
      virtual void print(std::ostream& ost=std::cout,int detail=0) const override;
      // accessors
      KKHIT const& hit() const { return kkhit_; }
      KKMAT const& mat() const { return kkmat_; }
    private:
      KKHIT kkhit_; // associated hit
      KKMAT kkmat_; // associated material
  };

  template <class KTRAJ> KKMHit<KTRAJ>::KKMHit(THITPTR const& thit, PKTRAJ const& reftraj) : kkhit_(thit,reftraj),
    kkmat_(thit->detCrossing(), reftraj, kkhit_.poca(), thit->isActive()) {}

  template <class KTRAJ> bool KKMHit<KTRAJ>::process(KKDATA& kkdata,TDir tdir) {
    // process in a fixed order to make material caching work
    bool retval(true);
    bool hitfirst = (tdir == TDir::forwards && kkhit_.time() < kkmat_.time()) ||
      (tdir == TDir::backwards && kkhit_.time() > kkmat_.time());
    if(hitfirst) {
      retval &= kkhit_.process(kkdata,tdir);
      if(kkmat_.detXing().use_count() > 0) retval &= kkmat_.process(kkdata,tdir);
    } else { 
      if(kkmat_.detXing().use_count() > 0) retval &= kkmat_.process(kkdata,tdir);
      retval &= kkhit_.process(kkdata,tdir);
    }
    KKEffBase::setStatus(tdir,KKEffBase::processed);
    return retval;
  }

  template <class KTRAJ> bool KKMHit<KTRAJ>::update(PKTRAJ const& ref) {
    if(ref.range().infinite())throw std::invalid_argument("Invalid range");
    // update the hit first, then use the POCA from that to update the material
    bool retval(true);
    KKEffBase::updateStatus();
    retval &= kkhit_.update(ref);
    if(kkmat_.detXing().use_count() > 0) retval &= kkmat_.update(ref,kkhit_.poca());
    return retval;
  }

  template <class KTRAJ> void KKMHit<KTRAJ>::print(std::ostream& ost, int detail) const {
    ost << "KKMHit " << static_cast<KKEff<KTRAJ> const&>(*this) << std::endl;
    hit().print(ost,detail);
    if(kkmat_.detXing().use_count() > 0) mat().print(ost,detail);
  }
  
  template <class KTRAJ> std::ostream& operator <<(std::ostream& ost, KKMHit<KTRAJ> const& kkmhit) {
    kkmhit.print(ost,0);
    return ost;
  }
}
#endif
