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
      typedef KKEff<KTRAJ> KKEFF;
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
      virtual float time() const override { return kkhit_.time(); }
      virtual unsigned nDOF() const override { return kkhit_.nDOF(); }
      virtual bool isActive() const override { return kkhit_.isActive(); }
      virtual void process(KKDATA& kkdata,TDir tdir) override;
      virtual float fitChi() const override { return kkhit_.fitChi(); }
      virtual float chisq(PDATA const& pdata) const override { return kkhit_.chisq(pdata); }
      virtual void update(PKTRAJ const& ref) override;
      virtual void update(PKTRAJ const& ref, MConfig const& mconfig) override;
      virtual void append(PKTRAJ& fit) override { return kkmat_.append(fit); }
      virtual void print(std::ostream& ost=std::cout,int detail=0) const override;
      // accessors
      KKHIT const& hit() const { return kkhit_; }
      KKMAT const& mat() const { return kkmat_; }
    private:
      KKHIT kkhit_; // associated hit
      KKMAT kkmat_; // associated material
  };

  template <class KTRAJ> KKMHit<KTRAJ>::KKMHit(THITPTR const& thit, PKTRAJ const& pktraj) : kkhit_(thit,pktraj),
    kkmat_(thit->detCrossing(), pktraj, thit->isActive()) { update(pktraj); }

  template <class KTRAJ> void KKMHit<KTRAJ>::process(KKDATA& kkdata,TDir tdir) {
    // process in a fixed order to make material caching work
    bool hitfirst = (tdir == TDir::forwards && kkhit_.time() < kkmat_.time()) ||
      (tdir == TDir::backwards && kkhit_.time() > kkmat_.time());
    if(hitfirst) {
      kkhit_.process(kkdata,tdir);
      kkmat_.process(kkdata,tdir);
    } else { 
      kkmat_.process(kkdata,tdir);
      kkhit_.process(kkdata,tdir);
    }
    KKEffBase::setStatus(tdir,KKEffBase::processed);
  }

  template <class KTRAJ> void KKMHit<KTRAJ>::update(PKTRAJ const& pktraj) {
    if(pktraj.range().infinite())throw std::invalid_argument("Invalid range");
    // update the hit first, then use that to update the material 
    KKEffBase::updateStatus();
    kkhit_.update(pktraj);
    kkmat_.setTime(kkhit_.time());
    kkmat_.update(pktraj);
  }
  
  template <class KTRAJ> void KKMHit<KTRAJ>::update(PKTRAJ const& pktraj, MConfig const& mconfig) {
    KKEffBase::updateStatus();
    kkhit_.update(pktraj,mconfig);
    kkmat_.setTime(kkhit_.time());
    kkmat_.update(pktraj,mconfig);
  }

  template <class KTRAJ> void KKMHit<KTRAJ>::print(std::ostream& ost, int detail) const {
    ost << "KKMHit " << static_cast<KKEff<KTRAJ> const&>(*this) << std::endl;
    hit().print(ost,detail);
    mat().print(ost,detail);
  }
  
  template <class KTRAJ> std::ostream& operator <<(std::ostream& ost, KKMHit<KTRAJ> const& kkmhit) {
    kkmhit.print(ost,0);
    return ost;
  }
}
#endif
