#ifndef KinKal_KKMHit_hh
#define KinKal_KKMHit_hh
//
//  Class defining a hit with associated material
//
#include "KinKal/KKHit.hh"
#include "KinKal/KKMat.hh"
#include "KinKal/KKEff.hh"
#include <stdexcept>
namespace KinKal {
  template <class KTRAJ> class KKMHit : public KKEff<KTRAJ> {
    public:
      typedef KKHit<KTRAJ> KKHIT;
      typedef KKMat<KTRAJ> KKMAT;
      typedef PKTraj<KTRAJ> PKTRAJ;
      typedef typename KTRAJ::PDATA PDATA;
      typedef KKData<PDATA::PDim()> KKDATA;
      KKMHit(KKHIT const& kkhit, KKMAT const& kkmat) : kkhit_(kkhit), kkmat_(kkmat) {}
      KKMHit(THit const& thit, PKTRAJ const& reftraj);
      // override the interface
      virtual double time() const override { return kkhit_.time(); }
      virtual unsigned nDOF() const override { return kkhit_.nDOF(); }
      virtual bool isActive() const override { return kkhit_.isActive(); }
      virtual bool process(KKDATA& kkdata,TDir tdir) override;
      virtual double chisq(PDATA const& pars) const override { return kkhit_.chisq(pars); }
      virtual bool update(PKTRAJ const& ref) override;
      virtual bool append(PKTRAJ& fit) override { return kkmat_.append(fit); }
      // accessors
      KKHIT const& hit() const { return kkhit_; }
      KKMAT const& mat() const { return kkmat_; }
    private:
      KKHIT kkhit_; // associated hit
      KKMAT kkmat_; // associated material
  };

  template <class KTRAJ> KKMHit<KTRAJ>::KKMHit(THit const& thit, PKTRAJ const& reftraj) : kkhit_(thit,reftraj),
    kkmat_(*thit.material(),kkhit_.poca(),thit.isActive()) {}

  template <class KTRAJ> bool KKMHit<KTRAJ>::process(KKDATA& kkdata,TDir tdir) {
    // process in a fixed order to make material caching work
    bool retval(true);
    bool hitfirst = (tdir == TDir::forwards && kkhit_.time() < kkmat_.time()) ||
      (tdir == TDir::backwards && kkhit_.time() > kkmat_.time());
    if(hitfirst) {
	retval &= kkhit_.process(kkdata,tdir);
	retval &= kkmat_.process(kkdata,tdir);
    } else { 
      retval &= kkmat_.process(kkdata,tdir);
      retval &= kkhit_.process(kkdata,tdir);
    }
    return retval;
  }

  template <class KTRAJ> bool KKMHit<KTRAJ>::update(PKTRAJ const& ref) {
    // update the hit first, then use the POCA from that to update the material
    bool retval(true);
    retval &= kkhit_.update(ref);
    retval &= kkmat_.update(kkhit_.poca());
    return retval;
  }
}
#endif
