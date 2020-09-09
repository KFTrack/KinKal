#ifndef KinKal_MaterialHit_hh
#define KinKal_MaterialHit_hh
//
//  Class defining a hit with associated material
//
#include "KinKal/Hit.hh"
#include "KinKal/Material.hh"
#include "KinKal/Effect.hh"
#include <stdexcept>
#include <ostream>
#include <memory>

namespace KinKal {
  template <class KTRAJ> class MaterialHit : public Effect<KTRAJ> {
    public:
      using KKEFF = Effect<KTRAJ>;
      using KKHIT = Hit<KTRAJ>;
      using KKMAT = Material<KTRAJ>;
      using PKTRAJ = ParticleTrajectory<KTRAJ>;
      using THIT = DetectorHit<KTRAJ>;
      using THITPTR = std::shared_ptr<THIT>;
      MaterialHit(KKHIT& kkhit, KKMAT& kkmat) : kkhit_(kkhit), kkmat_(kkmat) {}
      MaterialHit(THITPTR const& thit, PKTRAJ const& reftraj);
      // override the interface
      virtual double time() const override { return kkhit_.time(); }
      virtual unsigned nDOF() const override { return kkhit_.nDOF(); }
      virtual bool isActive() const override { return kkhit_.isActive(); }
      virtual void process(FitData& kkdata,TimeDir tdir) override;
      virtual double fitChi() const override { return kkhit_.fitChi(); }
      virtual double chisq(Parameters const& pdata) const override { return kkhit_.chisq(pdata); }
      virtual void update(PKTRAJ const& ref) override;
      virtual void update(PKTRAJ const& ref, MetaIterConfig const& miconfig) override;
      virtual void append(PKTRAJ& fit) override { return kkmat_.append(fit); }
      virtual void print(std::ostream& ost=std::cout,int detail=0) const override;
      // accessors
      KKHIT const& hit() const { return kkhit_; }
      KKMAT const& mat() const { return kkmat_; }
    private:
      KKHIT kkhit_; // associated hit
      KKMAT kkmat_; // associated material
  };

  template <class KTRAJ> MaterialHit<KTRAJ>::MaterialHit(THITPTR const& thit, PKTRAJ const& pktraj) : kkhit_(thit,pktraj),
    kkmat_(thit->detXingPtr(), pktraj, thit->isActive()) { update(pktraj); }

  template <class KTRAJ> void MaterialHit<KTRAJ>::process(FitData& kkdata,TimeDir tdir) {
    // process in a fixed order to make material caching work
    bool hitfirst = (tdir == TimeDir::forwards && kkhit_.time() < kkmat_.time()) ||
      (tdir == TimeDir::backwards && kkhit_.time() > kkmat_.time());
    if(hitfirst) {
      kkhit_.process(kkdata,tdir);
      kkmat_.process(kkdata,tdir);
    } else { 
      kkmat_.process(kkdata,tdir);
      kkhit_.process(kkdata,tdir);
    }
    KKEFF::setStatus(tdir,KKEFF::processed);
  }

  template <class KTRAJ> void MaterialHit<KTRAJ>::update(PKTRAJ const& pktraj) {
    if(pktraj.range().infinite())throw std::invalid_argument("Invalid range");
    // update the hit first, then use that to update the material 
    KKEFF::updateStatus();
    kkhit_.update(pktraj);
    kkmat_.setTime(kkhit_.time());
    kkmat_.update(pktraj);
  }
  
  template <class KTRAJ> void MaterialHit<KTRAJ>::update(PKTRAJ const& pktraj, MetaIterConfig const& miconfig) {
    KKEFF::updateStatus();
    kkhit_.update(pktraj,miconfig);
    kkmat_.setTime(kkhit_.time());
    kkmat_.update(pktraj,miconfig);
  }

  template <class KTRAJ> void MaterialHit<KTRAJ>::print(std::ostream& ost, int detail) const {
    ost << "MaterialHit " << static_cast<Effect<KTRAJ> const&>(*this) << std::endl;
    hit().print(ost,detail);
    mat().print(ost,detail);
  }
  
  template <class KTRAJ> std::ostream& operator <<(std::ostream& ost, MaterialHit<KTRAJ> const& kkmhit) {
    kkmhit.print(ost,0);
    return ost;
  }
}
#endif
