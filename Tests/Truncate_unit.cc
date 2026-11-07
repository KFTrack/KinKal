//
// Test Track::truncate: fit a ToyMC track, extrapolate it past its fitted range in each direction, truncate it inside
// the outermost extrapolated domain and, after extrapolating again, halfway back to the last active hit or material (with
// field correction this removes domains and their DomainWalls), and check the trajectory, the domains, the effects and a
// later extrapolation. Then make the last hits of a track inactive and check that a cut before them, but after every
// active hit and material, is accepted.
//
#include "KinKal/Trajectory/LoopHelix.hh"
#include "KinKal/Trajectory/CentralHelix.hh"
#include "KinKal/Trajectory/ParticleTrajectory.hh"
#include "KinKal/General/BFieldMap.hh"
#include "KinKal/Fit/Config.hh"
#include "KinKal/Fit/Track.hh"
#include "KinKal/Examples/StrawXing.hh"
#include "KinKal/Examples/SimpleWireHit.hh"
#include "KinKal/Fit/Measurement.hh"
#include "KinKal/Fit/DomainWall.hh"
#include "KinKal/Tests/ToyMC.hh"
#include <iostream>
#include <string>
#include <cmath>
#include <algorithm>
#include <limits>
#include <vector>

using namespace KinKal;

namespace {
  // extend in steps of step_ until the trajectory reaches tstop_ in the requested direction
  struct ReachTime {
    double tstop_;
    double step_ = 1.0;
    template <class KTRAJ> bool needsExtrapolation(ParticleTrajectory<KTRAJ> const& ptraj, TimeDir tdir) const {
      return tdir == TimeDir::forwards ? ptraj.range().end() < tstop_ : ptraj.range().begin() > tstop_;
    }
    double maxDt() const { return 1000.0; }
    double maxDtStep() const { return step_; }
    double dpTolerance() const { return 1.0e-4; }
  };

  int nfail(0);
  void check(bool ok, std::string const& what) {
    if(!ok){ ++nfail; std::cout << "FAILED: " << what << std::endl; }
  }

  template <class KTRAJ> bool contiguous(Track<KTRAJ> const& trk) {
    double prev = std::nan("");
    for(auto const& dom : trk.domains()){
      if(!std::isnan(prev) && dom->begin() != prev) return false;
      prev = dom->end();
    }
    return true;
  }

  template <class KTRAJ> bool isWall(std::unique_ptr<Effect<KTRAJ>> const& eff) {
    return dynamic_cast<DomainWall<KTRAJ> const*>(eff.get()) != nullptr;
  }

  // the outermost time of the active hits and materials (DomainWalls excluded): a cut must stay beyond it
  template <class KTRAJ> double lastActive(Track<KTRAJ> const& trk, bool fwd) {
    double edge = fwd ? trk.fitTraj().range().begin() : trk.fitTraj().range().end();
    for(auto const& eff : trk.effects())
      if(eff->active() && !isWall(eff)) edge = fwd ? std::max(edge,eff->time()) : std::min(edge,eff->time());
    return edge;
  }

  // truncate at cut and check the trajectory, the domains beyond the cut (removed with the DomainWalls into them), the effects and the state
  template <class KTRAJ> void checkCut(Track<KTRAJ>& trk, TimeDir tdir, double cut, std::string const& label) {
    bool fwd = tdir == TimeDir::forwards;
    auto const& ftraj = trk.fitTraj();
    auto before = ftraj.range();
    auto removed = [fwd,cut](Domain const& dom){ return fwd ? dom.begin() >= cut : dom.end() <= cut; };
    size_t npieces = ftraj.pieces().size(), neff = trk.effects().size(), ndom = trk.domains().size(), ndomcut(0), nwallcut(0);
    for(auto const& dom : trk.domains()) ndomcut += removed(*dom);
    for(auto const& eff : trk.effects()){
      auto dw = dynamic_cast<DomainWall<KTRAJ> const*>(eff.get());
      if(dw != nullptr) nwallcut += removed(fwd ? dw->nextDomain() : dw->prevDomain());
    }
    auto pos = ftraj.position3(cut);
    auto mom = ftraj.momentum3(cut);
    auto stat = trk.fitStatus().status_;
    check(trk.truncate(tdir,cut),label + " truncation");
    check((fwd ? ftraj.range().end() : ftraj.range().begin()) == cut,label + " trajectory ends at the cut");
    check((fwd ? ftraj.range().begin() : ftraj.range().end()) == (fwd ? before.begin() : before.end()),label + " other end unchanged");
    check(ftraj.pieces().size() <= npieces,label + " no pieces added");
    check(nwallcut == ndomcut,label + " one DomainWall into each domain beyond the cut");
    check(trk.effects().size() == neff - nwallcut && trk.domains().size() == ndom - ndomcut,label + " only the domains beyond the cut and their DomainWalls removed");
    for(auto const& eff : trk.effects()){
      auto dw = dynamic_cast<DomainWall<KTRAJ> const*>(eff.get());
      if(dw != nullptr) check(!removed(fwd ? dw->nextDomain() : dw->prevDomain()),label + " no DomainWall into a removed domain");
    }
    check((ftraj.position3(cut) - pos).R() == 0.0 && (ftraj.momentum3(cut) - mom).R() == 0.0,label + " state at the cut unchanged");
    check(trk.fitStatus().status_ == stat,label + " fit status unchanged");
    if(!trk.domains().empty()){
      auto const& outer = fwd ? *trk.domains().crbegin() : *trk.domains().cbegin();
      check(fwd ? outer->end() <= cut : outer->begin() >= cut,label + " no domain beyond the cut");
      check(contiguous(trk),label + " domains contiguous after truncation");
    }
  }

  template <class KTRAJ> void testDirection(Track<KTRAJ>& trk, TimeDir tdir, std::string const& label, bool bfcorr) {
    bool fwd = tdir == TimeDir::forwards;
    auto const& ftraj = trk.fitTraj();
    double tend = fwd ? ftraj.range().end() : ftraj.range().begin();
    ReachTime xtest{fwd ? tend + 20.0 : tend - 20.0};
    check(trk.extrapolate(tdir,xtest),label + " extrapolation");
    tend = fwd ? ftraj.range().end() : ftraj.range().begin();
    check(fwd ? tend >= xtest.tstop_ : tend <= xtest.tstop_,label + " extrapolation reaches the target time");
    if(bfcorr) check(!trk.domains().empty(),label + " extrapolation adds domains");
    double edge = lastActive(trk,fwd);
    // refused: at the last active hit or material, outside the trajectory, non-finite
    auto before = ftraj.range();
    size_t neff = trk.effects().size(), ndom = trk.domains().size();
    check(!trk.truncate(tdir,edge),label + " refuses a cut at the last active effect");
    check(!trk.truncate(tdir,fwd ? tend + 1.0 : tend - 1.0),label + " refuses a cut outside the trajectory");
    check(!trk.truncate(tdir,std::nan("")),label + " refuses a non-finite cut");
    check(ftraj.range().begin() == before.begin() && ftraj.range().end() == before.end() && trk.effects().size() == neff && trk.domains().size() == ndom,
        label + " refused cuts change nothing");
    // 1: inside the outermost domain (no domain removed), else halfway to the end
    double cut = 0.5*(edge + tend);
    if(!trk.domains().empty()){
      auto const& outer = fwd ? *trk.domains().crbegin() : *trk.domains().cbegin();
      if(fwd ? (outer->begin() >= edge && outer->end() > outer->begin()) : (outer->end() <= edge && outer->end() > outer->begin()))
        cut = std::clamp(outer->mid(),std::min(edge,tend),std::max(edge,tend));
    }
    checkCut(trk,tdir,cut,label + " outermost domain");
    // 2: extrapolate again, then cut halfway between the last active hit or material and the end: with field correction this
    // removes several domains
    check(trk.extrapolate(tdir,xtest),label + " extrapolation after truncation");
    check(fwd ? ftraj.range().end() >= xtest.tstop_ : ftraj.range().begin() <= xtest.tstop_,label + " re-extrapolation reaches the target time");
    check(contiguous(trk),label + " domains contiguous after re-extrapolation");
    tend = fwd ? ftraj.range().end() : ftraj.range().begin();
    cut = 0.5*(edge + tend);
    if(bfcorr){
      size_t nbeyond(0);
      for(auto const& dom : trk.domains()) nbeyond += fwd ? dom->begin() >= cut : dom->end() <= cut;
      check(nbeyond >= 2,label + " several domains beyond the halfway cut");
    }
    checkCut(trk,tdir,cut,label + " halfway");
    // a later extrapolation continues from the cut
    check(trk.extrapolate(tdir,xtest),label + " extrapolation after the second truncation");
    check(fwd ? ftraj.range().end() >= xtest.tstop_ : ftraj.range().begin() <= xtest.tstop_,label + " second re-extrapolation reaches the target time");
    check(contiguous(trk),label + " domains contiguous after the second re-extrapolation");
  }

  // the fit configuration of the Tests schedule files (MCAmbig.txt, MCAmbigFixedField.txt): convergence and divergence
  // limits loose enough for a seed smeared by 10 sigma, no end processing
  Config fitConfig(bool bfcorr) {
    Config config;
    config.maxniter_ = 10;
    config.dwt_ = 1.0e6;
    config.convdchisq_ = 1.0;
    config.divdchisq_ = 1.0e8;
    config.pdchisq_ = 1.0e8;
    config.divgap_ = 100.0;
    config.tol_ = 1.0e-4;
    config.minndof_ = 5;
    config.bfcorr_ = bfcorr;
    config.ends_ = false;
    for(double temp : {2.0,1.0,0.0}){
      MetaIterConfig miconfig(temp);
      miconfig.addUpdater(StrawXingConfig(0.3,5.0,10.0,true));
      config.schedule_.push_back(miconfig);
    }
    return config;
  }

  template <class KTRAJ> void testTraj(DVEC const& sigmas, bool bfcorr) {
    std::string label = KTRAJ::trajName() + (bfcorr ? " bfcorr" : " no bfcorr");
    // the uniform field of the Tests fit units: the ToyMC fits diverge in a gradient field (LoopHelixFit --Bgrad, KinKal v3.9.1);
    // with bfcorr the extrapolation still adds maxDtStep domains, so the domain truncation is exercised
    UniformBFieldMap bfield(VEC3(0.0,0.0,1.0));
    KKTest::ToyMC<KTRAJ> toy(bfield,105.0,-1,3000.0,123421,40,true,true,0.25,0.511);
    toy.setTolerance(1.0e-5);
    Config config = fitConfig(bfcorr);
    // a few events, so a seed whose fit fails does not decide the test
    unsigned ntested(0);
    for(unsigned ievt=0; ievt < 10 && ntested < 3; ++ievt){
      typename KKTest::ToyMC<KTRAJ>::HITCOL hits;
      typename KKTest::ToyMC<KTRAJ>::EXINGCOL xings;
      ParticleTrajectory<KTRAJ> tptraj;
      toy.simulateParticle(tptraj,hits,xings,true);
      double tmid = tptraj.range().mid();
      auto const& midhel = tptraj.nearestPiece(tmid);
      auto seedmom = midhel.momentum4(tmid);
      auto seedpos = midhel.position4(tmid);
      KTRAJ seed(seedpos,seedmom,midhel.charge(),bfield.fieldVect(seedpos.Vect()),TimeRange(hits.front()->time(),hits.back()->time()));
      toy.createSeed(seed,sigmas,10.0);
      Track<KTRAJ> trk(config,bfield,seed,hits,xings);
      if(!trk.fitStatus().usable()){
        std::cout << label << " event " << ievt << " skipped, fit not usable: " << trk.fitStatus() << std::endl;
        continue;
      }
      ++ntested;
      std::string elabel = label + " event " + std::to_string(ievt);
      testDirection(trk,TimeDir::forwards,elabel + " forwards",bfcorr);
      testDirection(trk,TimeDir::backwards,elabel + " backwards",bfcorr);
    }
    check(ntested == 3,label + " has 3 usable fits in 10 events");
    std::cout << label << ": " << ntested << " fits tested" << std::endl;
  }
  // The last NINACTIVE hits of a ToyMC track (no scintillator hit) are made inactive, and the straw crossings after the
  // last active hit are dropped, so only inactive effects follow it. After one long extrapolation step past the inactive
  // hits, a cut after every active hit and material but before an inactive hit must be accepted,
  // keep the inactive effects with usable reference pieces, and a cut at the last active effect must still be refused.
  template <class KTRAJ> void testInactive(DVEC const& sigmas, bool bfcorr) {
    using WIREHIT = SimpleWireHit<KTRAJ>;
    static const size_t NINACTIVE = 4;
    std::string label = KTRAJ::trajName() + (bfcorr ? " bfcorr" : " no bfcorr") + " inactive hits";
    UniformBFieldMap bfield(VEC3(0.0,0.0,1.0));
    KKTest::ToyMC<KTRAJ> toy(bfield,105.0,-1,3000.0,123421,40,true,false,0.25,0.511);
    toy.setTolerance(1.0e-5);
    Config config = fitConfig(bfcorr);
    unsigned ntested(0);
    for(unsigned ievt=0; ievt < 10 && ntested < 3; ++ievt){
      typename KKTest::ToyMC<KTRAJ>::HITCOL hits;
      typename KKTest::ToyMC<KTRAJ>::EXINGCOL xings;
      ParticleTrajectory<KTRAJ> tptraj;
      toy.simulateParticle(tptraj,hits,xings,true);
      if(hits.size() < 10) continue;
      std::sort(hits.begin(),hits.end(),[](auto const& a, auto const& b){ return a->time() < b->time(); });
      double tactive = hits[hits.size()-NINACTIVE-1]->time();
      bool wire(true);
      for(size_t ih = hits.size()-NINACTIVE; ih < hits.size(); ++ih){
        auto whit = std::dynamic_pointer_cast<WIREHIT>(hits[ih]);
        if(!whit){ wire = false; break; }
        // the same wire, drift and TOT, in the inactive state
        typename WIREHIT::PCA pca(tptraj,whit->wire(),CAHint(whit->time(),whit->time()),whit->precision());
        hits[ih] = std::make_shared<WIREHIT>(bfield,pca,WireHitState(WireHitState::inactive),whit->minDOCA(),whit->driftVelocity(),
            whit->timeVariance(),whit->tot(),whit->totVariance(),whit->cellRadius(),whit->id());
      }
      if(!wire) continue;
      xings.erase(std::remove_if(xings.begin(),xings.end(),[tactive](auto const& xing){ return xing->time() > tactive; }),xings.end());
      double tmid = tptraj.range().mid();
      auto const& midhel = tptraj.nearestPiece(tmid);
      auto seedpos = midhel.position4(tmid);
      KTRAJ seed(seedpos,midhel.momentum4(tmid),midhel.charge(),bfield.fieldVect(seedpos.Vect()),TimeRange(hits.front()->time(),tactive));
      toy.createSeed(seed,sigmas,10.0);
      Track<KTRAJ> trk(config,bfield,seed,hits,xings);
      std::string elabel = label + " event " + std::to_string(ievt);
      if(!trk.fitStatus().usable()){
        std::cout << elabel << " skipped, fit not usable: " << trk.fitStatus() << std::endl;
        continue;
      }
      ++ntested;
      auto const& ftraj = trk.fitTraj();
      ReachTime xtest{hits.back()->time() + 5.0, 50.0};
      check(trk.extrapolate(TimeDir::forwards,xtest),elabel + " extrapolation past the inactive hits");
      check(ftraj.range().end() > hits.back()->time(),elabel + " the trajectory reaches past the inactive hits");
      // after the extrapolation: the last active hit or material (DomainWalls do not block a cut) and the inactive hits
      double lastactive = lastActive(trk,true);
      std::vector<double> tinact;
      for(auto const& eff : trk.effects())
        if(!eff->active() && dynamic_cast<Measurement<KTRAJ> const*>(eff.get()) != nullptr) tinact.push_back(eff->time());
      check(tinact.size() == NINACTIVE,elabel + " the inactive hits are effects of the fit");
      std::sort(tinact.begin(),tinact.end());
      auto next = std::upper_bound(tinact.begin(),tinact.end(),lastactive);
      if(next == tinact.end()){
        check(false,elabel + " an inactive hit lies beyond every active hit and material");
        std::cout << "   last active hit or material " << lastactive << " last active hit " << tactive
          << " inactive hits " << (tinact.empty() ? 0.0 : tinact.front()) << " .. " << (tinact.empty() ? 0.0 : tinact.back()) << std::endl;
        continue;
      }
      double cut = 0.5*(lastactive + *next);
      size_t neff = trk.effects().size();
      check(!trk.truncate(TimeDir::forwards,lastactive),elabel + " refuses a cut at the last active effect");
      check(trk.truncate(TimeDir::forwards,cut),elabel + " truncation before the inactive hits");
      check(ftraj.range().end() == cut,elabel + " trajectory ends at the cut");
      check(trk.effects().size() == neff,elabel + " effects kept");
      for(auto const& eff : trk.effects()){
        if(eff->active() || eff->time() <= cut) continue;
        auto meas = dynamic_cast<Measurement<KTRAJ> const*>(eff.get());
        check(meas != nullptr && meas->hit()->refTrajPtr() && std::isfinite(meas->hit()->referenceTrajectory().momentum(eff->time())),
            elabel + " an inactive hit beyond the cut keeps a usable reference piece");
      }
      if(!trk.domains().empty()) check((*trk.domains().crbegin())->end() <= cut && contiguous(trk),elabel + " domains end at the cut");
      check(trk.fitStatus().usable(),elabel + " fit status still usable");
    }
    check(ntested == 3,label + " has 3 usable fits in 10 events");
    std::cout << label << ": " << ntested << " fits tested" << std::endl;
  }
}

int main(int, char**) {
  testTraj<LoopHelix>(DVEC(0.5, 0.5, 0.5, 0.5, 0.005, 0.5),true);
  testTraj<LoopHelix>(DVEC(0.5, 0.5, 0.5, 0.5, 0.005, 0.5),false);
  testTraj<CentralHelix>(DVEC(0.5, 0.003, 0.00001, 3.0 , 0.004, 0.1),true);
  testTraj<CentralHelix>(DVEC(0.5, 0.003, 0.00001, 3.0 , 0.004, 0.1),false);
  testInactive<LoopHelix>(DVEC(0.5, 0.5, 0.5, 0.5, 0.005, 0.5),true);
  testInactive<LoopHelix>(DVEC(0.5, 0.5, 0.5, 0.5, 0.005, 0.5),false);
  testInactive<CentralHelix>(DVEC(0.5, 0.003, 0.00001, 3.0 , 0.004, 0.1),true);
  testInactive<CentralHelix>(DVEC(0.5, 0.003, 0.00001, 3.0 , 0.004, 0.1),false);
  std::cout << (nfail == 0 ? "Truncate test passed" : "Truncate test FAILED: " + std::to_string(nfail) + " checks") << std::endl;
  return nfail == 0 ? 0 : 1;
}
