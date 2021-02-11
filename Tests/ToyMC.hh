#ifndef KinKal_ToyMC_hh
#define KinKal_ToyMC_hh
//
//  Toy MC for fit and hit testing
//
#include "TRandom3.h"
#include "KinKal/Trajectory/Line.hh"
#include "KinKal/Trajectory/ParticleTrajectory.hh"
#include "KinKal/Trajectory/PiecewiseClosestApproach.hh"
#include "KinKal/Tests/SimpleWireHit.hh"
#include "KinKal/Detector/StrawXing.hh"
#include "KinKal/Detector/StrawMat.hh"
#include "KinKal/Tests/ScintHit.hh"
#include "KinKal/Detector/BFieldMap.hh"
#include "KinKal/Detector/BFieldUtils.hh"
#include "KinKal/General/Vectors.hh"
#include "KinKal/General/PhysicalConstants.h"

namespace KKTest {
  using namespace KinKal;
  template <class KTRAJ> class ToyMC {
    public:
      using PKTRAJ = ParticleTrajectory<KTRAJ>;
      using HIT = Hit<KTRAJ>;
      using HITPTR = std::shared_ptr<HIT>;
      using HITCOL = std::vector<HITPTR>;
      using EXING = ElementXing<KTRAJ>;
      using EXINGPTR = std::shared_ptr<EXING>;
      using EXINGCOL = std::vector<EXINGPTR>;
      using WIREHIT = SimpleWireHit<KTRAJ>;
      using WIREHITPTR = std::shared_ptr<WIREHIT>;
      using SCINTHIT = ScintHit<KTRAJ>;
      using SCINTHITPTR = std::shared_ptr<SCINTHIT>;
      using STRAWXING = StrawXing<KTRAJ>;
      using STRAWXINGPTR = std::shared_ptr<STRAWXING>;
      using PTCA = PiecewiseClosestApproach<KTRAJ,Line>;
      // create from aseed
      ToyMC(BFieldMap const& bfield, double mom, int icharge, double zrange, int iseed, unsigned nhits, bool simmat, bool lighthit, bool nulltime, double ambigdoca ,double simmass) : 
	bfield_(bfield), mom_(mom), icharge_(icharge),
	tr_(iseed), nhits_(nhits), simmat_(simmat), lighthit_(lighthit), nulltime_(nulltime), ambigdoca_(ambigdoca), simmass_(simmass),
	sprop_(0.8*CLHEP::c_light), sdrift_(0.065), 
	zrange_(zrange), rmax_(800.0), rstraw_(2.5), rwire_(0.025), wthick_(0.015), wlen_(1000.0), sigt_(3.0), ineff_(0.05),
	ttsig_(0.5), twsig_(10.0), shmax_(80.0), clen_(200.0), cprop_(0.8*CLHEP::c_light),
	osig_(10.0), ctmin_(0.5), ctmax_(0.8), tbuff_(0.01), tol_(1e-4), tprec_(1e-8),
	smat_(matdb_,rstraw_, wthick_,rwire_) {}

      // generate a straw at the given time.  direction and drift distance are random
      Line generateStraw(PKTRAJ const& traj, double htime);
      // create a seed by randomizing the parameters
      void createSeed(KTRAJ& seed,DVEC const& sigmas, double seedsmear);
      void extendTraj(PKTRAJ& pktraj,double htime);
      void createTraj(PKTRAJ& pktraj);
      void createScintHit(PKTRAJ const& pktraj, HITCOL& thits);
      void simulateParticle(PKTRAJ& pktraj,HITCOL& thits, EXINGCOL& dxings);
      double createStrawMaterial(PKTRAJ& pktraj, const EXING* sxing);
      // set functions, for special purposes
      void setInefficiency(double ineff) { ineff_ = ineff; }
      // accessors
      double shVar() const {return sigt_*sigt_;}
      double chVar() const {return ttsig_*ttsig_;}
      double rStraw() const { return rstraw_; }
      double zRange() const { return zrange_; }
      StrawMat const& strawMaterial() const { return smat_; }

    private:
      BFieldMap const& bfield_;
      MatEnv::MatDBInfo matdb_;
      double mom_;
      int icharge_;
      TRandom3 tr_; // random number generator
      unsigned nhits_; // number of hits to simulate
      bool simmat_, lighthit_, nulltime_;
      double ambigdoca_, simmass_;
      double sprop_; // propagation speed along straw
      double sdrift_; // drift speed inside straw
      double zrange_, rmax_, rstraw_; // tracker dimension
      double rwire_, wthick_, wlen_; // wire radius, thickness, length
      double sigt_; // drift time resolution in ns
      double ineff_; // hit inefficiency
      // time hit parameters
      double ttsig_, twsig_, shmax_, clen_, cprop_;
      double osig_, ctmin_, ctmax_;
      double tbuff_;
      double tol_; // tolerance on spatial accuracy for 
      double tprec_; // time precision on TCA
      StrawMat smat_; // straw material
    
  };

  template <class KTRAJ> Line ToyMC<KTRAJ>::generateStraw(PKTRAJ const& traj, double htime) {
    // start with the true helix position at this time
    auto hpos = traj.position4(htime);
    auto hdir = traj.direction(htime);
    // generate a random direction for the straw
    double eta = tr_.Uniform(-M_PI,M_PI);
    VEC3 sdir(cos(eta),sin(eta),0.0);
    // generate a random drift perp to this and the trajectory
    double rdrift = tr_.Uniform(-rstraw_,rstraw_);
    VEC3 drift = (sdir.Cross(hdir)).Unit();
    VEC3 dpos = hpos.Vect() + rdrift*drift;
    //  cout << "Generating hit at position " << dpos << endl;
    double dprop = tr_.Uniform(0.0,rmax_);
    VEC3 mpos = dpos + sdir*dprop;
    VEC3 vprop = sdir*sprop_;
    // measured time is after propagation and drift
    double tmeas = htime + dprop/sprop_ + fabs(rdrift)/sdrift_;
    // smear measurement time
    tmeas = tr_.Gaus(tmeas,sigt_);
    // construct the trajectory for this hit
    return Line(mpos,tmeas,vprop,wlen_);
  }

  template <class KTRAJ> void ToyMC<KTRAJ>::simulateParticle(PKTRAJ& pktraj,HITCOL& thits, EXINGCOL& dxings) {
    // create the seed first
    createTraj(pktraj);
    // divide time range
    double dt(0.0);
    if(nhits_ > 0) dt = (pktraj.range().range()-2*tbuff_)/(nhits_);
    VEC3 bsim;
    // create the hits (and associated materials)
    for(size_t ihit=0; ihit<nhits_; ihit++){
      double htime = tbuff_ + pktraj.range().begin() + ihit*dt;
      // extend the trajectory in the BFieldMap to this time
      extendTraj(pktraj,htime);
      // create the hit at this time
      auto tline = generateStraw(pktraj,htime);
      CAHint tphint(htime,htime);
      PTCA tp(pktraj,tline,tphint,tprec_);
      WireHitState::LRAmbig ambig(WireHitState::null);
      if(fabs(tp.doca())> ambigdoca_) ambig = tp.doca() < 0 ? WireHitState::left : WireHitState::right;
      WireHitState::Dimension dim(WireHitState::time);
      if(ambig == WireHitState::null){
	if(nulltime_)
	  dim = WireHitState::both;
	else
	  dim = WireHitState::distance;
      }
      // null variance based on doca cuttoff
      double rmax = std::min(ambigdoca_,rstraw_);
      double nullvar = (rmax*rmax)/3.0; // range is +- doca
      // null time shift
      double nulldt = 0.5*ambigdoca_/sdrift_; // the shift should be the average drift time over this distance
      WireHitState whstate(ambig, dim, nullvar, nulldt);
      // construct the hit from this trajectory
      if(tr_.Uniform(0.0,1.0) > ineff_){
	thits.push_back(std::make_shared<WIREHIT>(bfield_, tline, whstate, sdrift_, sigt_*sigt_, rstraw_));
      }
      dxings.push_back(std::make_shared<STRAWXING>(tp,smat_));
      // compute material effects and change trajectory accordingly
      if(simmat_){
	double defrac = createStrawMaterial(pktraj, dxings.back().get());
	// terminate if there is catastrophic energy loss
	if(fabs(defrac) > 0.1)break;
      }
    }
    if(lighthit_ && tr_.Uniform(0.0,1.0) > ineff_){
      createScintHit(pktraj,thits);
    }
    // extend if there are no hits
    if(thits.size() == 0) extendTraj(pktraj,pktraj.range().end());
  }

  template <class KTRAJ> double ToyMC<KTRAJ>::createStrawMaterial(PKTRAJ& pktraj, const EXING* sxing) {
    double desum = 0.0;
    double tstraw = sxing->crossingTime();
    auto const& endpiece = pktraj.nearestPiece(tstraw);
    double mom = endpiece.momentum(tstraw);
    auto endmom = endpiece.momentum4(tstraw);
    auto endpos = endpiece.position4(tstraw);
    std::array<double,3> dmom {0.0,0.0,0.0}, momvar {0.0,0.0,0.0};
    sxing->materialEffects(pktraj,TimeDir::forwards, dmom, momvar);
    for(int idir=0;idir<=MomBasis::phidir_; idir++) {
      auto mdir = static_cast<MomBasis::Direction>(idir);
      double momsig = sqrt(momvar[idir]);
      double dm;
      // generate a random effect given this variance and mean.  Note momEffect is scaled to momentum
      switch( mdir ) {
	case KinKal::MomBasis::perpdir_: case KinKal::MomBasis::phidir_ :
	  dm = tr_.Gaus(dmom[idir],momsig);
	  break;
	case KinKal::MomBasis::momdir_ :
	  dm = std::min(0.0,tr_.Gaus(dmom[idir],momsig));
	  desum += dm*mom;
	  break;
	default:
	  throw std::invalid_argument("Invalid direction");
      }
      auto dmvec = endpiece.direction(tstraw,mdir);
      dmvec *= dm*mom;
      endmom.SetCoordinates(endmom.Px()+dmvec.X(), endmom.Py()+dmvec.Y(), endmom.Pz()+dmvec.Z(),endmom.M());
    }
    // generate a new piece and append
    KTRAJ newend(endpos,endmom,endpiece.charge(),endpiece.bnom(),TimeRange(tstraw,pktraj.range().end()));
    //      newend.print(cout,1);
    pktraj.append(newend);
    return desum/mom;
  }

  template <class KTRAJ> void ToyMC<KTRAJ>::createScintHit(PKTRAJ const& pktraj, HITCOL& thits) {
    // create a ScintHit at the end, axis parallel to z
    // first, find the position at showermax_.
    VEC3 shmpos, hend, lmeas;
    double cstart = pktraj.range().end() + tbuff_;
    hend = pktraj.position3(cstart);
    double ltime = cstart + shmax_/pktraj.speed(cstart);
    shmpos = pktraj.position3(ltime); // true position at shower-max
    // smear the x-y position by the transverse variance.
    lmeas.SetX(tr_.Gaus(shmpos.X(),twsig_));
    lmeas.SetY(tr_.Gaus(shmpos.Y(),twsig_));
    // set the z position to the sensor plane (end of the crystal)
    lmeas.SetZ(hend.Z()+clen_);
    // set the measurement time to correspond to the light propagation from showermax_, smeared by the resolution
    double tmeas = tr_.Gaus(ltime+(lmeas.Z()-shmpos.Z())/cprop_,ttsig_);
    // create the ttraj for the light propagation
    VEC3 lvel(0.0,0.0,cprop_);
//    TimeRange trange(cstart,cstart+clen_/cprop_);
    Line lline(lmeas,tmeas,lvel,clen_);
    // then create the hit and add it; the hit has no material
    thits.push_back(std::make_shared<SCINTHIT>(lline, ttsig_*ttsig_, twsig_*twsig_));
  }

  template <class KTRAJ> void ToyMC<KTRAJ>::createSeed(KTRAJ& seed,DVEC const& sigmas,double seedsmear){
    auto& seedpar = seed.params();
    // create covariance
    for(size_t ipar=0; ipar < NParams(); ipar++){
      double perr = sigmas[ipar]*seedsmear;
      seedpar.covariance()[ipar][ipar] = perr*perr;
      seedpar.parameters()[ipar] += tr_.Gaus(0.0,perr);
    }
  }

  template <class KTRAJ> void ToyMC<KTRAJ>::extendTraj(PKTRAJ& pktraj,double htime) {
    ROOT::Math::SMatrix<double,3> bgrad;
    VEC3 pos,vel, dBdt;
    pos = pktraj.position3(htime);
    vel = pktraj.velocity(htime);
    dBdt = bfield_.fieldDeriv(pos,vel);
//    std::cout << "end time " << pktraj.back().range().begin() << " hit time " << htime << std::endl;
    if(dBdt.R() != 0.0){
      TimeRange prange(pktraj.back().range().begin(),pktraj.back().range().begin());
      prange.end() = BFieldUtils::rangeInTolerance(prange.begin(), bfield_, pktraj.back(), tol_);
      if(prange.end() > htime) {
	return;
      } else {
	prange.begin() = prange.end();
	do {
	  prange.end() = BFieldUtils::rangeInTolerance(prange.begin(), bfield_, pktraj.back(), tol_);
	  VEC4 pos = pktraj.position4(prange.begin());
	  MOM4 mom =  pktraj.momentum4(prange.begin());
	  VEC3 bf = bfield_.fieldVect(pktraj.position3(prange.mid()));
	  KTRAJ newend(pos,mom,pktraj.charge(),bf,prange);
	  pktraj.append(newend);
	  prange.begin() = prange.end();
	} while(prange.end() < htime);
      }
    }
  }

  template <class KTRAJ> void ToyMC<KTRAJ>::createTraj(PKTRAJ& pktraj) {
    // randomize the position and momentum
    double tphi = tr_.Uniform(-M_PI,M_PI);
    double tcost = tr_.Uniform(ctmin_,ctmax_);
    double tsint = sqrt(1.0-tcost*tcost);
    MOM4 tmomv(mom_*tsint*cos(tphi),mom_*tsint*sin(tphi),mom_*tcost,simmass_);
    double tmax = fabs(zrange_/(CLHEP::c_light*tcost));
    VEC4 torigin(tr_.Gaus(0.0,osig_), tr_.Gaus(0.0,osig_), tr_.Gaus(-0.5*zrange_,osig_),tr_.Uniform(-tmax,tmax));
    VEC3 bsim = bfield_.fieldVect(torigin.Vect());
    KTRAJ ktraj(torigin,tmomv,icharge_,bsim,TimeRange(torigin.T(),torigin.T()+tmax+tbuff_));
    pktraj = PKTRAJ(ktraj);
  }
}
#endif
