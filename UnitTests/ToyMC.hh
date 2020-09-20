#ifndef KinKal_ToyMC_hh
#define KinKal_ToyMC_hh
//
//  Toy MC for fit and hit testing
//
#include "TRandom3.h"
#include "Trajectory/Line.hh"
#include "Trajectory/ParticleTrajectory.hh"
#include "Trajectory/PieceClosestApproach.hh"
#include "Detector/StrawHit.hh"
#include "Detector/StrawXing.hh"
#include "Detector/StrawMat.hh"
#include "Detector/ScintHit.hh"
#include "Detector/BFieldMap.hh"
#include "Detector/BFieldUtils.hh"
#include "General/Vectors.hh"
#include "Detector/WireCell.hh"
#include "CLHEP/Units/PhysicalConstants.h"

namespace KKTest {
  using namespace KinKal;
  template <class KTRAJ> class ToyMC {
    public:
      using PKTRAJ = ParticleTrajectory<KTRAJ>;
      using DHIT = DetectorHit<KTRAJ>;
      using DHITPTR = std::shared_ptr<DHIT>;
      using DHITCOL = std::vector<DHITPTR>;
      using DXING = DetectorXing<KTRAJ>;
      using DXINGPTR = std::shared_ptr<DXING>;
      using DXINGCOL = std::vector<DXINGPTR>;
      using STRAWHIT = StrawHit<KTRAJ>;
      using STRAWHITPTR = std::shared_ptr<STRAWHIT>;
      using SCINTHIT = ScintHit<KTRAJ>;
      using SCINTHITPTR = std::shared_ptr<SCINTHIT>;
      using STRAWXING = StrawXing<KTRAJ>;
      using STRAWXINGPTR = std::shared_ptr<STRAWXING>;
      using PTCA = PieceClosestApproach<KTRAJ,Line>;
      // create from aseed
      ToyMC(BFieldMap const& bfield, double mom, int icharge, double zrange, int iseed, unsigned nhits, bool simmat, bool lighthit, double ambigdoca ,double simmass) : 
	bfield_(bfield), mom_(mom), icharge_(icharge),
	tr_(iseed), nhits_(nhits), simmat_(simmat), lighthit_(lighthit), ambigdoca_(ambigdoca), simmass_(simmass),
	sprop_(0.8*CLHEP::c_light), sdrift_(0.065), 
	zrange_(zrange), rmax_(800.0), rstraw_(2.5), rwire_(0.025), wthick_(0.015), sigt_(3.0), ineff_(0.1),
	momvar_(1.0), ttsig_(0.5), twsig_(10.0), shmax_(80.0), clen_(200.0), cprop_(0.8*CLHEP::c_light),
	osig_(10.0), ctmin_(0.5), ctmax_(0.8), tbuff_(0.01), tol_(1e-4), tprec_(1e-8),
	smat_(matdb_,rstraw_, wthick_,rwire_),
	cell_(sdrift_,sigt_*sigt_,rstraw_) {}

      // generate a straw at the given time.  direction and drift distance are random
      Line generateStraw(PKTRAJ const& traj, double htime);
      // create a seed by randomizing the parameters
      void createSeed(KTRAJ& seed);
      void extendTraj(PKTRAJ& pktraj,double htime);
      void createTraj(PKTRAJ& pktraj);
      void createScintHit(PKTRAJ const& pktraj, DHITCOL& thits);
      void simulateParticle(PKTRAJ& pktraj,DHITCOL& thits, DXINGCOL& dxings);
      double createStrawMaterial(PKTRAJ& pktraj, STRAWXING const& sxing);
      // set functions, for special purposes
      void setSeedVar(double momvar) { momvar_ = momvar; }
      void setSmearSeed(bool smear) { smearseed_ = smear; }
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
      bool simmat_, lighthit_, smearseed_;
      double ambigdoca_, simmass_;
      double sprop_; // propagation speed along straw
      double sdrift_; // drift speed inside straw
      double zrange_, rmax_, rstraw_; // tracker dimension
      double rwire_, wthick_;
      double sigt_; // drift time resolution in ns
      double ineff_; // hit inefficiency
      double momvar_; // seed randomization factors
      // time hit parameters
      double ttsig_, twsig_, shmax_, clen_, cprop_;
      double osig_, ctmin_, ctmax_;
      double tbuff_;
      double tol_; // tolerance on spatial accuracy for 
      double tprec_; // time precision on TCA
      StrawMat smat_; // straw material
      SimpleCell cell_;
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
    // range doesn't really matter
    TimeRange trange(tmeas-dprop/sprop_,tmeas+dprop/sprop_);
    // construct the trajectory for this hit
    return Line(mpos,vprop,tmeas,trange);
  }

  template <class KTRAJ> void ToyMC<KTRAJ>::simulateParticle(PKTRAJ& pktraj,DHITCOL& thits, DXINGCOL& dxings) {
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
      LRAmbig ambig(LRAmbig::null);
      if(fabs(tp.doca())> ambigdoca_) ambig = tp.doca() < 0 ? LRAmbig::left : LRAmbig::right;
      // construct the hit from this trajectory
      auto sxing = std::make_shared<STRAWXING>(tp,smat_);
      if(tr_.Uniform(0.0,1.0) > ineff_){
	thits.push_back(std::make_shared<STRAWHIT>(bfield_, tline, cell_,sxing,ambig));
      } else {
	dxings.push_back(sxing);
      }
      // compute material effects and change trajectory accordingly
      if(simmat_){
	double defrac = createStrawMaterial(pktraj, *sxing.get());
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

  template <class KTRAJ> double ToyMC<KTRAJ>::createStrawMaterial(PKTRAJ& pktraj, STRAWXING const& sxing) {
    double desum = 0.0;
    double tstraw = sxing.crossingTime();
    auto const& endpiece = pktraj.nearestPiece(tstraw);
    double mom = endpiece.momentum(tstraw);
    auto endmom = endpiece.momentum4(tstraw);
    auto endpos = endpiece.position4(tstraw);
    std::array<double,3> dmom {0.0,0.0,0.0}, momvar {0.0,0.0,0.0};
    sxing.materialEffects(pktraj,TimeDir::forwards, dmom, momvar);
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
      //	cout << "mom change dir " << MomBasis::directionName(mdir) << " mean " << dmom[idir]  << " +- " << momsig << " value " << dm  << endl;
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

  template <class KTRAJ> void ToyMC<KTRAJ>::createScintHit(PKTRAJ const& pktraj, DHITCOL& thits) {
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
    TimeRange trange(cstart,cstart+clen_/cprop_);
    Line lline(lmeas,lvel,tmeas,trange);
    // then create the hit and add it; the hit has no material
    thits.push_back(std::make_shared<SCINTHIT>(lline, ttsig_*ttsig_, twsig_*twsig_));
    // test
    //    cout << "cstart " << cstart << " pos " << hend << endl;
    //    cout << "shmax_ " << ltime << " pos " << shmpos  << endl;
    //    VEC3 lhpos;
    //    lline.position3(tmeas,lhpos);
    //    cout << "tmeas " <<  tmeas  << " pos " << lmeas  << " llinepos " << lhpos << endl;
    //    Residual lres;
    //    thits.back()->resid(pktraj,lres);
    //    cout << "ScintHit " << lres << endl;
    //    PTCA tpl(pktraj,lline);
    //    cout <<"Light PTCA ";
    //    tpl.print(cout,2);
  }

  template <class KTRAJ> void ToyMC<KTRAJ>::createSeed(KTRAJ& seed){
    auto& seedpar = seed.params();
    // propagate the momentum and position variances to parameter variances
    for(int idir=0;idir<MomBasis::ndir;idir++){
      DVEC pder = seed.momDeriv(seed.range().mid(),MomBasis::Direction(idir));
      // convert derivative vector to a Nx1 matrix
      ROOT::Math::SMatrix<double,NParams(),1> dPdm;
      dPdm.Place_in_col(pder,0,0);
      ROOT::Math::SMatrix<double, 1,1, ROOT::Math::MatRepSym<double,1> > MVar;
      MVar(0,0) = momvar_/(mom_*mom_);
      auto cpars =ROOT::Math::Similarity(dPdm,MVar);
      for(size_t ipar=0;ipar<NParams();ipar++)
	seedpar.covariance()[ipar][ipar] += cpars[ipar][ipar];
    }
    // smearing on T0 from momentum is too small
    size_t it0 = NParams()-1; // assumed t0 is the last index FIXME! this should be a constexpr
    seedpar.covariance()[it0][it0]*= 100; // smear factor should be adjustable FIXME!

    // now, randomize the parameters within those errors.  Don't include correlations
    if(smearseed_){
      for(unsigned ipar=0;ipar < NParams(); ipar++){
	double perr = sqrt(seedpar.covariance()[ipar][ipar]);
	seedpar.parameters()[ipar] += tr_.Gaus(0.0,perr);
      }
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
