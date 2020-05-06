//
//  Toy MC for fit and hit testing
//
#include "TRandom3.h"
#include "KinKal/TLine.hh"
#include "KinKal/PKTraj.hh"
#include "KinKal/TPoca.hh"
#include "KinKal/StrawHit.hh"
#include "KinKal/StrawMat.hh"
#include "KinKal/ScintHit.hh"
#include "KinKal/BField.hh"
#include "KinKal/Vectors.hh"
#include "KinKal/D2T.hh"
#include "CLHEP/Units/PhysicalConstants.h"

namespace KKTest {
  using namespace KinKal;
  template <class KTRAJ> class ToyMC {
    public:
      typedef THit<KTRAJ> THIT;
      typedef std::shared_ptr<THIT> THITPTR;
      typedef DXing<KTRAJ> DXING;
      typedef std::shared_ptr<DXING> DXINGPTR;
      typedef StrawHit<KTRAJ> STRAWHIT;
      typedef std::shared_ptr<STRAWHIT> STRAWHITPTR;
      typedef ScintHit<KTRAJ> SCINTHIT;
      typedef std::shared_ptr<SCINTHIT> SCINTHITPTR;
      typedef StrawXing<KTRAJ> STRAWXING;
      typedef std::shared_ptr<STRAWXING> STRAWXINGPTR;
      typedef std::vector<THITPTR> THITCOL;
      typedef std::vector<DXINGPTR> DXINGCOL;
      typedef PKTraj<KTRAJ> PKTRAJ;
      typedef TPoca<PKTRAJ,TLine> TPOCA;
      typedef typename KTRAJ::DVEC DVEC;     
      // create from aseed
      ToyMC(BField const& bfield, double mom, int icharge, double zrange, int iseed, unsigned nhits, bool simmat, bool lighthit, double ambigdoca ,double simmass) : 
	bfield_(bfield), mom_(mom), icharge_(icharge),
	tr_(iseed), nhits_(nhits), simmat_(simmat), lighthit_(lighthit), ambigdoca_(ambigdoca), simmass_(simmass),
	sprop_(0.8*CLHEP::c_light), sdrift_(0.065), 
	zrange_(zrange), rmax_(800.0), rstraw_(2.5), rwire_(0.025), wthick_(0.015), sigt_(3.0), ineff_(0.1),
	momvar_(1.0), ttsig_(0.5), twsig_(10.0), shmax_(80.0), clen_(200.0), cprop_(0.8*CLHEP::c_light),
	osig_(10.0), ctmin_(0.5), ctmax_(0.8), tbuff_(0.1), tol_(0.01),
	smat_(matdb_,rstraw_, wthick_,rwire_),
	d2t_(sdrift_,sigt_*sigt_,rstraw_) {}

      // generate a straw at the given time.  direction and drift distance are random
      TLine generateStraw(PKTRAJ const& traj, double htime);
      // create a seed by randomizing the parameters
      void createSeed(KTRAJ& seed);
      void extendTraj(PKTRAJ& pktraj,double htime);
      void createTraj(PKTRAJ& pktraj);
      void createScintHit(PKTRAJ const& pktraj, THITCOL& thits);
      void simulateParticle(PKTRAJ& pktraj,THITCOL& thits, DXINGCOL& dxings);
      double createStrawMaterial(PKTRAJ& pktraj, STRAWXING const& sxing);
      // set functions, for special purposes
      void setSeedVar(double momvar) { momvar_ = momvar; }
      void setSmearSeed(bool smear) { smearseed_ = smear; }
      // accessors
      double shVar() const {return sigt_*sigt_;}
      double chVar() const {return ttsig_*ttsig_;}
      double rStraw() const { return rstraw_; }
      double zRange() const { return zrange_; }
      StrawMat const& strawMaterial() const { return smat_; }

    private:
      BField const& bfield_;
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
      StrawMat smat_; // straw material
      CVD2T d2t_;
  };

  template <class KTRAJ> TLine ToyMC<KTRAJ>::generateStraw(PKTRAJ const& traj, double htime) {
    // start with the true helix position at this time
    Vec4 hpos; hpos.SetE(htime);
    traj.position(hpos);
    Vec3 hdir = traj.direction(htime);
    // generate a random direction for the straw
    double eta = tr_.Uniform(-M_PI,M_PI);
    Vec3 sdir(cos(eta),sin(eta),0.0);
    // generate a random drift perp to this and the trajectory
    double rdrift = tr_.Uniform(-rstraw_,rstraw_);
    Vec3 drift = (sdir.Cross(hdir)).Unit();
    Vec3 dpos = hpos.Vect() + rdrift*drift;
    //  cout << "Generating hit at position " << dpos << endl;
    double dprop = tr_.Uniform(0.0,rmax_);
    Vec3 mpos = dpos + sdir*dprop;
    Vec3 vprop = sdir*sprop_;
    // measured time is after propagation and drift
    double tmeas = htime + dprop/sprop_ + fabs(rdrift)/sdrift_;
    // smear measurement time
    tmeas = tr_.Gaus(tmeas,sigt_);
    // range doesn't really matter
    TRange trange(tmeas-dprop/sprop_,tmeas+dprop/sprop_);
    // construct the trajectory for this hit
    return TLine(mpos,vprop,tmeas,trange);
  }

  template <class KTRAJ> void ToyMC<KTRAJ>::simulateParticle(PKTRAJ& pktraj,THITCOL& thits, DXINGCOL& dxings) {
    // create the seed first
    createTraj(pktraj);
    // divide time range
    double dt(0.0);
    if(nhits_ > 0) dt = (pktraj.range().range()-2*tbuff_)/(nhits_);
    Vec3 bsim;
    // create the hits (and associated materials)
    for(size_t ihit=0; ihit<nhits_; ihit++){
      double htime = tbuff_ + pktraj.range().low() + ihit*dt;
      // extend the trajectory in the BField to this time
      extendTraj(pktraj,htime);
      // create the hit at this time
      auto tline = generateStraw(pktraj,htime);
      TPOCA tp(pktraj,tline);
      LRAmbig ambig(LRAmbig::null);
      if(fabs(tp.doca())> ambigdoca_) ambig = tp.doca() < 0 ? LRAmbig::left : LRAmbig::right;
      // construct the hit from this trajectory
      auto sxing = std::make_shared<STRAWXING>(tp,smat_);
      if(tr_.Uniform(0.0,1.0) > ineff_){
	thits.push_back(std::make_shared<STRAWHIT>(bfield_, tline, d2t_,sxing,ambig));
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
    if(thits.size() == 0) extendTraj(pktraj,pktraj.range().high());
  }

  template <class KTRAJ> double ToyMC<KTRAJ>::createStrawMaterial(PKTRAJ& pktraj, STRAWXING const& sxing) {
    double desum = 0.0;
    double tstraw = sxing.crossingTime();
    auto const& endpiece = pktraj.nearestPiece(tstraw);
    double mom = endpiece.momentumMag(tstraw);
    Mom4 endmom = endpiece.momentum(tstraw);
    Vec4 endpos; endpos.SetE(tstraw);
    endpiece.position(endpos);
    std::array<double,3> dmom = {0.0,0.0,0.0}, momvar {0.0,0.0,0.0};
    sxing.momEffects(pktraj,TDir::forwards, dmom, momvar);
    for(int idir=0;idir<=LocalBasis::phidir; idir++) {
      auto mdir = static_cast<LocalBasis::LocDir>(idir);
      double momsig = sqrt(momvar[idir]);
      double dm;
      // generate a random effect given this variance and mean.  Note momEffect is scaled to momentum
      switch( mdir ) {
	case KinKal::LocalBasis::perpdir: case KinKal::LocalBasis::phidir :
	  dm = tr_.Gaus(dmom[idir],momsig);
	  break;
	case KinKal::LocalBasis::momdir :
	  dm = std::min(0.0,tr_.Gaus(dmom[idir],momsig));
	  desum += dm*mom;
	  break;
	default:
	  throw std::invalid_argument("Invalid direction");
      }
      //	cout << "mom change dir " << LocalBasis::directionName(mdir) << " mean " << dmom[idir]  << " +- " << momsig << " value " << dm  << endl;
      Vec3 dmvec;
      DVEC pder;
      endpiece.momDeriv(tstraw,mdir,pder,dmvec);
      dmvec *= dm*mom;
      endmom.SetCoordinates(endmom.Px()+dmvec.X(), endmom.Py()+dmvec.Y(), endmom.Pz()+dmvec.Z(),endmom.M());
    }
    // generate a new piece and append
    KTRAJ newend(endpos,endmom,endpiece.charge(),endpiece.bnom(),TRange(tstraw,pktraj.range().high()));
    //      newend.print(cout,1);
    pktraj.append(newend);
    return desum/mom;
  }

  template <class KTRAJ> void ToyMC<KTRAJ>::createScintHit(PKTRAJ const& pktraj, THITCOL& thits) {
    // create a ScintHit at the end, axis parallel to z
    // first, find the position at showermax_.
    Vec3 shmpos, hend, lmeas;
    double cstart = pktraj.range().high() + tbuff_;
    hend = pktraj.position(cstart);
    double ltime = cstart + shmax_/pktraj.speed(cstart);
    shmpos = pktraj.position(ltime); // true position at shower-max
    // smear the x-y position by the transverse variance.
    lmeas.SetX(tr_.Gaus(shmpos.X(),twsig_));
    lmeas.SetY(tr_.Gaus(shmpos.Y(),twsig_));
    // set the z position to the sensor plane (end of the crystal)
    lmeas.SetZ(hend.Z()+clen_);
    // set the measurement time to correspond to the light propagation from showermax_, smeared by the resolution
    double tmeas = tr_.Gaus(ltime+(lmeas.Z()-shmpos.Z())/cprop_,ttsig_);
    // create the ttraj for the light propagation
    Vec3 lvel(0.0,0.0,cprop_);
    TRange trange(cstart,cstart+clen_/cprop_);
    TLine lline(lmeas,lvel,tmeas,trange);
    // then create the hit and add it; the hit has no material
    thits.push_back(std::make_shared<SCINTHIT>(lline, ttsig_*ttsig_, twsig_*twsig_));
    // test
    //    cout << "cstart " << cstart << " pos " << hend << endl;
    //    cout << "shmax_ " << ltime << " pos " << shmpos  << endl;
    //    Vec3 lhpos;
    //    lline.position(tmeas,lhpos);
    //    cout << "tmeas " <<  tmeas  << " pos " << lmeas  << " llinepos " << lhpos << endl;
    //    RESIDUAL lres;
    //    thits.back()->resid(pktraj,lres);
    //    cout << "ScintHit " << lres << endl;
    //    TPOCA tpl(pktraj,lline);
    //    cout <<"Light TPOCA ";
    //    tpl.print(cout,2);
  }

  template <class KTRAJ> void ToyMC<KTRAJ>::createSeed(KTRAJ& seed){
    auto& seedpar = seed.params();
    // propagate the momentum and position variances to parameter variances
    DVEC pder;
    Vec3 unit;
    for(int idir=0;idir<LocalBasis::ndir;idir++){
      seed.momDeriv(seed.range().mid(),LocalBasis::LocDir(idir),pder,unit);
      // convert derivative vector to a Nx1 matrix
      ROOT::Math::SMatrix<double,KTRAJ::NParams(),1> dPdm;
      dPdm.Place_in_col(pder,0,0);
      ROOT::Math::SMatrix<double, 1,1, ROOT::Math::MatRepSym<double,1> > MVar;
      MVar(0,0) = momvar_/(mom_*mom_);
      auto cpars =ROOT::Math::Similarity(dPdm,MVar);
      for(size_t ipar=0;ipar<KTRAJ::NParams();ipar++)
	seedpar.covariance()[ipar][ipar] += cpars[ipar][ipar];
    }
    // smearing on T0 from momentum is too small
    size_t it0 = KTRAJ::NParams()-1; // assumed t0 is the last index FIXME!
    seedpar.covariance()[it0][it0]*= 100;

    // now, randomize the parameters within those errors.  Don't include correlations
    if(smearseed_){
      for(unsigned ipar=0;ipar < 6; ipar++){
	double perr = sqrt(seedpar.covariance()[ipar][ipar]);
	seedpar.parameters()[ipar] += tr_.Gaus(0.0,perr);
      }
    }
  }

  template <class KTRAJ> void ToyMC<KTRAJ>::extendTraj(PKTRAJ& pktraj,double htime) {
    ROOT::Math::SMatrix<double,3> bgrad;
    Vec3 pos,vel, dBdt;
    pos = pktraj.position(htime);
    vel = pktraj.velocity(htime);
    dBdt = bfield_.fieldDeriv(pos,vel);
//    std::cout << "end time " << pktraj.back().range().low() << " hit time " << htime << std::endl;
    if(dBdt.R() != 0.0){
      TRange prange(pktraj.back().range().low(),pktraj.back().range().low());
      pktraj.back().rangeInTolerance(prange,bfield_, tol_);
      prange.low() = prange.high();
      do {
	pktraj.back().rangeInTolerance(prange,bfield_, tol_);
//	std::cout << " Range " << prange << std::endl;
	Vec4 pos; pos.SetE(prange.low());
	Mom4 mom =  pktraj.momentum(prange.low());
	pktraj.position(pos);
	Vec3 bf = bfield_.fieldVect(pos.Vect());
	KTRAJ newend(pos,mom,pktraj.charge(),bf,prange);
	pktraj.append(newend);
	prange.low() = prange.high();
      } while(prange.high() < htime);
    }
//    std::cout << "Extended traj " << pktraj << std::endl;
  }

  template <class KTRAJ> void ToyMC<KTRAJ>::createTraj(PKTRAJ& pktraj) {
    // randomize the position and momentum
    double tphi = tr_.Uniform(-M_PI,M_PI);
    double tcost = tr_.Uniform(ctmin_,ctmax_);
    double tsint = sqrt(1.0-tcost*tcost);
    Mom4 tmomv(mom_*tsint*cos(tphi),mom_*tsint*sin(tphi),mom_*tcost,simmass_);
    double tmax = fabs(zrange_/(CLHEP::c_light*tcost));
    Vec4 torigin(tr_.Gaus(0.0,osig_), tr_.Gaus(0.0,osig_), tr_.Gaus(-0.5*zrange_,osig_),tr_.Uniform(-tmax,tmax));
    Vec3 bsim = bfield_.fieldVect(torigin.Vect());
    KTRAJ ktraj(torigin,tmomv,icharge_,bsim,TRange(torigin.T(),torigin.T()+tmax+tbuff_));
    pktraj = PKTRAJ(ktraj);
  }
}


