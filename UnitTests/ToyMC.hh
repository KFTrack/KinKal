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
      ToyMC(BField const& bfield, float mom, int icharge, float zrange, int iseed, unsigned nhits, bool simmat, bool lighthit, float ambigdoca ,float simmass) : 
	bfield_(bfield), mom_(mom), icharge_(icharge),
	tr_(iseed), nhits_(nhits), simmat_(simmat), lighthit_(lighthit), ambigdoca_(ambigdoca), simmass_(simmass),
	sprop_(0.8*CLHEP::c_light), sdrift_(0.065), 
	zrange_(zrange), rmax_(800.0), rstraw_(2.5), rwire_(0.025), wthick_(0.015), sigt_(3.0), ineff_(0.1),
	momvar_(1.0), ttsig_(0.5), twsig_(10.0), shmax_(80.0), clen_(200.0), cprop_(0.8*CLHEP::c_light),
	osig_(3.0), ctmin_(0.5), ctmax_(0.8), tbuff_(0.1),
	smat_(matdb_,rstraw_, wthick_,rwire_),
	d2t_(sdrift_,sigt_*sigt_,rstraw_) {}

      // generate a straw at the given time.  direction and drift distance are random
      TLine generateStraw(PKTRAJ const& traj, float htime);
      // create a seed by randomizing the parameters
      void createSeed(KTRAJ& seed);
      void extendTraj(PKTRAJ& pktraj,double htime);
      void createTraj(PKTRAJ& pktraj);
      void createScintHit(PKTRAJ const& pktraj, THITCOL& thits);
      void simulateParticle(PKTRAJ& pktraj,THITCOL& thits, DXINGCOL& dxings);
      float createStrawMaterial(PKTRAJ& pktraj, STRAWXING const& sxing);
      // set functions, for special purposes
      void setSeedVar(float momvar) { momvar_ = momvar; }
      void setSmearSeed(bool smear) { smearseed_ = smear; }
      // accessors
      float shVar() const {return sigt_*sigt_;}
      float chVar() const {return ttsig_*ttsig_;}
      float rStraw() const { return rstraw_; }
      float zRange() const { return zrange_; }
      StrawMat const& strawMaterial() const { return smat_; }

    private:
      BField const& bfield_;
      MatEnv::MatDBInfo matdb_;
      float mom_;
      int icharge_;
      TRandom3 tr_; // random number generator
      unsigned nhits_; // number of hits to simulate
      bool simmat_, lighthit_, smearseed_;
      float ambigdoca_, simmass_;
      float sprop_; // propagation speed along straw
      float sdrift_; // drift speed inside straw
      float zrange_, rmax_, rstraw_; // tracker dimension
      float rwire_, wthick_;
      double sigt_; // drift time resolution in ns
      float ineff_; // hit inefficiency
      float momvar_; // seed randomization factors
      // time hit parameters
      float ttsig_, twsig_, shmax_, clen_, cprop_;
      float osig_, ctmin_, ctmax_;
      float tbuff_;
      StrawMat smat_; // straw material
      CVD2T d2t_;
  };

  template <class KTRAJ> TLine ToyMC<KTRAJ>::generateStraw(PKTRAJ const& traj, float htime) {
    // start with the true helix position at this time
    Vec4 hpos; hpos.SetE(htime);
    traj.position(hpos);
    Vec3 hdir; traj.direction(htime,hdir);
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
    //  cout << "Creating " << nhits_ << " hits " << endl;
    // divide time range
    double dt = (pktraj.range().range()-2*tbuff_)/(nhits_-1);
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
	float defrac = createStrawMaterial(pktraj, *sxing.get());
	// terminate if there is catastrophic energy loss
	if(fabs(defrac) > 0.1)break;
      }
    }
    if(lighthit_ && tr_.Uniform(0.0,1.0) > ineff_){
      createScintHit(pktraj,thits);
    }
  }

  template <class KTRAJ> float ToyMC<KTRAJ>::createStrawMaterial(PKTRAJ& pktraj, STRAWXING const& sxing) {
    double desum = 0.0;
    float tstraw = sxing.crossingTime();
    auto const& endpiece = pktraj.nearestPiece(tstraw);
    double mom = endpiece.momentum(tstraw);
    Mom4 endmom;
    endpiece.momentum(tstraw,endmom);
    Vec4 endpos; endpos.SetE(tstraw);
    endpiece.position(endpos);
    std::array<float,3> dmom = {0.0,0.0,0.0}, momvar {0.0,0.0,0.0};
    sxing.momEffects(pktraj,TDir::forwards, dmom, momvar);
    for(int idir=0;idir<=KInter::theta2; idir++) {
      auto mdir = static_cast<KInter::MDir>(idir);
      double momsig = sqrt(momvar[idir]);
      double dm;
      // generate a random effect given this variance and mean.  Note momEffect is scaled to momentum
      switch( mdir ) {
	case KinKal::KInter::theta1: case KinKal::KInter::theta2 :
	  dm = tr_.Gaus(dmom[idir],momsig);
	  break;
	case KinKal::KInter::momdir :
	  dm = std::min(0.0,tr_.Gaus(dmom[idir],momsig));
	  desum += dm*mom;
	  break;
	default:
	  throw std::invalid_argument("Invalid direction");
      }
      //	cout << "mom change dir " << KInter::directionName(mdir) << " mean " << dmom[idir]  << " +- " << momsig << " value " << dm  << endl;
      Vec3 dmvec;
      DVEC pder;
      endpiece.momDeriv(mdir,tstraw,pder,dmvec);
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
    float cstart = pktraj.range().high() + tbuff_;
    pktraj.position(cstart,hend);
    float ltime = cstart + shmax_/pktraj.speed(cstart);
    pktraj.position(ltime,shmpos); // true position at shower-max
    // smear the x-y position by the transverse variance.
    lmeas.SetX(tr_.Gaus(shmpos.X(),twsig_));
    lmeas.SetY(tr_.Gaus(shmpos.Y(),twsig_));
    // set the z position to the sensor plane (end of the crystal)
    lmeas.SetZ(hend.Z()+clen_);
    // set the measurement time to correspond to the light propagation from showermax_, smeared by the resolution
    float tmeas = tr_.Gaus(ltime+(lmeas.Z()-shmpos.Z())/cprop_,ttsig_);
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
    for(int idir=0;idir<KInter::ndir;idir++){
      seed.momDeriv(KInter::MDir(idir),seed.range().mid(),pder,unit);
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
    ROOT::Math::SMatrix<float,3> bgrad;
    Vec3 pos,vel, dBdt;
    pktraj.position(htime,pos);
    pktraj.velocity(htime,vel);
    bfield_.fieldDeriv(pos,vel,dBdt);
    if(dBdt.R() != 0.0){
      auto const& back = pktraj.back();
      float tend = back.range().low();
      float tstep = 0.001*back.bnom().R()/dBdt.R(); // how far before BField changes by 1/1000.  Should be a parameter FIXME!
//      std::cout << "TStep " << tstep << std::endl;
      if(tstep < htime-tend){
	while(tend < htime-tstep){
	  tend += tstep;
	  Vec3 bf;
	  Vec4 pos; pos.SetE(tend);
	  Mom4 mom;
	  pktraj.momentum(tend,mom);
	  pktraj.position(pos);
	  bfield_.fieldVect(pos.Vect(),bf);
//	  std::cout << "BField " << bf << std::endl;
	  KTRAJ newend(pos,mom,pktraj.charge(),bf,TRange(tend,pktraj.range().high()));
	  pktraj.append(newend);
	}
      }
    }
  }

  template <class KTRAJ> void ToyMC<KTRAJ>::createTraj(PKTRAJ& pktraj) {
    // randomize the position and momentum
//    Vec4 torigin(tr_.Gaus(0.0,osig_), tr_.Gaus(0.0,osig_), tr_.Gaus(-zrange_/2.1,osig_),tr_.Gaus(0.0,osig_));
    Vec4 torigin(tr_.Gaus(0.0,osig_), tr_.Gaus(0.0,osig_), tr_.Gaus(0.0,osig_),tr_.Gaus(0.0,osig_));
    double tphi = tr_.Uniform(-M_PI,M_PI);
    double tcost = tr_.Uniform(ctmin_,ctmax_);
    double tsint = sqrt(1.0-tcost*tcost);
    Mom4 tmomv(mom_*tsint*cos(tphi),mom_*tsint*sin(tphi),mom_*tcost,simmass_);
    Vec3 bsim;
    bfield_.fieldVect(torigin.Vect(),bsim);
    KTRAJ ktraj(torigin,tmomv,icharge_,bsim);
    pktraj = PKTRAJ(ktraj);
    Vec3 vel; pktraj.velocity(0.0,vel);
    pktraj.setRange(TRange(-0.5*zrange_/vel.Z()-tbuff_,0.5*zrange_/vel.Z()+tbuff_));
  }

}


