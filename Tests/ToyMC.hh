#ifndef KinKal_ToyMC_hh
#define KinKal_ToyMC_hh
//
//  Toy MC for fit and hit testing
//
#include <memory>
#include "TRandom3.h"
#include "KinKal/MatEnv/MatDBInfo.hh"
#include "KinKal/MatEnv/DetMaterial.hh"
#include "KinKal/MatEnv/SimpleFileFinder.hh"
#include "KinKal/Examples/CaloDistanceToTime.hh"
#include "KinKal/Trajectory/SensorLine.hh"
#include "KinKal/Trajectory/ParticleTrajectory.hh"
#include "KinKal/Trajectory/PiecewiseClosestApproach.hh"
#include "KinKal/Examples/SimpleWireHit.hh"
#include "KinKal/Detector/StrawXing.hh"
#include "KinKal/Detector/StrawMaterial.hh"
#include "KinKal/Examples/ScintHit.hh"
#include "KinKal/General/BFieldMap.hh"
#include "KinKal/General/Vectors.hh"
#include "KinKal/General/PhysicalConstants.h"

namespace KKTest {
  using namespace KinKal;
  template <class KTRAJ> class ToyMC {
    public:
      using PTRAJ = ParticleTrajectory<KTRAJ>;
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
      using PCA = PiecewiseClosestApproach<KTRAJ,SensorLine>;
      // create from aseed
      ToyMC(BFieldMap const& bfield, double mom, int icharge, double zrange, int iseed, unsigned nhits, bool simmat, bool scinthit, double ambigdoca ,double simmass) :
        bfield_(bfield), matdb_(sfinder_,MatEnv::DetMaterial::moyalmean), // use the moyal based eloss model
        mom_(mom), icharge_(icharge),
        tr_(iseed), nhits_(nhits), simmat_(simmat), scinthit_(scinthit), ambigdoca_(ambigdoca), simmass_(simmass),
        sprop_(0.8*CLHEP::c_light), sdrift_(0.065),
        zrange_(zrange), rstraw_(2.5), rwire_(0.025), wthick_(0.015), wlen_(1000.0), sigt_(3.0), sigtot_(7.0), ineff_(0.05),
        scitsig_(0.5), shPosSig_(10.0), shmax_(40.0), cz_(0.5*zrange_+50), clen_(200.0), cprop_(0.8*CLHEP::c_light),
        osig_(10.0), ctmin_(0.5), ctmax_(0.8), tol_(1e-5), tprec_(1e-8), t0off_(700.0),
        smat_(matdb_,rstraw_, wthick_, 3*wthick_, rwire_), miconfig_(0.0) {
          miconfig_.addUpdater(std::any(StrawXingConfig(1.0e6,1.0e6,1.0e6,false))); // updater to force exact straw xing material calculation
        }

      // generate a straw at the given time.  direction and drift distance are random
      SensorLine generateStraw(PTRAJ const& traj, double htime);
      // create a seed by randomizing the parameters
      void createSeed(KTRAJ& seed,DVEC const& sigmas, double seedsmear);
      void extendTraj(PTRAJ& ptraj,double htime);
      void createTraj(PTRAJ& ptraj);
      void createScintHit(PTRAJ& ptraj, HITCOL& thits);
      void simulateParticle(PTRAJ& ptraj,HITCOL& thits, EXINGCOL& dxings, bool addmat=true);
      double createStrawMaterial(PTRAJ& ptraj, EXING* sxing);
      // set functions, for special purposes
      void setInefficiency(double ineff) { ineff_ = ineff; }
      void setTolerance(double tol) { tol_ = tol; }
      // accessors
      double shVar() const {return sigt_*sigt_;}
      double chVar() const {return scitsig_*scitsig_;}
      double rStraw() const { return rstraw_; }
      double zRange() const { return zrange_; }
      double strawRadius() const { return rstraw_; }
      StrawMaterial const& strawMaterial() const { return smat_; }

    private:
      BFieldMap const& bfield_;
      MatEnv::SimpleFileFinder sfinder_;
      MatEnv::MatDBInfo matdb_;
      double mom_;
      int icharge_;
      TRandom3 tr_; // random number generator
      unsigned nhits_; // number of hits to simulate
      bool simmat_, scinthit_;
      double ambigdoca_, simmass_;
      double sprop_; // propagation speed along straw
      double sdrift_; // drift speed inside straw
      double zrange_, rstraw_; // tracker dimension
      double rwire_, wthick_, wlen_; // wire radius, thickness, length
      double sigt_; // drift time resolution in ns
      double sigtot_; // TOT drift time resolution (ns)
      double ineff_; // hit inefficiency
                     // time hit parameters
      double scitsig_, shPosSig_, shmax_, cz_, clen_, cprop_;
      double osig_, ctmin_, ctmax_;
      double tol_; // tolerance on momentum accuracy due to BField effects
      double tprec_; // time precision on TCA
      double t0off_; // t0 offset
      StrawMaterial smat_; // straw material
      MetaIterConfig miconfig_; // configuration used when calculating initial effects

  };

  template <class KTRAJ> SensorLine ToyMC<KTRAJ>::generateStraw(PTRAJ const& traj, double htime) {
    // start with the true helix position at this time
    auto hpos = traj.position4(htime);
    auto tdir = traj.direction(htime);
    // generate a random azimuth direction for the straw
    double azimuth = tr_.Uniform(-M_PI,M_PI);
    VEC3 sdir(cos(azimuth),sin(azimuth),0.0);
    // generate a random drift perp to this and the trajectory
    double rdrift = tr_.Uniform(-rstraw_,rstraw_);
    VEC3 driftdir = (sdir.Cross(tdir)).Unit();
    VEC3 dpos = hpos.Vect() + rdrift*driftdir;
    //  cout << "Generating hit at position " << dpos << endl;
    double dprop = tr_.Uniform(0.0,wlen_);
    VEC3 mpos = dpos + sdir*dprop;
    VEC3 vprop = sdir*sprop_;
    // measured time is after propagation and drift
    double tmeas = htime + dprop/sprop_ + fabs(rdrift)/sdrift_;
    // smear measurement time
    tmeas = tr_.Gaus(tmeas,sigt_);
    // construct the trajectory for this hit
    return SensorLine(mpos,tmeas,vprop,wlen_);
  }

  template <class KTRAJ> void ToyMC<KTRAJ>::simulateParticle(PTRAJ& ptraj,HITCOL& thits, EXINGCOL& dxings, bool addmat) {
    // create the seed first
    createTraj(ptraj);
    // divide time range
    double dt(0.0);
    if(nhits_ > 0) dt = (ptraj.range().range())/(nhits_);
    VEC3 bsim;
    // create the hits (and associated materials)
    for(size_t ihit=0; ihit<nhits_; ihit++){
      double htime = ptraj.range().begin() + ihit*dt;
      // extend the trajectory in the BFieldMap to this time
      extendTraj(ptraj,htime);
      // create the hit at this time
      auto tline = generateStraw(ptraj,htime);
      CAHint tphint(htime,htime);
      PCA tp(ptraj,tline,tphint,tprec_);
      //      std::cout << "doca " << tp.doca() << " sensor TOCA " << tp.sensorToca() - fabs(tp.doca())/sdrift_ << " particle TOCA " << tp.particleToca() << " hit time " << htime << std::endl;
      if(tr_.Uniform(0.0,1.0) > ineff_){
        WireHitState::State ambig(WireHitState::null);
        if(fabs(tp.doca())> ambigdoca_) ambig = tp.doca() < 0 ? WireHitState::left : WireHitState::right;
        WireHitState whstate(ambig);
        double mindoca = std::min(ambigdoca_,rstraw_);
        // true TOT is the drift time
        double tot = tr_.Gaus(tp.deltaT(),sigtot_);
        thits.push_back(std::make_shared<WIREHIT>(bfield_, tp, whstate, mindoca, sdrift_, sigt_*sigt_, tot, sigtot_*sigtot_,  rstraw_, ihit));
      }
      // compute material effects and change trajectory accordingly
      auto xing = std::make_shared<STRAWXING>(tp,smat_);
      if(addmat) dxings.push_back(xing);
      if(simmat_){
        double defrac = createStrawMaterial(ptraj, xing.get());
        // terminate if there is catastrophic energy loss
        if(fabs(defrac) > 0.1)break;
      }
    }
    if(scinthit_ && tr_.Uniform(0.0,1.0) > ineff_){
      createScintHit(ptraj,thits);
    }
    // extend if there are no hits
    if(thits.size() == 0) extendTraj(ptraj,ptraj.range().end());
  }

  template <class KTRAJ> double ToyMC<KTRAJ>::createStrawMaterial(PTRAJ& ptraj, EXING* sxing) {
    double desum = 0.0;
    double tstraw = sxing->time();
    auto const& endtraj = ptraj.nearestTraj(tstraw);
    auto const& endpiece = *endtraj;
    sxing->updateReference(ptraj);
    sxing->updateState(miconfig_,true);
    double mom = endpiece.momentum(tstraw);
    auto endmom = endpiece.momentum4(tstraw);
    auto endpos = endpiece.position4(tstraw);
    std::array<double,3> dmom {0.0,0.0,0.0}, momvar {0.0,0.0,0.0};
    sxing->materialEffects(dmom, momvar);
    for(int idir=0;idir<=MomBasis::phidir_; idir++) {
      auto mdir = static_cast<MomBasis::Direction>(idir);
      double momsig = sqrt(momvar[idir]);
      double dm;
      // generate a random effect given this variance and mean
      switch( mdir ) {
        case KinKal::MomBasis::perpdir_: case KinKal::MomBasis::phidir_ :
          dm = tr_.Gaus(dmom[idir],momsig);
          break;
        case KinKal::MomBasis::momdir_ :
          dm = std::min(0.0,tr_.Gaus(dmom[idir],momsig));
          desum += dm;
          break;
        default:
          throw std::invalid_argument("Invalid direction");
      }
      auto dmvec = endpiece.direction(tstraw,mdir);
      dmvec *= dm*mom;
      endmom.SetCoordinates(endmom.Px()+dmvec.X(), endmom.Py()+dmvec.Y(), endmom.Pz()+dmvec.Z(),endmom.M());
    }
    // generate a new piece and append
    VEC3 bnom = bfield_.fieldVect(endpos.Vect());
    KTRAJ newend(endpos,endmom,endpiece.charge(),bnom,TimeRange(tstraw,ptraj.range().end()));
    //      newend.print(cout,1);
    ptraj.append(newend);
    return desum/mom;
  }

  template <class KTRAJ> void ToyMC<KTRAJ>::createScintHit(PTRAJ& ptraj, HITCOL& thits) {
    // create a ScintHit at the end, axis parallel to z
    // first, find the position at showermax.
    VEC3 shmax,shmeas;
    double tend = thits.back()->time();
    // extend to the calorimeter z
    VEC3 pvel = ptraj.velocity(tend);
    VEC3 ppos = ptraj.position3(tend);
    double shstart = tend + (cz_-ppos.Z())/pvel.Z();
    // extend the trajectory to here
    extendTraj(ptraj,shstart);
    pvel = ptraj.velocity(shstart);
    // compute time at showermax
    double shmaxtime = shstart + shmax_/pvel.R();
    shmax = ptraj.position3(shmaxtime); // true position at shower-max
    // Compute measurement position: smear the x-y position by the transverse variance.
    shmeas.SetX(tr_.Gaus(shmax.X(),shPosSig_));
    shmeas.SetY(tr_.Gaus(shmax.Y(),shPosSig_));
    // set the z position to the sensor plane (forward end of the crystal)
    shmeas.SetZ(cz_+clen_);
    // set the measurement time to correspond to the light propagation from showermax_, smeared by the time resolution
    double tmeas = tr_.Gaus(shmaxtime+(shmeas.Z() - shmax.Z())/cprop_,scitsig_);
    // create the ttraj for the light propagation
    VEC3 lvel(0.0,0.0,cprop_);
    SensorLine lline(shmeas, tmeas, lvel, clen_);
    // then create the hit and add it; the hit has no material
    CAHint tphint(tmeas,tmeas);
    PCA pca(ptraj,lline,tphint,tprec_);
    thits.push_back(std::make_shared<SCINTHIT>(pca, scitsig_*scitsig_, shPosSig_*shPosSig_));
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

  template <class KTRAJ> void ToyMC<KTRAJ>::extendTraj(PTRAJ& ptraj,double htime) {
    if(htime > ptraj.range().end()){
      VEC3 pos,vel, dBdt;
      pos = ptraj.position3(htime);
      vel = ptraj.velocity(htime);
      dBdt = bfield_.fieldDeriv(pos,vel);
      //    std::cout << "end time " << ptraj.back().range().begin() << " hit time " << htime << std::endl;
      if(dBdt.R() != 0.0){
        double tbeg = ptraj.range().end();
        while(tbeg < htime) {
          double tend = std::min(tbeg + bfield_.rangeInTolerance(ptraj.back(),tbeg,tol_),htime);
          double tmid = 0.5*(tbeg+tend);
          auto bf = bfield_.fieldVect(ptraj.position3(tmid));
          auto pos = ptraj.position4(tend);
          auto mom =  ptraj.momentum4(tend);
          TimeRange prange(tbeg,tend);
          KTRAJ newend(pos,mom,ptraj.charge(),bf,prange);
          // make sure phi0 stays continuous
          newend.syncPhi0(ptraj.back());
          ptraj.append(newend);
          tbeg = tend;
        }
      }
    }
  }

  template <class KTRAJ> void ToyMC<KTRAJ>::createTraj(PTRAJ& ptraj) {
    // randomize the position and momentum
    double tphi = tr_.Uniform(-M_PI,M_PI);
    double tcost = tr_.Uniform(ctmin_,ctmax_);
    double tsint = sqrt(1.0-tcost*tcost);
    MOM4 tmomv(mom_*tsint*cos(tphi),mom_*tsint*sin(tphi),mom_*tcost,simmass_);
    double tmax = fabs(zrange_/(CLHEP::c_light*tcost));
    double tbeg = tr_.Uniform(t0off_-tmax,t0off_+tmax);
    VEC4 torigin(tr_.Gaus(0.0,osig_), tr_.Gaus(0.0,osig_), tr_.Gaus(-0.5*zrange_,osig_),tbeg);
    VEC3 bsim = bfield_.fieldVect(torigin.Vect());
    KTRAJ ktraj(torigin,tmomv,icharge_,bsim,TimeRange(tbeg,tbeg+tmax));
    // set initial range
    double tend = tbeg + bfield_.rangeInTolerance(ktraj,tbeg,tol_);
    ktraj.setRange(TimeRange(tbeg,tend));
    ptraj.append(ktraj);
  }
}
#endif
