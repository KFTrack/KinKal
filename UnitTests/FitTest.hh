//
// ToyMC test of fitting an KTRAJ-based Track
//
#include "MatEnv/MatDBInfo.hh"
#include "MatEnv/DetMaterial.hh"
#include "KinKal/ParticleTrajectory.hh"
#include "KinKal/Line.hh"
#include "KinKal/StrawHit.hh"
#include "KinKal/StrawMat.hh"
#include "KinKal/ScintHit.hh"
#include "KinKal/Constraint.hh"
#include "KinKal/BFieldMap.hh"
#include "KinKal/Vectors.hh"
#include "KinKal/Config.hh"
#include "KinKal/Measurement.hh"
#include "KinKal/Material.hh"
#include "KinKal/BFieldEffect.hh"
#include "KinKal/DetectorXing.hh"
#include "KinKal/Track.hh"
#include "UnitTests/ToyMC.hh"
#include "UnitTests/HitInfo.hh"
#include "UnitTests/MaterialInfo.hh"
#include "UnitTests/BFieldInfo.hh"
#include "UnitTests/ParticleTrajectoryInfo.hh"
#include "CLHEP/Units/PhysicalConstants.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <getopt.h>
#include <typeinfo>
#include <vector>
#include <cmath>
#include <ctime>
#include <chrono>
#include <cfenv>
#include <memory>
#include <cstdlib>
#include <cstring>

#include "TH1F.h"
#include "TTree.h"
#include "TBranch.h"
#include "TSystem.h"
#include "THelix.h"
#include "TPolyLine3D.h"
#include "TAxis3D.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TVector3.h"
#include "TPolyLine3D.h"
#include "TPolyMarker3D.h"
#include "TLegend.h"
#include "TGraph.h"
#include "TRandom3.h"
#include "TH2F.h"
#include "TDirectory.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "Math/VectorUtil.h"
#include <limits>

using namespace MatEnv;
using namespace KinKal;
using namespace std;
// avoid confusion with root
using KinKal::Line;
void print_usage() {
  printf("Usage: FitTest  --momentum f --simparticle i --fitparticle i--charge i --nhits i --hres f --seed i -maxniter i --deweight f --ambigdoca f --nevents i --simmat i--fitmat i --ttree i --Bz f --dBx f --dBy f --dBz f--Bgrad f --tolerance f--TFile c --PrintBad i --PrintDetail i --ScintHit i --bfcorr i --invert i --Schedule a --ssmear i --constrainpar i\n");
}

// utility function to compute transverse distance between 2 similar trajectories.  Also
// return the difference in time at the comparison point
template <class KTRAJ>
double dTraj(KTRAJ const& kt1, KTRAJ const& kt2, double t1, double& t2) {
  double dt = numeric_limits<float>::max();
  VEC3 pos1 = kt1.position3(t1);
  VEC3 dir1 = kt1.direction(t1);
  t2 = t1;
  VEC3 delta, v2;
  unsigned maxniter(100);
  unsigned niter(0);
  while(fabs(dt) > 1e-5 && niter < maxniter){
    v2 = kt2.velocity(t2);
    delta = kt2.position3(t2) - pos1;
    dt = delta.Dot(v2)/v2.Mag2();
    t2 -= dt;
    niter++;
  }
  if(niter >= maxniter) cout << "traj iteration not converged, dt = " << dt << endl;
  return (delta.Cross(dir1)).R();
}

template <class KTRAJ>
int FitTest(int argc, char **argv) {
  struct KTRAJPars{
    Float_t pars_[NParams()];
    static std::string leafnames() {
      std::string names;
      for(size_t ipar=0;ipar<NParams();ipar++){
	names +=KTRAJ::paramName(static_cast<typename KTRAJ::ParamIndex>(ipar)) + string("/f");
	if(ipar < NParams()-1)names += ":";
      }
      return names;
    }
  };

  using KKEFF = Effect<KTRAJ>;
  using KKHIT = Measurement<KTRAJ>;
  using KKMAT = Material<KTRAJ>;
  using KKBF = BFieldEffect<KTRAJ>;
  using KKEND = TrackEnd<KTRAJ>;
  using PKTRAJ = ParticleTrajectory<KTRAJ>;
  using DHIT = DetectorHit<KTRAJ>;
  using DHITPTR = std::shared_ptr<DHIT>;
  using DHITCOL = std::vector<DHITPTR>;
  using DXING = DetectorXing<KTRAJ>;
  using DXINGPTR = std::shared_ptr<DXING>;
  using DXINGCOL = std::vector<DXINGPTR>;
  using KKTRK = KinKal::Track<KTRAJ>;
  using KKCONFIGPTR = std::shared_ptr<Config>;
  using STRAWHIT = StrawHit<KTRAJ>;
  using STRAWHITPTR = std::shared_ptr<STRAWHIT>;
  using SCINTHIT = ScintHit<KTRAJ>;
  using SCINTHITPTR = std::shared_ptr<SCINTHIT>;
  using CONSTRAINT = Constraint<KTRAJ>;
  using CONSTRAINTPTR = std::shared_ptr<CONSTRAINT>;
  using Clock = std::chrono::high_resolution_clock;
  using PMASK = std::array<bool,NParams()>; // parameter mask

  // enable throw on FPE; not working with clang!
  fetestexcept(FE_ALL_EXCEPT );
  // parameters
  int opt;
  double mom(105.0);
  double masses[5]={0.511,105.66,139.57, 493.68, 938.0};
  int isimmass(0), ifitmass(0), icharge(-1);
  double simmass, fitmass;
  unsigned maxniter(5);
  double dwt(1.0e6);
  unsigned nevents(100);
  bool ttree(true), printbad(false);
  string tfname("FitTest.root"), sfile("Schedule.txt");
  int detail(Config::minimal), invert(0);
  double ambigdoca(-1.0);// minimum doca to set ambiguity, default sets for all hits
  Config::BFCorr bfcorr(Config::nocorr);
  bool fitmat(true);
  vector<double> sigmas = { 3.0, 3.0, 3.0, 3.0, 0.1, 3.0}; // base sigmas for parameter plots
  BFieldMap *BF(0);
  double Bgrad(0.0), dBx(0.0), dBy(0.0), dBz(0.0), Bz(1.0);
  double zrange(3000);
  double tol(0.1);
  int iseed(123421);
  int conspar(-1), iprint(-1);
  unsigned nhits(40);
  unsigned nsteps(200); // steps for traj comparison
  bool simmat(true), lighthit(true), seedsmear(true);

  static struct option long_options[] = {
    {"momentum",     required_argument, 0, 'm' },
    {"simparticle",     required_argument, 0, 'S'  },
    {"fitparticle",     required_argument, 0, 'F'  },
    {"charge",     required_argument, 0, 'q'  },
    {"seed",     required_argument, 0, 's'  },
    {"hres",     required_argument, 0, 'h'  },
    {"nhits",     required_argument, 0, 'n'  },
    {"maxniter",     required_argument, 0, 'i'  },
    {"deweight",     required_argument, 0, 'w'  },
    {"simmat",     required_argument, 0, 'b'  },
    {"fitmat",     required_argument, 0, 'f'  },
    {"ambigdoca",     required_argument, 0, 'd'  },
    {"nevents",     required_argument, 0, 'N'  },
    {"ttree",     required_argument, 0, 'r'  },
    {"tolerance",     required_argument, 0, 't'  },
    {"TFile",     required_argument, 0, 'T'  },
    {"dBx",     required_argument, 0, 'x'  },
    {"dBy",     required_argument, 0, 'y'  },
    {"dBz",     required_argument, 0, 'Z'  },
    {"Bz",     required_argument, 0, 'z'  },
    {"Bgrad",     required_argument, 0, 'g'  },
    {"PrintBad",     required_argument, 0, 'P'  },
    {"PrintDetail",     required_argument, 0, 'D'  },
    {"ScintHit",     required_argument, 0, 'L'  },
    {"bfcorr",     required_argument, 0, 'B'  },
    {"invert",     required_argument, 0, 'I'  },
    {"Schedule",     required_argument, 0, 'u'  },
    {"seedsmear",     required_argument, 0, 'M' },
    {"constrainpar",     required_argument, 0, 'c' },
    {"iprint",     required_argument, 0, 'p' },
    {NULL, 0,0,0}
  };

  int long_index =0;
  while ((opt = getopt_long_only(argc, argv,"",
	  long_options, &long_index )) != -1) {
    switch (opt) {
      case 'm' : mom = atof(optarg);
		 break;
      case 'S' : isimmass = atoi(optarg);
		 break;
      case 'F' : ifitmass = atoi(optarg);
		 break;
      case 'q' : icharge = atoi(optarg);
		 break;
      case 'z' : Bz = atof(optarg);
		 break;
      case 'n' : nhits = atoi(optarg);
		 break;
      case 's' : iseed = atoi(optarg);
		 break;
      case 'i' : maxniter = atoi(optarg);
		 break;
      case 'w' : dwt = atof(optarg);
		 break;
      case 'b' : simmat = atoi(optarg);
		 break;
      case 'f' : fitmat = atoi(optarg);
		 break;
      case 'L' : lighthit = atoi(optarg);
		 break;
      case 'B' : bfcorr = Config::BFCorr(atoi(optarg));
		 break;
      case 'r' : ttree = atoi(optarg);
		 break;
      case 'd' : ambigdoca = atof(optarg);
		 break;
      case 'M' : seedsmear = atoi(optarg);
		 break;
      case 'N' : nevents = atoi(optarg);
		 break;
      case 'x' : dBx = atof(optarg);
		 break;
      case 'y' : dBy = atof(optarg);
		 break;
      case 'Z' : dBz = atof(optarg);
		 break;
      case 'g' : Bgrad = atof(optarg);
		 break;
      case 'P' : printbad = atoi(optarg);
		 break;
      case 'p' : iprint = atoi(optarg);
		 break;
      case 'D' : detail = atoi(optarg);
		 break;
      case 'c' : conspar = atoi(optarg);
		 break;
      case 'I' : invert = atoi(optarg);
		 break;
      case 't' : tol = atof(optarg);
		 break;
      case 'T' : tfname = optarg;
		 break;
      case 'u' : sfile = optarg;
		 break;
      default: print_usage();
	       exit(EXIT_FAILURE);
    }
  }

  // construct BFieldMap
  VEC3 bnom;
  if(Bgrad != 0){
    BF = new GradientBFieldMap(Bz-0.5*Bgrad,Bz+0.5*Bgrad,-0.5*zrange,0.5*zrange);
    bnom = BF->fieldVect(VEC3(0.0,0.0,0.0));
  } else {
    VEC3 bsim(dBx,dBy,Bz+dBz);
    BF = new UniformBFieldMap(bsim);
    bnom = VEC3(0.0,0.0,Bz);
  }
  // create ToyMC
  simmass = masses[isimmass];
  fitmass = masses[ifitmass];
  KKTest::ToyMC<KTRAJ> toy(*BF, mom, icharge, zrange, iseed, nhits, simmat, lighthit, ambigdoca, simmass );
  // generate hits
  DHITCOL thits; // this program shares hit ownership with Track
  DXINGCOL dxings; // this program shares det xing ownership with Track
  PKTRAJ tptraj;
  toy.simulateParticle(tptraj, thits, dxings);
  toy.setSmearSeed(seedsmear);
  cout << "True initial " << tptraj.front() << endl;
//  cout << "vector of hit points " << thits.size() << endl;
//  cout << "True " << tptraj << endl;
  double startmom = tptraj.momentum(tptraj.range().begin());
  double endmom = tptraj.momentum(tptraj.range().end());
  VEC3 end, bend;
  bend = tptraj.front().direction(tptraj.range().end());
  end = tptraj.back().direction(tptraj.range().end());
  double angle = ROOT::Math::VectorUtil::Angle(bend,end);
  cout << "total momentum change = " << endmom-startmom << " total angle change = " << angle << endl;
  // create the fit seed by randomizing the parameters at the middle.  Overrwrite to use the fit BFieldMap
  auto const& midhel = tptraj.nearestPiece(0.0);
  auto seedmom = midhel.momentum4(0.0);
  seedmom.SetM(fitmass);
  // buffer the seed range
  TimeRange seedrange(tptraj.range().begin()-0.5,tptraj.range().end()+0.5);
  KTRAJ seedtraj(midhel.position4(0.0),seedmom,midhel.charge(),bnom,seedrange);
  if(invert) seedtraj.invertCT(); // for testing wrong propagation direction
  toy.createSeed(seedtraj);
  cout << "Seed params " << seedtraj.params().parameters() <<" covariance " << endl << seedtraj.params().covariance() << endl;
  // Create the Track from these hits
  //
  KKCONFIGPTR configptr = make_shared<Config>(*BF);
  configptr->dwt_ = dwt;
  configptr->maxniter_ = maxniter;
  configptr->bfcorr_ = bfcorr;
  configptr->addmat_ = fitmat;
  configptr->tol_ = tol;
  configptr->plevel_ = (Config::printLevel)detail;
  // read the schedule from the file
  string fullfile;
  if(strncmp(sfile.c_str(),"/",1) == 0) {
    fullfile = string(sfile);
  } else {
    if(const char* source = std::getenv("PACKAGE_SOURCE")){
      fullfile = string(source) + string("/UnitTests/") + string(sfile);
    } else {
      cout << "PACKAGE_SOURCE not defined" << endl;
      return -1;
    }
  }
  std::ifstream ifs (fullfile, std::ifstream::in);
  if ( (ifs.rdstate() & std::ifstream::failbit ) != 0 ){
    std::cerr << "Error opening " << fullfile << std::endl;
    return -1;
  }

  string line;
  unsigned nmiter(0);
  while (getline(ifs,line)){ 
    if(strncmp(line.c_str(),"#",1)!=0){
      istringstream ss(line);
      MetaIterConfig mconfig(ss);
      mconfig.miter_ = nmiter++;
      configptr->schedule_.push_back(mconfig);
    }
  }
  cout << *configptr << endl;
  // if requested, constrain a parameter
  PMASK mask = {false};
  if(conspar >= 0 && conspar < (int)NParams()){
    mask[conspar] = true;
    auto const& front = tptraj.front();
    // take the true parameters but the seed covariance
    Parameters cparams = front.params();
    cparams.covariance() = seedtraj.params().covariance();
    thits.push_back(std::make_shared<CONSTRAINT>(front.range().mid(),cparams,mask));
  }
// create and fit the track
  KKTRK kktrk(configptr,seedtraj,thits,dxings);
//  kktrk.print(cout,detail);
  TFile fitfile((KTRAJ::trajName() + tfname).c_str(),"RECREATE");
  // tree variables
  KTRAJPars ftpars_, btpars_, spars_, ffitpars_, ffiterrs_, bfitpars_, bfiterrs_;
  float chisq_, btmom_, ftmom_, ffmom_, bfmom_, ffmomerr_, bfmomerr_, chiprob_;
  float fft_,eft_;
  int ndof_, niter_, status_, igap_, nmeta_, nkkbf_, nkkhit_, nkkmat_;
  float maxgap_, avgap_;

  if(nevents <=0 ){
    // draw the fit result
    TCanvas* pttcan = new TCanvas("pttcan","PieceKTRAJ",1000,1000);
    auto const& fptraj = kktrk.fitTraj();
    unsigned np = fptraj.range().range()*fptraj.speed(fptraj.range().mid());
    TPolyLine3D* fitpl = new TPolyLine3D(np);
    fitpl->SetLineColor(kBlue);
    fitpl->SetLineStyle(kSolid);
    double ts = fptraj.range().range()/(np-1);
    for(unsigned ip=0;ip<np;ip++){
      double tp = fptraj.range().begin() + ip*ts;
      VEC3 ppos = fptraj.position3(tp);
      fitpl->SetPoint(ip,ppos.X(),ppos.Y(),ppos.Z());
    }
    fitpl->Draw();
// now draw the truth
    TPolyLine3D* ttpl = new TPolyLine3D(np);
    ttpl->SetLineColor(kGreen);
    ttpl->SetLineStyle(kDashDotted);
    ts = tptraj.range().range()/(np-1);
    for(unsigned ip=0;ip<np;ip++){
      double tp = tptraj.range().begin() + ip*ts;
      VEC3 ppos = tptraj.position3(tp);
      ttpl->SetPoint(ip,ppos.X(),ppos.Y(),ppos.Z());
    }
    ttpl->Draw();
    // draw the hits
    std::vector<TPolyLine3D*> htpls;
    for(auto const& thit : thits) {
      TPolyLine3D* line = new TPolyLine3D(2);
      VEC3 plow, phigh;
      STRAWHITPTR shptr = std::dynamic_pointer_cast<STRAWHIT> (thit); 
      SCINTHITPTR lhptr = std::dynamic_pointer_cast<SCINTHIT> (thit);
      if(shptr.use_count() > 0){
	auto const& tline = shptr->wire();
	plow = tline.position3(tline.range().begin());
	phigh = tline.position3(tline.range().end());
	line->SetLineColor(kRed);
      } else if (lhptr.use_count() > 0){
	auto const& tline = lhptr->sensorAxis();
	plow = tline.position3(tline.range().begin());
	phigh = tline.position3(tline.range().end());
	line->SetLineColor(kCyan);
      }
      line->SetPoint(0,plow.X(),plow.Y(), plow.Z());
      line->SetPoint(1,phigh.X(),phigh.Y(), phigh.Z());
      line->Draw();
      htpls.push_back(line);
    }

    // draw the origin and axes
    TAxis3D* rulers = new TAxis3D();
    rulers->GetXaxis()->SetAxisColor(kBlue);
    rulers->GetXaxis()->SetLabelColor(kBlue);
    rulers->GetYaxis()->SetAxisColor(kCyan);
    rulers->GetYaxis()->SetLabelColor(kCyan);
    rulers->GetZaxis()->SetAxisColor(kOrange);
    rulers->GetZaxis()->SetLabelColor(kOrange);
    rulers->Draw();
    pttcan->Write();

  } else {
    TTree* ftree(0);
    KKHIV hinfovec;
    KKBFIV bfinfovec;
    KKMIV minfovec;
    KTIV tinfovec;
    if(ttree){
      ftree = new TTree("fit","fit");
      ftree->Branch("ftpars.", &ftpars_,KTRAJPars::leafnames().c_str());
      ftree->Branch("btpars.", &btpars_,KTRAJPars::leafnames().c_str());
      ftree->Branch("spars.", &spars_,KTRAJPars::leafnames().c_str());
      ftree->Branch("ffpars.", &ffitpars_,KTRAJPars::leafnames().c_str());
      ftree->Branch("fferrs.", &ffiterrs_,KTRAJPars::leafnames().c_str());
      ftree->Branch("bfpars.", &bfitpars_,KTRAJPars::leafnames().c_str());
      ftree->Branch("bferrs.", &bfiterrs_,KTRAJPars::leafnames().c_str());
      ftree->Branch("chisq", &chisq_,"chisq/F");
      ftree->Branch("ndof", &ndof_,"ndof/I");
      ftree->Branch("nkkbf", &nkkbf_,"nkkbf/I");
      ftree->Branch("nkkmat", &nkkmat_,"nkkmat/I");
      ftree->Branch("nkkhit", &nkkhit_,"nkkhit/I");
      ftree->Branch("chiprob", &chiprob_,"chiprob/F");
      ftree->Branch("niter", &niter_,"niter/I");
      ftree->Branch("nmeta", &nmeta_,"nmeta/I");
      ftree->Branch("status", &status_,"status/I");
      ftree->Branch("ftmom", &ftmom_,"ftmom/F");
      ftree->Branch("btmom", &btmom_,"btmom/F");
      ftree->Branch("ffmom", &ffmom_,"ffmom/F");
      ftree->Branch("bfmom", &bfmom_,"bfmom/F");
      ftree->Branch("ffmomerr", &ffmomerr_,"ffmomerr/F");
      ftree->Branch("bfmomerr", &bfmomerr_,"bfmomerr/F");
      ftree->Branch("fft", &fft_,"fft/F");
      ftree->Branch("eft", &eft_,"eft/F");
      ftree->Branch("maxgap", &maxgap_,"maxgap/F");
      ftree->Branch("avgap", &avgap_,"avgap/F");
      ftree->Branch("igap", &igap_,"igap/I");
      ftree->Branch("hinfovec",&hinfovec);
      ftree->Branch("bfinfovec",&bfinfovec);
      ftree->Branch("minfovec",&minfovec);
      ftree->Branch("tinfovec",&tinfovec);
    }
    // now repeat this to gain statistics
    vector<TH1F*> fdp(NParams());
    vector<TH1F*> bdp(NParams());
    vector<TH1F*> fpull(NParams());
    vector<TH1F*> bpull(NParams());
    vector<TH1F*> fiterrh(NParams());
    TH1F* hniter = new TH1F("niter", "Total Iterations", 50,-0.5,49.5);
    TH1F* hnmeta = new TH1F("nmeta", "Meta Iterations", 10,-0.5,9.5);
    TH1F* hnfail = new TH1F("nfail", "Failed Iterations", 50,-0.5,49.5);
    TH1F* hndiv = new TH1F("ndiv", "Diverged Iterations", 50,-0.5,49.5);
    hnfail->SetLineColor(kRed);
    hndiv->SetLineColor(kOrange);
    TH1F* statush = new TH1F("statush", "Fit Status", 10,-0.5,9.5);
    TH1F* ndof = new TH1F("ndof", "N Degree of Freedom", 100,-0.5,99.5);
    TH1F* chisq = new TH1F("chisq", "Chisquared", 100,0,100);
    TH1F* chisqndof = new TH1F("chisqndof", "Chisquared per NDOF", 100,0,10.0);
    TH1F* chisqprob = new TH1F("chisqprob", "Chisquared probability", 100,0,1.0);
    TH1F* logchisqprob = new TH1F("logchisqprob", "Log10 of Chisquared probability", 100,-10,0.0);
    string htitle, hname;
    TH2F* corravg = new TH2F("corravg","Average correlation matrix magnitudes",NParams(),-0.5,NParams()-0.5,NParams(), -0.5,NParams()-0.5);
    TAxis* xax = corravg->GetXaxis();
    TAxis* yax = corravg->GetYaxis();
    double nsig(10.0);
    double pscale = nsig/sqrt(nhits);
    for(size_t ipar=0;ipar< NParams(); ipar++){
      auto tpar = static_cast<typename KTRAJ::ParamIndex>(ipar);
      hname = string("fd") + KTRAJ::paramName(tpar);
      htitle = string("Front #Delta ") + KTRAJ::paramTitle(tpar);
      fdp[ipar] = new TH1F(hname.c_str(),htitle.c_str(),100,-pscale*sigmas[ipar],pscale*sigmas[ipar]);
      hname = string("bd") + KTRAJ::paramName(tpar);
      htitle = string("Back #Delta ") + KTRAJ::paramTitle(tpar);
      bdp[ipar] = new TH1F(hname.c_str(),htitle.c_str(),100,-pscale*sigmas[ipar],pscale*sigmas[ipar]);
      hname = string("fp") + KTRAJ::paramName(tpar);
      htitle = string("Front Pull ") + KTRAJ::paramTitle(tpar);
      fpull[ipar] = new TH1F(hname.c_str(),htitle.c_str(),100,-nsig,nsig);
      hname = string("bp") + KTRAJ::paramName(tpar);
      htitle = string("Back Pull ") + KTRAJ::paramTitle(tpar);
      bpull[ipar] = new TH1F(hname.c_str(),htitle.c_str(),100,-nsig,nsig);
      hname = string("e") + KTRAJ::paramName(tpar);
      htitle = string("Error ") + KTRAJ::paramTitle(tpar);
      fiterrh[ipar] = new TH1F(hname.c_str(),htitle.c_str(),100,0.0,pscale*sigmas[ipar]);
      xax->SetBinLabel(ipar+1,KTRAJ::paramName(tpar).c_str());
      yax->SetBinLabel(ipar+1,KTRAJ::paramName(tpar).c_str());
    }
    TH1F* fmompull = new TH1F("fmompull","Front Momentum Pull;#Delta P/#sigma _{p}",100,-nsig,nsig);
    TH1F* bmompull = new TH1F("bmompull","Back Momentum Pull;#Delta P/#sigma _{p}",100,-nsig,nsig);
    double duration (0.0);
    unsigned nfail(0), ndiv(0);

    configptr->plevel_ = Config::none;
    for(unsigned ievent=0;ievent<nevents;ievent++){
    if( (ievent % iprint) == 0)
      cout << "event " << ievent << endl;
    // create a random true initial helix with hits and material interactions from this.  This also handles BFieldMap inhomogeneity truth tracking
      PKTRAJ tptraj;
      thits.clear();
      dxings.clear();
      toy.simulateParticle(tptraj,thits,dxings);
      double tmid = tptraj.range().mid();
      auto const& midhel = tptraj.nearestPiece(tmid);
      auto seedmom = midhel.momentum4(tmid);
      seedmom.SetM(fitmass);
      TimeRange seedrange(tptraj.range().begin()-0.5,tptraj.range().end()+0.5);
      KTRAJ seedtraj(midhel.position4(tmid),seedmom,midhel.charge(),bnom,seedrange);
      if(invert)seedtraj.invertCT();
      toy.createSeed(seedtraj);
  // if requested, constrain a parameter
      if(conspar >= 0 && conspar < (int)NParams()){
	auto const& front = tptraj.front();
	Parameters cparams = front.params();
	cparams.covariance() = seedtraj.params().covariance();
	thits.push_back(std::make_shared<CONSTRAINT>(front.range().mid(),cparams,mask));
      }
      auto start = Clock::now();
      KKTRK kktrk(configptr,seedtraj,thits,dxings);
      auto stop = Clock::now();
      duration += std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start).count();
      auto const& fstat = kktrk.fitStatus();
      if(fstat.status_ == FitStatus::failed)nfail++;
      if(fstat.status_ == FitStatus::diverged)ndiv++;
      niter_ = 0;
      for(auto const& fstat: kktrk.history()){
	if(fstat.status_ != FitStatus::unfit)niter_++;
      }
      // reset some fit parameters, to signal failed filts
      chiprob_ = -1.0;
      maxgap_ = avgap_ = -1;
      igap_ = -1;
      // fill effect information
      nkkbf_ = 0; nkkhit_ = 0; nkkmat_ = 0;
      // accumulate chisquared info
      chisq_ = fstat.chisq_;
      ndof_ = fstat.ndof_;
      niter_ = fstat.iter_;
      nmeta_ = fstat.miter_;
      status_ = fstat.status_;
      chiprob_ = fstat.prob_;
      hinfovec.clear();
      bfinfovec.clear();
      minfovec.clear();
      tinfovec.clear();
      statush->Fill(fstat.status_);
      if(fstat.usable()){
	// truth parameters, front and back
	double ttlow = tptraj.range().begin();
	double tthigh = tptraj.range().end();
	KTRAJ const& fttraj = tptraj.nearestPiece(ttlow);
	KTRAJ const& bttraj = tptraj.nearestPiece(tthigh);
	for(size_t ipar=0;ipar<6;ipar++){
	  spars_.pars_[ipar] = seedtraj.params().parameters()[ipar];
	  ftpars_.pars_[ipar] = fttraj.params().parameters()[ipar];
	  btpars_.pars_[ipar] = bttraj.params().parameters()[ipar];
	}
	ftmom_ = tptraj.momentum(ttlow);
	btmom_ = tptraj.momentum(tthigh);
	ndof->Fill(fstat.ndof_);
	chisq->Fill(fstat.chisq_);
	chisqndof->Fill(fstat.chisq_/fstat.ndof_);
	chisqprob->Fill(chiprob_);
	if(chiprob_ > 0.0) logchisqprob->Fill(log10(chiprob_));
	hniter->Fill(niter_);
	hnmeta->Fill(nmeta_);
	for(auto const& eff: kktrk.effects()) {
	  const KKHIT* kkhit = dynamic_cast<const KKHIT*>(eff.get());
	  const KKBF* kkbf = dynamic_cast<const KKBF*>(eff.get());
	  const KKMAT* kkmat = dynamic_cast<const KKMAT*>(eff.get());
	  if(kkhit != 0){
	    nkkhit_++;
	    HitInfo hinfo;
	    hinfo.active_ = kkhit->isActive();
	    hinfo.time_ = kkhit->time();
	    hinfo.fitchi_ = kkhit->fitChi();
	    const STRAWHIT* strawhit = dynamic_cast<const STRAWHIT*>(kkhit->detectorHit().get());
	    const SCINTHIT* scinthit = dynamic_cast<const SCINTHIT*>(kkhit->detectorHit().get());
	    const CONSTRAINT* constraint = dynamic_cast<const CONSTRAINT*>(kkhit->detectorHit().get());
	    if(strawhit != 0){
	      hinfo.type_ = HitInfo::straw;
	      hinfo.resid_ = strawhit->refResidual().value();
	      hinfo.residvar_ = strawhit->refResidual().variance();
	    } else if(scinthit != 0){
	      hinfo.type_ = HitInfo::scint;
	      hinfo.resid_ = scinthit->refResidual().value();
	      hinfo.residvar_ = scinthit->refResidual().variance();
	    } else if(constraint != 0){
	      hinfo.type_ = HitInfo::constraint;
	      hinfo.resid_ = constraint->chi(constraint->constraintParameters());
	      hinfo.residvar_ = 1.0;
	    } else {
	      hinfo.type_ = HitInfo::unknown;
	      hinfo.resid_ =  0.0;
	      hinfo.residvar_ = 1.0;
	    }
	    hinfovec.push_back(hinfo);
	  }
	  if(kkmat != 0){
	    nkkmat_++;
	    KinKal::MaterialInfo minfo;
	    minfo.time_ = kkmat->time();
	    minfo.active_ = kkmat->isActive();
	    minfo.nxing_ = kkmat->detXing().matXings().size();
	    std::array<double,3> dmom = {0.0,0.0,0.0}, momvar = {0.0,0.0,0.0};
	    kkmat->detXing().materialEffects(kkmat->refKTraj(),TimeDir::forwards, dmom, momvar);
	    minfo.dmomf_ = dmom[MomBasis::momdir_];
	    minfo.momvar_ = momvar[MomBasis::momdir_];
	    minfo.perpvar_ = momvar[MomBasis::perpdir_];
	    minfovec.push_back(minfo);
	  }
	  if(kkbf != 0){
	    nkkbf_++;
	    BFieldInfo bfinfo;
	    bfinfo.active_ = kkbf->isActive();
	    bfinfo.dp_ = kkbf->time();
	    bfinfo.dp_ = kkbf->deltaP().R();
	    bfinfo.range_ = kkbf->range().range();
	    bfinfovec.push_back(bfinfo);
	  }
	}

	// step through the fit traj and compare to the truth
	auto const& fptraj = kktrk.fitTraj();
	double dt = fptraj.range().range()/nsteps;
	for(unsigned istep=0;istep < nsteps;istep++){
	  double tstep = fptraj.range().begin()+dt*istep;
	  double ttrue;
	  double dperp = dTraj(fptraj,tptraj,tstep,ttrue);
	  ParticleTrajectoryInfo ktinfo;
	  ktinfo.time_ = tstep;
	  ktinfo.dperp_ = dperp;
	  ktinfo.dt_= tstep-ttrue;
	  tinfovec.push_back(ktinfo);
	}
	// compare parameters at the first traj of both true and fit
	// correct the true parameters in case the BFieldMap isn't nominal
	// correct the sampling time for the t0 difference
	double ftlow,fthigh;
	dTraj(tptraj,fptraj,ttlow,ftlow);
	dTraj(tptraj,fptraj,tthigh,fthigh);
	KTRAJ fftraj(fptraj.measurementState(ftlow),fptraj.charge(),tptraj.bnom(ttlow),fptraj.nearestPiece(ftlow).range());
	KTRAJ bftraj(fptraj.measurementState(fthigh),fptraj.charge(),tptraj.bnom(tthigh),fptraj.nearestPiece(fthigh).range());
	// fit parameters
	auto const& ffpars = fftraj.params();
	auto const& bfpars = bftraj.params();
	double maxgap, avgap;
	size_t igap;
	fptraj.gaps(maxgap, igap, avgap);
	maxgap_ = maxgap;
	avgap_ = avgap;
	igap_ = igap;

	Parameters ftpars, btpars;
	ftpars = fttraj.params();
	btpars = bttraj.params();

	// accumulate parameter difference and pull
	vector<double> fcerr(6,0.0), bcerr(6,0.0);
	for(size_t ipar=0;ipar< NParams(); ipar++){
	  fcerr[ipar] = sqrt(ffpars.covariance()[ipar][ipar]);
	  bcerr[ipar] = sqrt(bfpars.covariance()[ipar][ipar]);
	  fdp[ipar]->Fill(ffpars.parameters()[ipar]-ftpars.parameters()[ipar]);
	  bdp[ipar]->Fill(bfpars.parameters()[ipar]-btpars.parameters()[ipar]);
	  fpull[ipar]->Fill((ffpars.parameters()[ipar]-ftpars.parameters()[ipar])/fcerr[ipar]);
	  bpull[ipar]->Fill((bfpars.parameters()[ipar]-btpars.parameters()[ipar])/bcerr[ipar]);
	  fiterrh[ipar]->Fill(fcerr[ipar]);
	}
	// accumulate average correlation matrix
	auto const& cov = ffpars.covariance();
	//    auto cormat = cov;
	for(size_t ipar=0; ipar <NParams();ipar++){
	  for(size_t jpar=ipar;jpar < NParams(); jpar++){
	    double corr = cov[ipar][jpar]/(fcerr[ipar]*fcerr[jpar]);
	    //	cormat[ipar][jpar] = corr;
	    corravg->Fill(ipar,jpar,fabs(corr));
	  }
	}
	// extract fit parameters and errors
	for(size_t ipar=0;ipar<6;ipar++){
	  ffitpars_.pars_[ipar] = fftraj.params().parameters()[ipar];
	  bfitpars_.pars_[ipar] = bftraj.params().parameters()[ipar];
	  ffiterrs_.pars_[ipar] = sqrt(fftraj.params().covariance()(ipar,ipar));
	  bfiterrs_.pars_[ipar] = sqrt(bftraj.params().covariance()(ipar,ipar));
	}
	ffmom_ = fptraj.momentum(ftlow);
	bfmom_ = fptraj.momentum(fthigh);
	ffmomerr_ = sqrt(fptraj.momentumVar(ftlow));
	bfmomerr_ = sqrt(fptraj.momentumVar(fthigh));
	fft_ = fptraj.range().begin();
	eft_ = fptraj.range().end();
	fmompull->Fill((ffmom_-ftmom_)/ffmomerr_);
	bmompull->Fill((bfmom_-btmom_)/bfmomerr_);
	// state space parameter difference and errors
	//	ParticleStateMeasurement tslow = tptraj.state(tlow);
	//	ParticleStateMeasurement tshigh = tptraj.state(thigh);
	//	ParticleStateMeasurement slow = fptraj.measurementState(tlow);
	//	ParticleStateMeasurement shigh = fptraj.measurementState(thigh);

	// test
      } else if(printbad){
	cout << "Bad Fit event " << ievent << " status " << kktrk.fitStatus() << endl;
	cout << "True Traj " << tptraj << endl;
	cout << "Seed Traj " << seedtraj << endl;
	kktrk.print(cout,detail);
      }
      if(ttree)ftree->Fill();
    }
    cout << nfail << " Failed fits and " << ndiv << " Diverged fits " << endl;
    hnfail->Fill(nfail);
    hndiv->Fill(ndiv);
    cout <<"Time/fit = " << duration/double(nevents) << " Nanoseconds " << endl;
    // fill canvases
    TCanvas* fdpcan = new TCanvas("fdpcan","fdpcan",800,600);
    fdpcan->Divide(3,2);
    for(size_t ipar=0;ipar<NParams();++ipar){
      fdpcan->cd(ipar+1);
      fdp[ipar]->Fit("gaus","q");
    }
    fdpcan->Write();
    TCanvas* bdpcan = new TCanvas("bdpcan","bdpcan",800,600);
    bdpcan->Divide(3,2);
    for(size_t ipar=0;ipar<NParams();++ipar){
      bdpcan->cd(ipar+1);
      bdp[ipar]->Fit("gaus","q");
    }
    bdpcan->Write();
    TCanvas* fpullcan = new TCanvas("fpullcan","fpullcan",800,600);
    fpullcan->Divide(3,3);
    for(size_t ipar=0;ipar<NParams();++ipar){
      fpullcan->cd(ipar+1);
      fpull[ipar]->Fit("gaus","q");
    }
    fpullcan->cd(NParams()+1);
    fmompull->Fit("gaus","q");
    fpullcan->Write();
    TCanvas* bpullcan = new TCanvas("bpullcan","bpullcan",800,600);
    bpullcan->Divide(3,3);
    for(size_t ipar=0;ipar<NParams();++ipar){
      bpullcan->cd(ipar+1);
      bpull[ipar]->Fit("gaus","q");
    }
    bpullcan->cd(NParams()+1);
    bmompull->Fit("gaus","q");
    bpullcan->Write();
    TCanvas* perrcan = new TCanvas("perrcan","perrcan",800,600);
    perrcan->Divide(3,2);
    for(size_t ipar=0;ipar<NParams();++ipar){
      perrcan->cd(ipar+1);
      fiterrh[ipar]->Draw();
    }
    perrcan->Write();
    TCanvas* corrcan = new TCanvas("corrcan","corrcan",600,600);
    corrcan->Divide(1,1);
    corrcan->cd(1);
    corravg->Scale(1.0/double(nevents));
    corravg->SetStats(0);
    gPad->SetLogz();
    corravg->Draw("colorztext0");
    corrcan->Write();

    TCanvas* statuscan = new TCanvas("statuscan","statuscan",800,600);
    statuscan->Divide(3,2);
    statuscan->cd(1);
    statush->Draw();
    statuscan->cd(2);
    hniter->Draw();
    hnfail->Draw("same");
    hndiv->Draw("same");
    statuscan->cd(3);
    ndof->Draw();
    statuscan->cd(4);
    chisqndof->Draw();
    statuscan->cd(5);
    chisqprob->Draw();
    statuscan->cd(6);
    logchisqprob->Draw();
    statuscan->Write();
  }
  fitfile.Write();
  fitfile.Close();

  if (kktrk.fitStatus().status_ != KinKal::FitStatus::converged)
  {
    exit(1);
  }

  exit(EXIT_SUCCESS);
}
