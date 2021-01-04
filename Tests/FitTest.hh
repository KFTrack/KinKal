//
// ToyMC test of fitting an KTRAJ-based Track
//
#include "KinKal/MatEnv/MatDBInfo.hh"
#include "KinKal/MatEnv/DetMaterial.hh"
#include "KinKal/General/Vectors.hh"
#include "KinKal/Trajectory/ParticleTrajectory.hh"
#include "KinKal/Trajectory/Line.hh"
#include "KinKal/Tests/SimpleWireHit.hh"
#include "KinKal/Detector/StrawMat.hh"
#include "KinKal/Detector/ParameterHit.hh"
#include "KinKal/Tests/ScintHit.hh"
#include "KinKal/Detector/ElementXing.hh"
#include "KinKal/Detector/BFieldMap.hh"
#include "KinKal/Fit/Config.hh"
#include "KinKal/Fit/HitConstraint.hh"
#include "KinKal/Fit/Material.hh"
#include "KinKal/Fit/BFieldEffect.hh"
#include "KinKal/Fit/Track.hh"
#include "KinKal/Tests/ToyMC.hh"
#include "KinKal/Tests/HitInfo.hh"
#include "KinKal/Tests/MaterialInfo.hh"
#include "KinKal/Tests/BFieldInfo.hh"
#include "KinKal/Tests/ParticleTrajectoryInfo.hh"
#include "KinKal/General/PhysicalConstants.h"

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
#include "TFitResult.h"
#include "Math/VectorUtil.h"
#include <limits>

using namespace MatEnv;
using namespace KinKal;
using namespace std;
// avoid confusion with root
using KinKal::Line;
void print_usage() {
  printf("Usage: FitTest  --momentum f --simparticle i --fitparticle i--charge i --nhits i --hres f --seed i --maxniter i --deweight f --ambigdoca f --nevents i --simmat i--fitmat i --ttree i --Bz f --dBx f --dBy f --dBz f--Bgrad f --tolerance f --TFilesuffix c --PrintBad i --PrintDetail i --ScintHit i --nulltime i--bfcorr i --invert i --Schedule a --ssmear i --constrainpar i\n");
}

// utility function to compute transverse distance between 2 similar trajectories.  Also
// return the difference in time at the comparison point
template <class KTRAJ>
double dTraj(KTRAJ const& kt1, KTRAJ const& kt2, double t1, double& t2) {
  double dt = numeric_limits<float>::max();
  VEC3 pos1 = kt1.position3(t1);
  VEC3 dir1 = kt1.direction(t1);
  t2 = t1;
  unsigned maxniter(100);
  unsigned niter(0);
  VEC3 pos2;
//  double delta;
  while(fabs(dt) > 1e-5 && niter < maxniter){
    VEC3 v2 = kt2.velocity(t2);
    pos2 = kt2.position3(t2);
// test
//    delta = pos2.Z() - pos1.Z();
//    dt = delta/v2.Z();
    VEC3 dp = kt2.position3(t2) - pos1;
    dt = dp.Dot(v2)/v2.Mag2();

t2 -= dt;
    niter++;
  }
  if(niter >= maxniter) cout << "traj iteration not converged, dt = " << dt << endl;
  return ((pos2-pos1).Cross(dir1)).R();
}

template <class KTRAJ>
int FitTest(int argc, char *argv[],const vector<double>& sigmas) {
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
  using KKHIT = HitConstraint<KTRAJ>;
  using KKMAT = Material<KTRAJ>;
  using KKBF = BFieldEffect<KTRAJ>;
  using KKEND = TrackEnd<KTRAJ>;
  using PKTRAJ = ParticleTrajectory<KTRAJ>;
  using MEAS = Hit<KTRAJ>;
  using MEASPTR = std::shared_ptr<MEAS>;
  using MEASCOL = std::vector<MEASPTR>;
  using EXING = ElementXing<KTRAJ>;
  using EXINGPTR = std::shared_ptr<EXING>;
  using EXINGCOL = std::vector<EXINGPTR>;
  using KKTRK = KinKal::Track<KTRAJ>;
  using KKCONFIGPTR = std::shared_ptr<Config>;
  using STRAWHIT = WireHit<KTRAJ>;
  using STRAWHITPTR = std::shared_ptr<STRAWHIT>;
  using SCINTHIT = ScintHit<KTRAJ>;
  using SCINTHITPTR = std::shared_ptr<SCINTHIT>;
  using PARHIT = ParameterHit<KTRAJ>;
  using PARHITPTR = std::shared_ptr<PARHIT>;
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
  unsigned maxniter(10);
  double dwt(1.0e6);
  unsigned nevents(100);
  bool ttree(true), printbad(false);
  string tfname(""), sfile("Schedule.txt");
  int detail(Config::minimal), invert(0);
  double ambigdoca(-1.0);// minimum doca to set ambiguity, default sets for all hits
  Config::BFCorr bfcorr(Config::nocorr);
  bool fitmat(true);
  BFieldMap *BF(0);
  double Bgrad(0.0), dBx(0.0), dBy(0.0), dBz(0.0), Bz(1.0);
  double zrange(3000);
  double tol(0.1);
  double seedposvar(3.0); //position variance for seed smearing
  int iseed(123421);
  int conspar(-1), iprint(-1);
  unsigned nhits(40);
  unsigned nsteps(200); // steps for traj comparison
  bool simmat(true), lighthit(true), seedsmear(true), nulltime(true);
  int retval(EXIT_SUCCESS);

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
    {"TFilesuffix",     required_argument, 0, 'T'  },
    {"dBx",     required_argument, 0, 'x'  },
    {"dBy",     required_argument, 0, 'y'  },
    {"dBz",     required_argument, 0, 'Z'  },
    {"Bz",     required_argument, 0, 'z'  },
    {"Bgrad",     required_argument, 0, 'g'  },
    {"PrintBad",     required_argument, 0, 'P'  },
    {"PrintDetail",     required_argument, 0, 'D'  },
    {"ScintHit",     required_argument, 0, 'L'  },
    {"NullTime",     required_argument, 0, 'v'  },
    {"bfcorr",     required_argument, 0, 'B'  },
    {"invert",     required_argument, 0, 'I'  },
    {"Schedule",     required_argument, 0, 'u'  },
    {"seedsmear",     required_argument, 0, 'M' },
    {"constrainpar",     required_argument, 0, 'c' },
    {"iprint",     required_argument, 0, 'p' },
    {NULL, 0,0,0}
  };

  int long_index =0;
  while ((opt = getopt_long(argc, argv,"",
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
      case 'v' : nulltime = atoi(optarg);
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
  KKTest::ToyMC<KTRAJ> toy(*BF, mom, icharge, zrange, iseed, nhits, simmat, lighthit, nulltime, ambigdoca, simmass );
  // generate hits
  MEASCOL thits; // this program shares hit ownership with Track
  EXINGCOL dxings; // this program shares det xing ownership with Track
  PKTRAJ tptraj;
  toy.simulateParticle(tptraj, thits, dxings);
  if(nevents < 0)cout << "True initial " << tptraj.front() << endl;
//  cout << "vector of hit points " << thits.size() << endl;
//  cout << "True " << tptraj << endl;
  double startmom = tptraj.momentum(tptraj.range().begin());
  double endmom = tptraj.momentum(tptraj.range().end());
  VEC3 end, bend;
  bend = tptraj.front().direction(tptraj.range().end());
  end = tptraj.back().direction(tptraj.range().end());
  double angle = ROOT::Math::VectorUtil::Angle(bend,end);
  if(nevents < 0)cout << "total momentum change = " << endmom-startmom << " total angle change = " << angle << endl;
  // create the fit seed by randomizing the parameters at the middle.  Overrwrite to use the fit BFieldMap
  auto const& midhel = tptraj.nearestPiece(0.0);
  auto seedmom = midhel.momentum4(0.0);
  seedmom.SetM(fitmass);
  // buffer the seed range
  TimeRange seedrange(tptraj.range().begin()-0.5,tptraj.range().end()+0.5);
  KTRAJ seedtraj(midhel.position4(0.0),seedmom,midhel.charge(),bnom,seedrange);
  if(invert) seedtraj.invertCT(); // for testing wrong propagation direction
  toy.createSeed(seedtraj,sigmas[NParams()],seedposvar,seedsmear);
  if(nevents < 0)cout << "Seed Traj " << seedtraj << endl;
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
      fullfile = string(source) + string("/Tests/") + string(sfile);
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
  if(nevents < 0)cout << *configptr << endl;
  // if requested, constrain a parameter
  PMASK mask = {false};
  if(conspar >= 0 && conspar < (int)NParams()){
    cout << "Constraining parameter " << conspar << endl;
    mask[conspar] = true;
    auto const& front = tptraj.front();
    // take the true parameters but the seed covariance
    Parameters cparams = front.params();
    cparams.covariance() = seedtraj.params().covariance();
    thits.push_back(std::make_shared<PARHIT>(front.range().begin(),cparams,mask));
  }
// create and fit the track
  KKTRK kktrk(configptr,seedtraj,thits,dxings);
//  kktrk.print(cout,detail);
  TFile fitfile((KTRAJ::trajName() + string("FitTest") + tfname + string(".root")).c_str(),"RECREATE");
  // tree variables
  KTRAJPars ftpars_, mtpars_, btpars_, spars_, ffitpars_, ffiterrs_, mfitpars_, mfiterrs_, bfitpars_, bfiterrs_;
  float chisq_, btmom_, mtmom_, ftmom_, ffmom_, mfmom_, bfmom_, ffmomerr_, mfmomerr_, bfmomerr_, chiprob_;
  float fft_,mft_, bft_;
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
    if (kktrk.fitStatus().status_ != KinKal::Status::converged)retval = -1;
  } else {
    TTree* ftree(0);
    KKHIV hinfovec;
    KKBFIV bfinfovec;
    KKMIV minfovec;
    KTIV tinfovec;
    if(ttree){
      ftree = new TTree("fit","fit");
      ftree->Branch("ftpars.", &ftpars_,KTRAJPars::leafnames().c_str());
      ftree->Branch("mtpars.", &mtpars_,KTRAJPars::leafnames().c_str());
      ftree->Branch("btpars.", &btpars_,KTRAJPars::leafnames().c_str());
      ftree->Branch("spars.", &spars_,KTRAJPars::leafnames().c_str());
      ftree->Branch("ffpars.", &ffitpars_,KTRAJPars::leafnames().c_str());
      ftree->Branch("fferrs.", &ffiterrs_,KTRAJPars::leafnames().c_str());
      ftree->Branch("mfpars.", &mfitpars_,KTRAJPars::leafnames().c_str());
      ftree->Branch("mferrs.", &mfiterrs_,KTRAJPars::leafnames().c_str());
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
      ftree->Branch("mtmom", &mtmom_,"mtmom/F");
      ftree->Branch("btmom", &btmom_,"btmom/F");
      ftree->Branch("ffmom", &ffmom_,"ffmom/F");
      ftree->Branch("mfmom", &mfmom_,"mfmom/F");
      ftree->Branch("bfmom", &bfmom_,"bfmom/F");
      ftree->Branch("ffmomerr", &ffmomerr_,"ffmomerr/F");
      ftree->Branch("mfmomerr", &mfmomerr_,"mfmomerr/F");
      ftree->Branch("bfmomerr", &bfmomerr_,"bfmomerr/F");
      ftree->Branch("fft", &fft_,"fft/F");
      ftree->Branch("mft", &mft_,"mft/F");
      ftree->Branch("bft", &bft_,"bft/F");
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
    vector<TH1F*> mdp(NParams());
    vector<TH1F*> bdp(NParams());
    vector<TH1F*> fpull(NParams());
    vector<TH1F*> mpull(NParams());
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
//    double pscale = nsig/sqrt(nhits);
    double pscale =  nsig;
    for(size_t ipar=0;ipar< NParams(); ipar++){
      auto tpar = static_cast<typename KTRAJ::ParamIndex>(ipar);
      hname = string("fd") + KTRAJ::paramName(tpar);
      htitle = string("Front #Delta ") + KTRAJ::paramTitle(tpar);
      fdp[ipar] = new TH1F(hname.c_str(),htitle.c_str(),100,-pscale*sigmas[ipar],pscale*sigmas[ipar]);
      hname = string("md") + KTRAJ::paramName(tpar);
      htitle = string("Mid #Delta ") + KTRAJ::paramTitle(tpar);
      mdp[ipar] = new TH1F(hname.c_str(),htitle.c_str(),100,-pscale*sigmas[ipar],pscale*sigmas[ipar]);
      hname = string("bd") + KTRAJ::paramName(tpar);
      htitle = string("Back #Delta ") + KTRAJ::paramTitle(tpar);
      bdp[ipar] = new TH1F(hname.c_str(),htitle.c_str(),100,-pscale*sigmas[ipar],pscale*sigmas[ipar]);
      hname = string("fp") + KTRAJ::paramName(tpar);
      htitle = string("Front Pull ") + KTRAJ::paramTitle(tpar);
      fpull[ipar] = new TH1F(hname.c_str(),htitle.c_str(),100,-nsig,nsig);
      hname = string("mp") + KTRAJ::paramName(tpar);
      htitle = string("Mid Pull ") + KTRAJ::paramTitle(tpar);
      mpull[ipar] = new TH1F(hname.c_str(),htitle.c_str(),100,-nsig,nsig);
      hname = string("bp") + KTRAJ::paramName(tpar);
      htitle = string("Back Pull ") + KTRAJ::paramTitle(tpar);
      bpull[ipar] = new TH1F(hname.c_str(),htitle.c_str(),100,-nsig,nsig);
      hname = string("e") + KTRAJ::paramName(tpar);
      htitle = string("Error ") + KTRAJ::paramTitle(tpar);
      fiterrh[ipar] = new TH1F(hname.c_str(),htitle.c_str(),100,0.0,pscale*sigmas[ipar]);
      xax->SetBinLabel(ipar+1,KTRAJ::paramName(tpar).c_str());
      yax->SetBinLabel(ipar+1,KTRAJ::paramName(tpar).c_str());
    }
    TH1F* fmomres = new TH1F("fmomres","Front Momentum Resolution;P_{reco}-P_{true}(MeV/c)",100,-pscale*sigmas[NParams()],pscale*sigmas[NParams()]);
    TH1F* mmomres = new TH1F("mmomres","Mid Momentum Resolution;P_{reco}-P_{true}(MeV/c)",100,-pscale*sigmas[NParams()],pscale*sigmas[NParams()]);
    TH1F* bmomres = new TH1F("bmomres","Back Momentum Resolution;P_{reco}-P_{true}(MeV/c)",100,-pscale*sigmas[NParams()],pscale*sigmas[NParams()]);
    TH1F* fmompull = new TH1F("fmompull","Front Momentum Pull;#Delta P/#sigma _{p}",100,-nsig,nsig);
    TH1F* mmompull = new TH1F("mmompull","Mid Momentum Pull;#Delta P/#sigma _{p}",100,-nsig,nsig);
    TH1F* bmompull = new TH1F("bmompull","Back Momentum Pull;#Delta P/#sigma _{p}",100,-nsig,nsig);
    double duration (0.0);
    unsigned nfail(0), ndiv(0);

    configptr->plevel_ = Config::none;
    for(unsigned ievent=0;ievent<nevents;ievent++){
      if( (ievent % iprint) == 0) cout << "event " << ievent << endl;
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
      toy.createSeed(seedtraj,10*sigmas[NParams()],seedposvar,seedsmear);
  // if requested, constrain a parameter
      if(conspar >= 0 && conspar < (int)NParams()){
	auto const& front = tptraj.front();
	Parameters cparams = front.params();
	cparams.covariance() = seedtraj.params().covariance();
	thits.push_back(std::make_shared<PARHIT>(front.range().mid(),cparams,mask));
      }
      auto start = Clock::now();
      KKTRK kktrk(configptr,seedtraj,thits,dxings);
      auto stop = Clock::now();
      duration += std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start).count();
      auto const& fstat = kktrk.fitStatus();
      if(fstat.status_ == Status::failed)nfail++;
      if(fstat.status_ == Status::diverged)ndiv++;
      niter_ = 0;
      for(auto const& fstat: kktrk.history()){
	if(fstat.status_ != Status::unfit)niter_++;
      }
      // reset some fit parameters, to signal failed filts
      chiprob_ = -1.0;
      maxgap_ = avgap_ = -1;
      igap_ = -1;
      // fill effect information
      nkkbf_ = 0; nkkhit_ = 0; nkkmat_ = 0;
      // accumulate chisquared info
      chisq_ = fstat.chisq_.chisq();
      ndof_ = fstat.chisq_.nDOF();
      niter_ = fstat.iter_;
      nmeta_ = fstat.miter_;
      status_ = fstat.status_;
      chiprob_ = fstat.chisq_.probability();
      hinfovec.clear();
      bfinfovec.clear();
      minfovec.clear();
      tinfovec.clear();
      statush->Fill(fstat.status_);
      if(fstat.usable()){
	// truth parameters, front and back
	double ttlow = tptraj.range().begin();
	double ttmid = tptraj.range().mid();
	double tthigh = tptraj.range().end();
	KTRAJ const& fttraj = tptraj.nearestPiece(ttlow);
	KTRAJ const& mttraj = tptraj.nearestPiece(ttmid);
	KTRAJ const& bttraj = tptraj.nearestPiece(tthigh);
	for(size_t ipar=0;ipar<6;ipar++){
	  spars_.pars_[ipar] = seedtraj.params().parameters()[ipar];
	  ftpars_.pars_[ipar] = fttraj.params().parameters()[ipar];
	  mtpars_.pars_[ipar] = mttraj.params().parameters()[ipar];
	  btpars_.pars_[ipar] = bttraj.params().parameters()[ipar];
	}
	ftmom_ = tptraj.momentum(ttlow);
	mtmom_ = tptraj.momentum(ttmid);
	btmom_ = tptraj.momentum(tthigh);
	ndof->Fill(ndof_);
	chisq->Fill(chisq_);
	chisqndof->Fill(fstat.chisq_.chisqPerNDOF());
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
	    hinfo.active_ = kkhit->active();
	    hinfo.time_ = kkhit->time();
	    hinfo.chisq_ = kkhit->chisq().chisq();
	    hinfo.ndof_ = kkhit->chisq().nDOF();
	    hinfo.ambig_ = -1000;
	    hinfo.dim_ = -1000;
	    const STRAWHIT* strawhit = dynamic_cast<const STRAWHIT*>(kkhit->hit().get());
	    const SCINTHIT* scinthit = dynamic_cast<const SCINTHIT*>(kkhit->hit().get());
	    const PARHIT* constraint = dynamic_cast<const PARHIT*>(kkhit->hit().get());
	    if(strawhit != 0){
	      hinfo.ambig_ = strawhit->hitState().lrambig_;
	      hinfo.dim_ = strawhit->hitState().dimension_;
	      // straw hits can have multiple residuals
	      hinfo.type_ = HitInfo::strawtime;
	      hinfo.resid_ = strawhit->residual(WireHitState::time).value();
	      hinfo.residvar_ = strawhit->residual(WireHitState::time).variance();
	      hinfovec.push_back(hinfo);
	      //
	      hinfo.type_ = HitInfo::strawdistance;
	      hinfo.resid_ = strawhit->residual(WireHitState::distance).value();
	      hinfo.residvar_ = strawhit->residual(WireHitState::distance).variance();
	      hinfovec.push_back(hinfo);
	    } else if(scinthit != 0){
	      hinfo.type_ = HitInfo::scint;
	      hinfo.resid_ = scinthit->residual().value();
	      hinfo.residvar_ = scinthit->residual().variance();
	      hinfovec.push_back(hinfo);
	    } else if(constraint != 0){
	      hinfo.type_ = HitInfo::constraint;
	      hinfo.resid_ = sqrt(constraint->chisq().chisq());
	      hinfo.residvar_ = 1.0;
	      hinfovec.push_back(hinfo);
	    } else {
	      hinfo.type_ = HitInfo::unknown;
	      hinfo.resid_ =  0.0;
	      hinfo.residvar_ = 1.0;
	    }
	  }
	  if(kkmat != 0){
	    nkkmat_++;
	    KinKal::MaterialInfo minfo;
	    minfo.time_ = kkmat->time();
	    minfo.active_ = kkmat->active();
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
	    bfinfo.active_ = kkbf->active();
	    bfinfo.time_ = kkbf->time();
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
	double ftlow,ftmid,fthigh;
	dTraj(tptraj,fptraj,ttlow,ftlow);
	dTraj(tptraj,fptraj,ttmid,ftmid);
	dTraj(tptraj,fptraj,tthigh,fthigh);
	KTRAJ fftraj(fptraj.stateEstimate(ftlow),fptraj.charge(),tptraj.bnom(ttlow),fptraj.nearestPiece(ftlow).range());
	KTRAJ mftraj(fptraj.stateEstimate(ftmid),fptraj.charge(),tptraj.bnom(ttmid),fptraj.nearestPiece(ftmid).range());
	KTRAJ bftraj(fptraj.stateEstimate(fthigh),fptraj.charge(),tptraj.bnom(tthigh),fptraj.nearestPiece(fthigh).range());
	// fit parameters
	auto const& ffpars = fftraj.params();
	auto const& mfpars = mftraj.params();
	auto const& bfpars = bftraj.params();
	double maxgap, avgap;
	size_t igap;
	fptraj.gaps(maxgap, igap, avgap);
	maxgap_ = maxgap;
	avgap_ = avgap;
	igap_ = igap;

	Parameters ftpars, mtpars, btpars;
	ftpars = fttraj.params();
	mtpars = mttraj.params();
	btpars = bttraj.params();

	// accumulate parameter difference and pull
	vector<double> fcerr(6,0.0), mcerr(6,0.0), bcerr(6,0.0);
	for(size_t ipar=0;ipar< NParams(); ipar++){
	  fcerr[ipar] = sqrt(ffpars.covariance()[ipar][ipar]);
	  mcerr[ipar] = sqrt(mfpars.covariance()[ipar][ipar]);
	  bcerr[ipar] = sqrt(bfpars.covariance()[ipar][ipar]);
	  fdp[ipar]->Fill(ffpars.parameters()[ipar]-ftpars.parameters()[ipar]);
	  mdp[ipar]->Fill(mfpars.parameters()[ipar]-mtpars.parameters()[ipar]);
	  bdp[ipar]->Fill(bfpars.parameters()[ipar]-btpars.parameters()[ipar]);
	  fpull[ipar]->Fill((ffpars.parameters()[ipar]-ftpars.parameters()[ipar])/fcerr[ipar]);
	  mpull[ipar]->Fill((mfpars.parameters()[ipar]-mtpars.parameters()[ipar])/mcerr[ipar]);
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
	  mfitpars_.pars_[ipar] = mftraj.params().parameters()[ipar];
	  bfitpars_.pars_[ipar] = bftraj.params().parameters()[ipar];
	  ffiterrs_.pars_[ipar] = sqrt(fftraj.params().covariance()(ipar,ipar));
	  mfiterrs_.pars_[ipar] = sqrt(mftraj.params().covariance()(ipar,ipar));
	  bfiterrs_.pars_[ipar] = sqrt(bftraj.params().covariance()(ipar,ipar));
	}
	ffmom_ = fptraj.momentum(ftlow);
	mfmom_ = fptraj.momentum(ftmid);
	bfmom_ = fptraj.momentum(fthigh);
	ffmomerr_ = sqrt(fptraj.momentumVar(ftlow));
	mfmomerr_ = sqrt(fptraj.momentumVar(ftmid));
	bfmomerr_ = sqrt(fptraj.momentumVar(fthigh));
	fft_ = fptraj.range().begin();
	mft_ = fptraj.range().mid();
	bft_ = fptraj.range().end();
	fmomres->Fill((ffmom_-ftmom_));
	mmomres->Fill((mfmom_-mtmom_));
	bmomres->Fill((bfmom_-btmom_));
	fmompull->Fill((ffmom_-ftmom_)/ffmomerr_);
	mmompull->Fill((mfmom_-mtmom_)/mfmomerr_);
	bmompull->Fill((bfmom_-btmom_)/bfmomerr_);
	// state space parameter difference and errors
	//	ParticleStateEstimate tslow = tptraj.state(tlow);
	//	ParticleStateEstimate tshigh = tptraj.state(thigh);
	//	ParticleStateEstimate slow = fptraj.stateEstimate(tlow);
	//	ParticleStateEstimate shigh = fptraj.stateEstimate(thigh);

	// test
      } else if(printbad){
	cout << "Bad Fit event " << ievent << " status " << kktrk.fitStatus() << endl;
	cout << "True Traj " << tptraj << endl;
	cout << "Seed Traj " << seedtraj << endl;
	kktrk.print(cout,detail);
      }
      if(ttree)ftree->Fill();
    }
// Test fit success
    cout << nfail << " Failed fits and " << ndiv << " Diverged fits " << endl;
    hnfail->Fill(nfail);
    hndiv->Fill(ndiv);
    if(float(nfail+ndiv)/float(nevents)> 0.1){
      retval = -2;
    }
    cout <<"Time/fit = " << duration/double(nevents) << " Nanoseconds " << endl;
    // fill canvases
    TCanvas* fdpcan = new TCanvas("fdpcan","fdpcan",800,600);
    fdpcan->Divide(3,3);
    for(size_t ipar=0;ipar<NParams();++ipar){
      fdpcan->cd(ipar+1);
      fdp[ipar]->Fit("gaus","q");
    }
    fdpcan->cd(NParams()+1);
    // Test momentum resolution
    TFitResultPtr ffitr = fmomres->Fit("gaus","qS");
    if(fabs(ffitr->Parameter(1))/ffitr->Error(1) > 10.0 || ffitr->Parameter(2) > 2.0*sigmas[NParams()] ){
      cout << "Front momentum resolution out of tolerance "
      << ffitr->Parameter(1) << " +- " << ffitr->Error(1) << " sigma " << ffitr->Parameter(2) << endl;
      retval=-3;	
    }

    fdpcan->Write();
    TCanvas* mdpcan = new TCanvas("mdpcan","mdpcan",800,600);
    mdpcan->Divide(3,3);
    for(size_t ipar=0;ipar<NParams();++ipar){
      mdpcan->cd(ipar+1);
      mdp[ipar]->Fit("gaus","q");
    }
    mdpcan->cd(NParams()+1);
    // Test momentum resolution
    TFitResultPtr mfitr = mmomres->Fit("gaus","qS");
    if(fabs(mfitr->Parameter(1))/mfitr->Error(1) > 10.0 || mfitr->Parameter(2) > 2.0*sigmas[NParams()] ){
      cout << "Mid momentum resolution out of tolerance "
      << mfitr->Parameter(1) << " +- " << mfitr->Error(1) << " sigma " << mfitr->Parameter(2) << endl;
      retval=-3;	
    }
    mdpcan->Write();
    TCanvas* bdpcan = new TCanvas("bdpcan","bdpcan",800,600);
    bdpcan->Divide(3,3);
    for(size_t ipar=0;ipar<NParams();++ipar){
      bdpcan->cd(ipar+1);
      bdp[ipar]->Fit("gaus","q");
    }
    bdpcan->cd(NParams()+1);
    TFitResultPtr bfitr = bmomres->Fit("gaus","qS");
    if(fabs(bfitr->Parameter(1))/bfitr->Error(1) > 10.0 || bfitr->Parameter(2) > 2.0*sigmas[NParams()] ){
      cout << "Back momentum resolution out of tolerance "
      << bfitr->Parameter(1) << " +- " << bfitr->Error(1) << " sigma " << bfitr->Parameter(2) << endl;
      retval=-3;	
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
    TCanvas* mpullcan = new TCanvas("mpullcan","mpullcan",800,600);
    mpullcan->Divide(3,3);
    for(size_t ipar=0;ipar<NParams();++ipar){
      mpullcan->cd(ipar+1);
      mpull[ipar]->Fit("gaus","q");
    }
    mpullcan->cd(NParams()+1);
    fmompull->Fit("gaus","q");
    mpullcan->Write();
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
  cout << "Exiting with status " << retval << endl;
  exit(retval);
}
