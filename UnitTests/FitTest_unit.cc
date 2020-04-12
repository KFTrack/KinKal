//
// ToyMC test of fitting an KTRAJ-based KKTrk
//
#include "MatEnv/MatDBInfo.hh"
#include "MatEnv/DetMaterial.hh"
#include "KinKal/PKTraj.hh"
#include "KinKal/LHelix.hh"
#include "KinKal/TLine.hh"
#include "KinKal/TPoca.hh"
#include "KinKal/StrawHit.hh"
#include "KinKal/StrawMat.hh"
#include "KinKal/LightHit.hh"
#include "KinKal/BField.hh"
#include "KinKal/Vectors.hh"
#include "KinKal/KKConfig.hh"
#include "KinKal/KKHit.hh"
#include "KinKal/DXing.hh"
#include "KinKal/KKTrk.hh"
#include "CLHEP/Units/PhysicalConstants.h"

#include <iostream>
#include <getopt.h>
#include <typeinfo>
#include <vector>
#include <cmath>
#include <cfenv>
#include <memory>

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

using namespace MatEnv;
using namespace KinKal;
using namespace std;
// avoid confusion with root
using KinKal::TLine;
// define the typedefs: to change to a different trajectory implementation, just change this line
typedef LHelix KTRAJ;
typedef PKTraj<KTRAJ> PKTRAJ;
typedef KKTrk<KTRAJ> KKTRK;
typedef shared_ptr<KKConfig> KKCONFIGPTR;
typedef THit<KTRAJ> THIT;
typedef std::shared_ptr<THIT> THITPTR;
typedef DXing<KTRAJ> DXING;
typedef std::shared_ptr<DXING> DXINGPTR;
typedef StrawHit<KTRAJ> STRAWHIT;
typedef std::shared_ptr<STRAWHIT> STRAWHITPTR;
typedef LightHit<KTRAJ> LIGHTHIT;
typedef std::shared_ptr<LIGHTHIT> LIGHTHITPTR;
typedef StrawXing<KTRAJ> STRAWXING;
typedef shared_ptr<STRAWXING> STRAWXINGPTR;
typedef vector<THITPTR> THITCOL;
typedef vector<DXINGPTR> DXINGCOL;
typedef Residual<KTRAJ> RESIDUAL;
typedef TPoca<PKTRAJ,TLine> TPOCA;

// ugly global variables
float zrange(3000.0), rmax(800.0); // tracker dimension
float sprop(0.8*CLHEP::c_light), sdrift(0.065), rstraw(2.5);
float ambigdoca(-1.0);// minimum doca to set ambiguity, default sets for all hits
float sigt(3); // drift time resolution in ns
float ineff(0.1); // hit inefficiency
float tbuff(0.1);
int iseed(124223);
unsigned nhits(40);
float escale(5.0);
vector<float> sigmas = { 3.0, 3.0, 3.0, 3.0, 0.1, 3.0}; // base sigmas for parameters (per hit!)
bool simmat(true), lighthit(true), updatehits(false);
  // time hit parameters
float ttvar(0.25), twvar(100.0), shmax(80.0), vlight(0.8*CLHEP::c_light), clen(200.0);
// define the BF
Vec3 bnom(0.0,0.0,1.0);
double Bgrad(0.0), By(0.0);
BField* BF(0);
TRandom* TR = new TRandom3(iseed);
CVD2T d2t(sdrift,sigt*sigt,rstraw);

void print_usage() {
  printf("Usage: FitTest  --momentum f --costheta f --azimuth f --particle i --charge i --zrange f --nhits i --hres f --seed i --escale f --nmeta i --maxniter i --maxtemp f--ambigdoca f --ntries i --convdchisq f --simmat i --ttree i --By f --Bgrad f --TFile c --ineff f --PrintBad i --PrintDetail i --LightHit i --UpdateHits i\n");
}

struct KTRAJPars{
  Float_t pars_[KTRAJ::NParams()];
  static std::string leafnames() {
    std::string names;
    for(size_t ipar=0;ipar<KTRAJ::NParams();ipar++){
      names +=KTRAJ::paramName(static_cast<KTRAJ::ParamIndex>(ipar)) + string("/f");
      if(ipar < KTRAJ::NParams()-1)names += ":";
    }
    return names;
  }
};

struct KKHitInfo {
  Float_t resid_, residvar_, chiref_, chifit_;
  static std::string leafnames() { return std::string("resid/f:residvar/f:chiref/f:chifit/f"); }
};


// helper function
KinKal::TLine GenerateStraw(PKTRAJ const& traj, double htime) {
  // start with the true helix position at this time
  Vec4 hpos; hpos.SetE(htime);
  traj.position(hpos);
  Vec3 hdir; traj.direction(htime,hdir);
  // generate a random direction for the straw
  double eta = TR->Uniform(-M_PI,M_PI);
  Vec3 sdir(cos(eta),sin(eta),0.0);
  // generate a random drift perp to this and the trajectory
  double rdrift = TR->Uniform(-rstraw,rstraw);
  Vec3 drift = (sdir.Cross(hdir)).Unit();
  Vec3 dpos = hpos.Vect() + rdrift*drift;
//  cout << "Generating hit at position " << dpos << endl;
  double dprop = TR->Uniform(0.0,rmax);
  Vec3 mpos = dpos + sdir*dprop;
  Vec3 vprop = sdir*sprop;
  // measured time is after propagation and drift
  double tmeas = htime + dprop/sprop + fabs(rdrift)/sdrift;
  // smear measurement time
  tmeas = TR->Gaus(tmeas,sigt);
  // range doesn't really matter
  TRange trange(tmeas-dprop/sprop,tmeas+dprop/sprop);
  // construct the trajectory for this hit
  return TLine(mpos,vprop,tmeas,trange);
}

void createSeed(KTRAJ& seed){
  auto& seedpar = seed.params();
  seedpar.covariance() = ROOT::Math::SMatrixIdentity();
  for(unsigned ipar=0;ipar < 6; ipar++){
    double perr = sigmas[ipar]*escale/sqrt(nhits);
    seedpar.parameters()[ipar] += TR->Gaus(0.0,perr);
    seedpar.covariance()[ipar][ipar] *= perr*perr;
  }
}

double createHits(PKTRAJ& plhel,StrawMat const& smat, THITCOL& thits, DXINGCOL& dxings) {
  //  cout << "Creating " << nhits << " hits " << endl;
  // divide time range
  double dt = (plhel.range().range()-2*tbuff)/(nhits-1);
  double desum(0.0);
  double dscatsum(0.0);
  for(size_t ihit=0; ihit<nhits; ihit++){
    double htime = tbuff + plhel.range().low() + ihit*dt;
    auto tline = GenerateStraw(plhel,htime);
    TPoca<PKTRAJ,TLine> tp(plhel,tline);
    LRAmbig ambig(LRAmbig::null);
    if(fabs(tp.doca())> ambigdoca) ambig = tp.doca() < 0 ? LRAmbig::left : LRAmbig::right;
    // construct the hit from this trajectory
    auto sxing = std::make_shared<STRAWXING>(tp,smat);
    if(TR->Uniform(0.0,1.0) > ineff){
      thits.push_back(std::make_shared<STRAWHIT>(*BF, tline, d2t,sxing,ambig));
    } else {
      dxings.push_back(sxing);
    }
    // compute material effects and change trajectory accordingly
    if(simmat){
      auto const& endpiece = plhel.nearestPiece(tp.particleToca());
      double mom = endpiece.momentum(tp.particleToca());
      Mom4 endmom;
      endpiece.momentum(tp.particleToca(),endmom);
      Vec4 endpos; endpos.SetE(tp.particleToca());
      endpiece.position(endpos);
      std::array<float,3> dmom = {0.0,0.0,0.0}, momvar {0.0,0.0,0.0};
      sxing->momEffects(plhel,TDir::forwards, dmom, momvar);
      for(int idir=0;idir<=KInter::theta2; idir++) {
	auto mdir = static_cast<KInter::MDir>(idir);
	double momsig = sqrt(momvar[idir]);
	double dm;
	// generate a random effect given this variance and mean.  Note momEffect is scaled to momentum
	switch( mdir ) {
	  case KinKal::KInter::theta1: case KinKal::KInter::theta2 :
	    dm = TR->Gaus(dmom[idir],momsig);
	    dscatsum += dm;
	    break;
	  case KinKal::KInter::momdir :
	    dm = std::min(0.0,TR->Gaus(dmom[idir],momsig));
	    desum += dm*mom;
	    break;
	  default:
	    throw std::invalid_argument("Invalid direction");
	}
//	cout << "mom change dir " << KInter::directionName(mdir) << " mean " << dmom[idir]  << " +- " << momsig << " value " << dm  << endl;

	Vec3 dmvec;
	endpiece.dirVector(mdir,tp.particleToca(),dmvec);
	dmvec *= dm*mom;
	endmom.SetCoordinates(endmom.Px()+dmvec.X(), endmom.Py()+dmvec.Y(), endmom.Pz()+dmvec.Z(),endmom.M());
      }
	// terminate if there is catastrophic energy loss
      if(fabs(desum)/mom > 0.1)break;
      // generate a new piece and append
      KTRAJ newend(endpos,endmom,endpiece.charge(),bnom,TRange(tp.particleToca(),plhel.range().high()));
      plhel.append(newend);
    }
  }
  if(lighthit && TR->Uniform(0.0,1.0) > ineff){
    // create a LightHit at the end, axis parallel to z
    // first, find the position at showermax.
    Vec3 shmpos, hend, lmeas;
    float cstart = plhel.range().high() + 0.5;
    plhel.position(cstart,hend);
    float ltime = cstart + shmax/plhel.speed(cstart);
    plhel.position(ltime,shmpos); // true position at shower-max
    // smear the x-y position by the transverse variance.
    lmeas.SetX(TR->Gaus(shmpos.X(),sqrt(twvar)));
    lmeas.SetY(TR->Gaus(shmpos.Y(),sqrt(twvar)));
    // set the z position to the sensor plane (end of the crystal)
    lmeas.SetZ(hend.Z()+clen);
    // set the measurement time to correspond to the light propagation from showermax, smeared by the resolution
    float tmeas = TR->Gaus(ltime+(lmeas.Z()-shmpos.Z())/vlight,sqrt(ttvar));
    // create the ttraj for the light propagation
    Vec3 lvel(0.0,0.0,vlight);
    TLine lline(lmeas,lvel,tmeas);
    // then create the hit and add it; the hit has no material
    thits.push_back(std::make_shared<LIGHTHIT>(lline, ttvar, twvar));
 // test
//    cout << "cstart " << cstart << " pos " << hend << endl;
//    cout << "shmax " << ltime << " pos " << shmpos  << endl;
//    Vec3 lhpos;
//    lline.position(tmeas,lhpos);
//    cout << "tmeas " <<  tmeas  << " pos " << lmeas  << " llinepos " << lhpos << endl;
//    RESIDUAL lres;
//    thits.back()->resid(plhel,lres);
//    cout << "LightHit " << lres << endl;
//    TPOCA tpl(plhel,lline);
//    cout <<"Light TPOCA ";
//    tpl.print(cout,2);
  }

//  cout << "Total energy loss " << desum << " scattering " << dscatsum << endl;
  return desum;
}

int main(int argc, char **argv) {
// enable throw on FPE
  fetestexcept(FE_ALL_EXCEPT );

  int opt;
  double mom(105.0), cost(0.7), phi(0.5);
  double masses[5]={0.511,105.66,139.57, 493.68, 938.0};
  int imass(0), icharge(-1);
  double pmass;
  unsigned maxniter(10);
  unsigned nmeta(2);
  float maxtemp(0.0);
  unsigned ntries(1000);
  double convdchisq(0.1);
  bool ttree(true), printbad(false);
  string tfname("FitTest.root");
  int detail(0);

  static struct option long_options[] = {
    {"momentum",     required_argument, 0, 'm' },
    {"costheta",     required_argument, 0, 'c'  },
    {"azimuth",     required_argument, 0, 'a'  },
    {"particle",     required_argument, 0, 'p'  },
    {"charge",     required_argument, 0, 'q'  },
    {"zrange",     required_argument, 0, 'z'  },
    {"seed",     required_argument, 0, 's'  },
    {"hres",     required_argument, 0, 'h'  },
    {"nhits",     required_argument, 0, 'n'  },
    {"ineff",     required_argument, 0, 'i'  },
    {"escale",     required_argument, 0, 'e'  },
    {"maxniter",     required_argument, 0, 'x'  },
    {"nmeta",     required_argument, 0, 'M'  },
    {"maxtemp",     required_argument, 0, 'X'  },
    {"simmat",     required_argument, 0, 'b'  },
    {"ambigdoca",     required_argument, 0, 'd'  },
    {"ntries",     required_argument, 0, 't'  },
    {"convdchisq",     required_argument, 0, 'C'  },
    {"ttree",     required_argument, 0, 'r'  },
    {"TFile",     required_argument, 0, 'T'  },
    {"By",     required_argument, 0, 'y'  },
    {"Bgrad",     required_argument, 0, 'g'  },
    {"PrintBad",     required_argument, 0, 'P'  },
    {"PrintDetail",     required_argument, 0, 'D'  },
    {"LightHit",     required_argument, 0, 'L'  },
    {"UpdateHits",     required_argument, 0, 'U'  },
  };

  int long_index =0;
  while ((opt = getopt_long_only(argc, argv,"",
	  long_options, &long_index )) != -1) {
    switch (opt) {
      case 'm' : mom = atof(optarg);
		 break;
      case 'c' : cost = atof(optarg);
		 break;
      case 'a' : phi = atof(optarg);
		 break;
      case 'p' : imass = atoi(optarg);
		 break;
      case 'q' : icharge = atoi(optarg);
		 break;
      case 'z' : zrange = atof(optarg);
		 break;
      case 'h' : sigt = atof(optarg);
		 break;
      case 'i' : ineff = atof(optarg);
		 break;
      case 'n' : nhits = atoi(optarg);
		 break;
      case 's' : iseed = atoi(optarg);
		 break;
      case 'e' : escale = atof(optarg);
		 break;
      case 'M' : nmeta = atoi(optarg);
		 break;
      case 'x' : maxniter = atoi(optarg);
		 break;
      case 'X' : maxtemp = atof(optarg);
		 break;
      case 'b' : simmat = atoi(optarg);
		 break;
      case 'L' : lighthit = atoi(optarg);
		 break;
      case 'U' : updatehits = atoi(optarg);
		 break;
      case 'r' : ttree = atoi(optarg);
		 break;
      case 'd' : ambigdoca = atof(optarg);
		 break;
      case 't' : ntries = atoi(optarg);
		 break;
      case 'C' : convdchisq= atof(optarg);
		 break;
      case 'y' : By = atof(optarg);
		 break;
      case 'g' : Bgrad = atof(optarg);
		 break;
      case 'P' : printbad = atoi(optarg);
		 break;
      case 'D' : detail = atoi(optarg);
		 break;
      case 'T' : tfname = optarg;
		 break;
      default: print_usage();
	       exit(EXIT_FAILURE);
    }
  }
  // construct BField
  if(Bgrad != 0){
    BF = new GradBField(1.0-0.5*zrange*Bgrad,1.0+0.5*zrange*Bgrad,-0.5*zrange,0.5*zrange);
    Vec3 bn;
    BF->fieldVect(bn,Vec3(0.0,0.0,0.0));
    bnom = bn;
  } else {
    bnom = Vec3(0.0,By,1.0);
    BF = new UniformBField(bnom);
  }

  MatDBInfo matdbinfo;
  const DetMaterial* wallmat = matdbinfo.findDetMaterial("straw-wall");
  const DetMaterial* gasmat = matdbinfo.findDetMaterial("straw-gas");
  const DetMaterial* wiremat = matdbinfo.findDetMaterial("straw-wire");
  float rwire(0.025), wthick(0.015);
  StrawMat smat(rstraw,wthick,rwire, *wallmat, *gasmat, *wiremat);

  pmass = masses[imass];
  Vec4 origin(0.0,0.0,0.0,0.0);
  float sint = sqrt(1.0-cost*cost);
  Mom4 momv(mom*sint*cos(phi),mom*sint*sin(phi),mom*cost,pmass);
  KTRAJ lhel(origin,momv,icharge,bnom);
  cout << "True initial " << lhel << endl;
  PKTRAJ plhel(lhel);
  // truncate the range according to the Z range
  Vec3 vel; plhel.velocity(0.0,vel);
  plhel.setRange(TRange(-0.5*zrange/vel.Z()-tbuff,0.5*zrange/vel.Z()+tbuff));
  // generate hits
  THITCOL thits; // this program shares hit ownership with KKTrk
  DXINGCOL dxings; // this program shares det xing ownership with KKTrk
  createHits(plhel,smat, thits, dxings);
//  cout << "vector of hit points " << thits.size() << endl;
//  cout << "True " << plhel << endl;
  double startmom = plhel.momentum(plhel.range().low());
  double endmom = plhel.momentum(plhel.range().high());
  Vec3 end, bend;
  plhel.front().direction(plhel.range().high(),bend);
  plhel.back().direction(plhel.range().high(),end);
  double angle = ROOT::Math::VectorUtil::Angle(bend,end);
  cout << "total momentum change = " << startmom-endmom << " total angle change = " << angle << endl;
  // create the fit seed by randomizing the parameters at the middle
  auto seedhel = plhel.nearestPiece(0.0);
  createSeed(seedhel);
  cout << "Seed params " << seedhel.params().parameters() <<" covariance " << endl << seedhel.params().covariance() << endl;
  // Create the KKTrk from these hits
  KKCONFIGPTR configptr = make_shared<KKConfig>(*BF);
  configptr->dwt_ = 1.0e6;
  configptr->convdchisq_ = convdchisq;
  configptr->maxniter_ = maxniter;
  // add schedule; MC-truth based ambiguity
  MConfig mconfig;
  mconfig.processmat_ = mconfig.processbfield_ =  true;
  mconfig.updatehits_ = updatehits;
  mconfig.hitupdateparams_.push_back(make_any<WHUParams>(ambigdoca,100.0)); // 1st parameter turns off drift, 2nd says accept all hits
  float tstep = maxtemp/(std::max(nmeta,(unsigned)1));
  float temp = maxtemp;
  for(unsigned imeta = 0; imeta< nmeta+1; imeta++){
    mconfig.temp_ = temp;
    temp -= tstep;
    configptr->schedule_.push_back(mconfig);
  }
   // cout << *configptr << endl;
// create and fit the track
  KKTRK kktrk(configptr,seedhel,thits,dxings);
  kktrk.print(cout,detail);
  TFile fitfile(tfname.c_str(),"RECREATE");
  // tree variables
  KTRAJPars ftpars_, etpars_, spars_, ffitpars_, ffiterrs_, efitpars_, efiterrs_;
  float chisq_, etmom_, ftmom_, ffmom_, efmom_, tde_, chiprob_;
  float fft_,eft_;
  int ndof_, niter_, status_;
  if(ntries <=0 ){
    // draw the fit result
    TCanvas* pttcan = new TCanvas("pttcan","PieceKTRAJ",1000,1000);
    auto const& fithel = kktrk.fitTraj();
    unsigned np = fithel.range().range()*fithel.speed(fithel.range().mid());
    TPolyLine3D* fitpl = new TPolyLine3D(np);
    fitpl->SetLineColor(kBlack);
    fitpl->SetLineStyle(kSolid);
    double ts = fithel.range().range()/(np-1);
    for(unsigned ip=0;ip<np;ip++){
      double tp = fithel.range().low() + ip*ts;
      Vec3 ppos;
      fithel.position(tp,ppos);
      fitpl->SetPoint(ip,ppos.X(),ppos.Y(),ppos.Z());
    }
    fitpl->Draw();
// now draw the truth
    TPolyLine3D* thelpl = new TPolyLine3D(np);
    thelpl->SetLineColor(kGreen);
    thelpl->SetLineStyle(kDotted);
    ts = plhel.range().range()/(np-1);
    for(unsigned ip=0;ip<np;ip++){
      double tp = plhel.range().low() + ip*ts;
      Vec3 ppos;
      plhel.position(tp,ppos);
      thelpl->SetPoint(ip,ppos.X(),ppos.Y(),ppos.Z());
    }
    thelpl->Draw();

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
    vector<KKHitInfo> hinfo;
    if(ttree){
      ftree = new TTree("fit","fit");
      ftree->Branch("ftpars.", &ftpars_,KTRAJPars::leafnames().c_str());
      ftree->Branch("etpars.", &etpars_,KTRAJPars::leafnames().c_str());
      ftree->Branch("spars.", &spars_,KTRAJPars::leafnames().c_str());
      ftree->Branch("ffpars.", &ffitpars_,KTRAJPars::leafnames().c_str());
      ftree->Branch("fferrs.", &ffiterrs_,KTRAJPars::leafnames().c_str());
      ftree->Branch("efpars.", &efitpars_,KTRAJPars::leafnames().c_str());
      ftree->Branch("eferrs.", &efiterrs_,KTRAJPars::leafnames().c_str());
      ftree->Branch("chisq", &chisq_,"chisq/F");
      ftree->Branch("ndof", &ndof_,"ndof/I");
      ftree->Branch("chiprob", &chiprob_,"chiprob/F");
      ftree->Branch("niter", &niter_,"niter/I");
      ftree->Branch("status", &status_,"status/I");
      ftree->Branch("ftmom", &ftmom_,"ftmom/F");
      ftree->Branch("etmom", &etmom_,"etmom/F");
      ftree->Branch("ffmom", &ffmom_,"ffmom/F");
      ftree->Branch("efmom", &efmom_,"efmom/F");
      ftree->Branch("fft", &fft_,"fft/F");
      ftree->Branch("eft", &eft_,"eft/F");
      ftree->Branch("tde", &tde_,"tde/F");
      ftree->Branch("hinfo",&hinfo);
    }
    // now repeat this to gain statistics
    vector<TH1F*> fdpgenh(KTRAJ::NParams());
    vector<TH1F*> bdpgenh(KTRAJ::NParams());
    vector<TH1F*> fdpullgenh(KTRAJ::NParams());
    vector<TH1F*> bdpullgenh(KTRAJ::NParams());
    vector<TH1F*> fiterrh(KTRAJ::NParams());
    TH1F* niter = new TH1F("niter", "N Iterations", 100,0,100);
    TH1F* ndof = new TH1F("ndof", "N Degree of Freedom", 100,0,100);
    TH1F* chisq = new TH1F("chisq", "Chisquared", 100,0,100);
    TH1F* chisqndof = new TH1F("chisqndof", "Chisquared per NDOF", 100,0,10.0);
    TH1F* chisqprob = new TH1F("chisqprob", "Chisquared probability", 100,0,1.0);
    TH1F* logchisqprob = new TH1F("logchisqprob", "Chisquared probability", 100,-10,0.0);
    string htitle, hname;
    TH2F* corravg = new TH2F("corravg","Average correlation matrix magnitudes",KTRAJ::NParams(),-0.5,KTRAJ::NParams()-0.5,KTRAJ::NParams(), -0.5,KTRAJ::NParams()-0.5);
    TAxis* xax = corravg->GetXaxis();
    TAxis* yax = corravg->GetYaxis();
    double nsig(10.0);
    double pscale = nsig/sqrt(nhits);
    for(size_t ipar=0;ipar< KTRAJ::NParams(); ipar++){
      auto tpar = static_cast<KTRAJ::ParamIndex>(ipar);
      hname = string("fd") + KTRAJ::paramName(tpar);
      htitle = string("Front #Delta ") + KTRAJ::paramTitle(tpar);
      fdpgenh[ipar] = new TH1F(hname.c_str(),htitle.c_str(),100,-pscale*sigmas[ipar],pscale*sigmas[ipar]);
      hname = string("bd") + KTRAJ::paramName(tpar);
      htitle = string("Back #Delta ") + KTRAJ::paramTitle(tpar);
      bdpgenh[ipar] = new TH1F(hname.c_str(),htitle.c_str(),100,-pscale*sigmas[ipar],pscale*sigmas[ipar]);
      hname = string("fp") + KTRAJ::paramName(tpar);
      htitle = string("Front Pull ") + KTRAJ::paramTitle(tpar);
      fdpullgenh[ipar] = new TH1F(hname.c_str(),htitle.c_str(),100,-nsig,nsig);
      hname = string("bp") + KTRAJ::paramName(tpar);
      htitle = string("Back Pull ") + KTRAJ::paramTitle(tpar);
      bdpullgenh[ipar] = new TH1F(hname.c_str(),htitle.c_str(),100,-nsig,nsig);
      hname = string("e") + KTRAJ::paramName(tpar);
      htitle = string("Error ") + KTRAJ::paramTitle(tpar);
      fiterrh[ipar] = new TH1F(hname.c_str(),htitle.c_str(),100,0.0,pscale*sigmas[ipar]);
      xax->SetBinLabel(ipar+1,KTRAJ::paramName(tpar).c_str());
      yax->SetBinLabel(ipar+1,KTRAJ::paramName(tpar).c_str());
    }
    for(unsigned itry=0;itry<ntries;itry++){
      // randomize the helix
      Vec4 torigin(TR->Gaus(0.0,3.0), TR->Gaus(0.0,3.0), TR->Gaus(0.0,3.0),TR->Gaus(0.0,3.0));
      double tphi = TR->Uniform(-M_PI,M_PI);
      double tcost = TR->Uniform(0.5,0.8);
      double tsint = sqrt(1.0-tcost*tcost);
      Mom4 tmomv(mom*tsint*cos(tphi),mom*tsint*sin(tphi),mom*tcost,pmass);
      KTRAJ tlhel(torigin,tmomv,icharge,bnom);
      PKTRAJ tplhel(tlhel);
      Vec3 vel; tplhel.velocity(0.0,vel);
      tplhel.setRange(TRange(-0.5*zrange/vel.Z()-tbuff,0.5*zrange/vel.Z()+tbuff));
      thits.clear();
      dxings.clear();
      tde_ = createHits(tplhel,smat, thits,dxings);
      auto seedhel = tplhel.nearestPiece(0.0);
      createSeed(seedhel);
      KKTRK kktrk(configptr,seedhel,thits,dxings);
      // compare parameters at the first traj of both true and fit
      auto const& tpars = tplhel.front().params();
      auto const& fpars = kktrk.fitTraj().front().params();
      auto const& btpars = tplhel.back().params();
      auto const& bfpars = kktrk.fitTraj().back().params();
     // momentum
      // accumulate parameter difference and pull
      vector<double> cerr(6,0.0), bcerr(6,0.0);
      for(size_t ipar=0;ipar< KTRAJ::NParams(); ipar++){
	cerr[ipar] = sqrt(fpars.covariance()[ipar][ipar]);
	bcerr[ipar] = sqrt(bfpars.covariance()[ipar][ipar]);
	fdpgenh[ipar]->Fill(fpars.parameters()[ipar]-tpars.parameters()[ipar]);
	bdpgenh[ipar]->Fill(bfpars.parameters()[ipar]-btpars.parameters()[ipar]);
	fdpullgenh[ipar]->Fill((fpars.parameters()[ipar]-tpars.parameters()[ipar])/cerr[ipar]);
	bdpullgenh[ipar]->Fill((bfpars.parameters()[ipar]-btpars.parameters()[ipar])/bcerr[ipar]);
	fiterrh[ipar]->Fill(cerr[ipar]);
      }
      // accumulate average correlation matrix
      auto const& cov = fpars.covariance();
      //    auto cormat = cov;
      for(unsigned ipar=0; ipar <KTRAJ::NParams();ipar++){
	for(unsigned jpar=ipar;jpar < KTRAJ::NParams(); jpar++){
	  double corr = cov[ipar][jpar]/(cerr[ipar]*cerr[jpar]);
	  //	cormat[ipar][jpar] = corr;
	  corravg->Fill(ipar,jpar,fabs(corr));
	}
      }
      // accumulate chisquared info
      auto const& fstat = kktrk.fitStatus();
      chiprob_ = fstat.prob_; 
      niter->Fill(fstat.iter_);
      ndof->Fill(fstat.ndof_);
      chisq->Fill(fstat.chisq_);
      chisqndof->Fill(fstat.chisq_/fstat.ndof_);
      chisqprob->Fill(chiprob_);
      logchisqprob->Fill(log10(chiprob_));
      // fill tree
      for(size_t ipar=0;ipar<6;ipar++){
	spars_.pars_[ipar] = seedhel.params().parameters()[ipar];
	ftpars_.pars_[ipar] = tplhel.front().params().parameters()[ipar];
	etpars_.pars_[ipar] = tplhel.back().params().parameters()[ipar];
	ffitpars_.pars_[ipar] = kktrk.fitTraj().front().params().parameters()[ipar];
	efitpars_.pars_[ipar] = kktrk.fitTraj().back().params().parameters()[ipar];
	ffiterrs_.pars_[ipar] = sqrt(kktrk.fitTraj().front().params().covariance()(ipar,ipar));
	efiterrs_.pars_[ipar] = sqrt(kktrk.fitTraj().back().params().covariance()(ipar,ipar));
      }
      ftmom_ = tplhel.front().momentum(tplhel.range().low());
      etmom_ = tplhel.back().momentum(tplhel.range().high());
      ffmom_ = kktrk.fitTraj().front().momentum(kktrk.fitTraj().range().low());
      efmom_ = kktrk.fitTraj().back().momentum(kktrk.fitTraj().range().high());
      fft_ = kktrk.fitTraj().range().low();
      eft_ = kktrk.fitTraj().range().high();
      chisq_ = fstat.chisq_;
      ndof_ = fstat.ndof_;
      niter_ = fstat.iter_;
      status_ = fstat.status_;

      // test
      if(printbad && !kktrk.fitStatus().usable()){
	cout << "Bad Fit try " << itry << endl;
	cout << "True Traj " << tlhel << endl;
	cout << "Seed Traj " << seedhel << endl;
	kktrk.print(cout,detail);
      }
      if(ttree)ftree->Fill();
    }
    // fill canvases
    TCanvas* fdpcan = new TCanvas("fdpcan","fdpcan",800,600);
    fdpcan->Divide(3,2);
    for(size_t ipar=0;ipar<KTRAJ::NParams();++ipar){
      fdpcan->cd(ipar+1);
      fdpgenh[ipar]->Fit("gaus","q");
    }
    fdpcan->Write();
    TCanvas* bdpcan = new TCanvas("bdpcan","bdpcan",800,600);
    bdpcan->Divide(3,2);
    for(size_t ipar=0;ipar<KTRAJ::NParams();++ipar){
      bdpcan->cd(ipar+1);
      bdpgenh[ipar]->Fit("gaus","q");
    }
    bdpcan->Write();
    TCanvas* fpullcan = new TCanvas("fpullcan","fpullcan",800,600);
    fpullcan->Divide(3,2);
    for(size_t ipar=0;ipar<KTRAJ::NParams();++ipar){
      fpullcan->cd(ipar+1);
      fdpullgenh[ipar]->Fit("gaus","q");
    }
    fpullcan->Write();
    TCanvas* bpullcan = new TCanvas("bpullcan","bpullcan",800,600);
    bpullcan->Divide(3,2);
    for(size_t ipar=0;ipar<KTRAJ::NParams();++ipar){
      bpullcan->cd(ipar+1);
      bdpullgenh[ipar]->Fit("gaus","q");
    }
    bpullcan->Write();
    TCanvas* perrcan = new TCanvas("perrcan","perrcan",800,600);
    perrcan->Divide(3,2);
    for(size_t ipar=0;ipar<KTRAJ::NParams();++ipar){
      perrcan->cd(ipar+1);
      fiterrh[ipar]->Draw();
    }
    perrcan->Write();
    TCanvas* corrcan = new TCanvas("corrcan","corrcan",600,600);
    corrcan->Divide(1,1);
    corrcan->cd(1);
    corravg->Scale(1.0/float(ntries));
    corravg->SetStats(0);
    gPad->SetLogz();
    corravg->Draw("colorztext0");
    corrcan->Write();

    TCanvas* statuscan = new TCanvas("statuscan","statuscan",800,600);
    statuscan->Divide(3,2);
    statuscan->cd(1);
    niter->Draw();
    statuscan->cd(2);
    ndof->Draw();
    statuscan->cd(3);
    chisq->Draw();
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
  exit(EXIT_SUCCESS);
}
