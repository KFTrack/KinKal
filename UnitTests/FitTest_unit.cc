//
// ToyMC test of fitting an LHelix-based KKTrk
//
#include "MatEnv/MatDBInfo.hh"
#include "MatEnv/DetMaterial.hh"
#include "KinKal/PKTraj.hh"
#include "KinKal/LHelix.hh"
#include "KinKal/TLine.hh"
#include "KinKal/TPoca.hh"
#include "KinKal/StrawHit.hh"
#include "KinKal/StrawMat.hh"
#include "KinKal/Context.hh"
#include "KinKal/Vectors.hh"
#include "KinKal/KKHit.hh"
#include "KinKal/KTMIsect.hh"
#include "KinKal/KKTrk.hh"
#include "CLHEP/Units/PhysicalConstants.h"

#include <iostream>
#include <getopt.h>
#include <typeinfo>
#include <vector>
#include <cmath>

#include "TH1F.h"
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
#include "TF1.h"
#include "TDirectory.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "Math/VectorUtil.h"

using namespace KinKal;
using namespace std;
// avoid confusion with root
using KinKal::TLine;
typedef KinKal::PKTraj<LHelix> PLHelix;
typedef KinKal::KKTrk<LHelix> KKTRK;
// ugly global variables
double zrange(3000.0), rmax(800.0); // tracker dimension
double sprop(0.8*CLHEP::c_light), sdrift(0.065), rstraw(2.5);
double ambigdoca(0.5);// minimum doca to set ambiguity
double sigt(3); // drift time resolution in ns
double tbuff(0.1);
int iseed(124223);
unsigned nhits(40);
double escale(5.0);
vector<double> sigmas = { 3.0, 3.0, 3.0, 3.0, 0.03, 3.0}; // base sigmas for parameters (per hit!)
// define the context
UniformBField BF(1.0); // 1 Tesla
Context context(BF);
TRandom* TR = new TRandom3(iseed);
CVD2T d2t(sdrift,sigt*sigt);

void print_usage() {
  printf("Usage: FitTest  --momentum f --costheta f --azimuth f --particle i --charge i --zrange f --nhits i --hres f --seed i --escale f --maxniter f --ambigdoca f --ntries i --addmat i\n");
}

// helper function
KinKal::TLine GenerateStraw(PLHelix const& traj, double htime) {
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
  // find the ends where this reaches the cylinder
  double drho = dpos.Rho();
  double ddot = dpos.Dot(sdir);
  double shlen = sqrt(rmax*rmax + ddot*ddot - drho*drho); // straw half-length
  // choose the closest end to be the measurement end
  double rprop = (fabs(-ddot-shlen) < fabs(-ddot+shlen)) ?  -ddot-shlen : -ddot+shlen;
  // sign propagation velocity away from the measurement
  if(rprop>0){
    sdir *= -1.0;
    rprop *= -1.0;
  }
  Vec3 mpos = dpos - sdir*rprop;
  Vec3 vprop = sdir*sprop;
  // measured time is after propagation and drift
  double tmeas = hpos.T() + fabs(rprop)/sprop + fabs(rdrift)/sdrift;
  // smear measurement time
  tmeas = TR->Gaus(tmeas,sigt);
  // measurement time is the longest time
  TRange trange(tmeas-2.0*shlen/sprop,tmeas);
  // construct the trajectory for this hit
  return TLine(mpos,vprop,tmeas,trange);
}

void createSeed(LHelix& seed){
  auto& seedpar = seed.params();
  seedpar.covariance() = ROOT::Math::SMatrixIdentity();
  for(unsigned ipar=0;ipar < 6; ipar++){
    double perr = sigmas[ipar]*escale/sqrt(nhits);
    seedpar.parameters()[ipar] += TR->Gaus(0.0,perr);
    seedpar.covariance()[ipar][ipar] *= perr*perr;
  }
}

void createHits(PLHelix& plhel,StrawMat const& smat, std::vector<StrawHit>& shits,bool addmat) {
  //  cout << "Creating " << nhits << " hits " << endl;
  // divide time range
  double dt = (plhel.range().range()-2*tbuff)/(nhits-1);
  double desum(0.0);
  double dscatsum(0.0);
  for(size_t ihit=0; ihit<nhits; ihit++){
    double htime = tbuff + plhel.range().low() + ihit*dt;
    auto tline = GenerateStraw(plhel,htime);
    TDPoca<PLHelix,TLine> tp(plhel,tline);
    WireHit::LRAmbig ambig(WireHit::null);
    if(fabs(tp.doca())> ambigdoca) ambig = tp.doca() < 0 ? WireHit::left : WireHit::right;
    // construct the hit from this trajectory
    StrawHit sh(tline,context,d2t,smat,ambigdoca,ambig);
    shits.push_back(sh);
    // compute material effects and change trajectory accordingly
    bool simmat(true);
    if(addmat && simmat){
      std::vector<MIsect> misects;
      smat.intersect(tp,misects);
      auto const& endpiece = plhel.nearestPiece(tp.t0());
      KTMIsect<LHelix> ktmi(endpiece,tp.t0(),misects);
      double mom = endpiece.momentum(tp.t0());
      Mom4 endmom;
      endpiece.momentum(tp.t0(),endmom);
      Vec4 endpos; endpos.SetE(tp.t0());
      endpiece.position(endpos);
      for(int idir=0;idir<=KInter::theta2; idir++) {
	auto mdir = static_cast<KInter::MDir>(idir);
	double dmom, momvar, dm;
	ktmi.momEffect(TDir::forwards, mdir, dmom, momvar);
	// generate a random effect given this variance and mean.  Note momEffect is scaled to momentum
	switch( mdir ) {
	  case KinKal::KInter::theta1: case KinKal::KInter::theta2 :
	    dm = TR->Gaus(dmom,sqrt(momvar));
	    dscatsum += dm;
	    break;
	  case KinKal::KInter::momdir :
//	    dm = -TR->Landau(-dmom,sqrt(momvar)); // have to flip sign for root!
	    dm = dmom;
	    desum += dm*mom;
	    break;
	  default:
	    throw std::invalid_argument("Invalid direction");
	}
//	cout << "direction " << mdir << " doca " << tp.doca() << " dmom " << dmom << " +- " << sqrt(momvar) << " sample " << dm  << endl;
	Vec3 dmvec;
	endpiece.dirVector(mdir,tp.t0(),dmvec);
	dmvec *= dm*mom;
//	dmvec *= 0.0;
	endmom.SetCoordinates(endmom.Px()+dmvec.X(), endmom.Py()+dmvec.Y(), endmom.Pz()+dmvec.Z(),endmom.M());
      }
      // generate a new piece and append
      LHelix newend(endpos,endmom,endpiece.charge(),context,TRange(tp.t0(),plhel.range().high()));
      if(!plhel.append(newend))
	cout << "Error appending traj " << newend << endl;
    }
  }
//  cout << "Total energy loss " << desum << " scattering " << dscatsum << endl;
}

int main(int argc, char **argv) {
  int opt;
  double mom(105.0), cost(0.7), phi(0.5);
  double masses[5]={0.511,105.66,139.57, 493.68, 938.0};
  int imass(0), icharge(-1);
  double pmass;
  unsigned maxniter(10);
  unsigned ntries(1000);
  double mindchisq(0.1);
  bool addmat(true);

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
    {"escale",     required_argument, 0, 'e'  },
    {"maxniter",     required_argument, 0, 'x'  },
    {"addmat",     required_argument, 0, 'b'  },
    {"ambigdoca",     required_argument, 0, 'd'  },
    {"ntries",     required_argument, 0, 't'  },
    {"mindchisq",     required_argument, 0, 'i'  },
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
      case 'n' : nhits = atoi(optarg);
		 break;
      case 's' : iseed = atoi(optarg);
		 break;
      case 'e' : escale = atof(optarg);
		 break;
      case 'x' : maxniter = atoi(optarg);
		 break;
      case 'b' : addmat = atoi(optarg);
		 break;
      case 'd' : ambigdoca = atof(optarg);
		 break;
      case 't' : ntries = atoi(optarg);
		 break;
      case 'i' : mindchisq= atof(optarg);
		 break;
      default: print_usage();
	       exit(EXIT_FAILURE);
    }
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
  LHelix lhel(origin,momv,icharge,context);
  cout << "True initial " << lhel << endl;
  PLHelix plhel(lhel);
  // truncate the range according to the Z range
  Vec3 vel; plhel.velocity(0.0,vel);
  plhel.setRange(TRange(-0.5*zrange/vel.Z()-tbuff,0.5*zrange/vel.Z()+tbuff));
  // generate hits
  std::vector<StrawHit> shits; // this owns the hits
  std::vector<const THit*> thits; // this references them as THits
  createHits(plhel,smat, shits,addmat);
  for(auto& shit : shits) thits.push_back(&shit);
//  cout << "vector of hit points " << thits.size() << endl;
  cout << "True " << plhel << endl;
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
  Config config;
  config.dwt_ = 1.0e6;
  config.mindchisq_ = mindchisq;
  config.maxniter_ = maxniter;
  config.addmat_ = addmat;
  KKTRK kktrk(seedhel,thits,config);
  // fit the track
  kktrk.fit();
  auto const& effs = kktrk.effects();
  cout << "KKTrk " << kktrk.status() << endl;
  for(auto const& eff : effs) {
    cout << "Eff at time " << eff->time() << " status " << eff->status(TDir::forwards)  << " " << eff->status(TDir::backwards);
    auto ihit = dynamic_cast<const KKHit<LHelix>*>(eff.get());
    auto imhit = dynamic_cast<const KKMHit<LHelix>*>(eff.get());
    if(ihit != 0){
      cout << " Hit status " << ihit->poca().status() << " doca " << ihit->poca().doca() << ihit->refResid() << endl;
    } else if(imhit != 0){
      cout << " MHit status " << imhit->hit().poca().status() << " doca " << imhit->hit().poca().doca() << imhit->hit().refResid() << endl;
    } else
      cout << endl;
  }
  TFile fitfile("FitTest.root","RECREATE");
  if(ntries <=0 ){
  // draw the fit result
    TCanvas* pttcan = new TCanvas("pttcan","PieceLHelix",1000,1000);
    auto const& fithel = kktrk.fitTraj();
    unsigned np = fithel.range().range()*fithel.speed(fithel.range().mid());
    TPolyLine3D* fitpl = new TPolyLine3D(np);
    fitpl->SetLineColor(kBlack);
    fitpl->SetLineStyle(kSolid);
    cout << "Fit Result " << fithel << endl;
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
    // now repeat this to gain statistics
    vector<TH1F*> dpgenh(LHelix::NParams());
    vector<TH1F*> dpullgenh(LHelix::NParams());
    vector<TH1F*> fiterrh(LHelix::NParams());
    TH1F* niter = new TH1F("niter", "N Iterations", 100,0,100);
    TH1F* ndof = new TH1F("ndof", "N Degree of Freedom", 100,0,100);
    TH1F* chisq = new TH1F("chisq", "Chisquared", 100,0,100);
    TH1F* chisqndof = new TH1F("chisqndof", "Chisquared per NDOF", 100,0,10.0);
    TH1F* chisqprob = new TH1F("chisqprob", "Chisquared probability", 100,0,1.0);
    TH1F* logchisqprob = new TH1F("logchisqprob", "Chisquared probability", 100,-10,0.0);
    string htitle, hname;
    TH2F* corravg = new TH2F("corravg","Average correlation matrix magnitudes",LHelix::NParams(),-0.5,LHelix::NParams()-0.5,LHelix::NParams(), -0.5,LHelix::NParams()-0.5);
    TAxis* xax = corravg->GetXaxis();
    TAxis* yax = corravg->GetYaxis();
    double nsig(10.0);
    double pscale = nsig/sqrt(nhits);
    for(size_t ipar=0;ipar< LHelix::NParams(); ipar++){
      auto tpar = static_cast<LHelix::ParamIndex>(ipar);
      hname = string("d") + LHelix::paramName(tpar);
      htitle = string("#Delta ") + LHelix::paramTitle(tpar);
      dpgenh[ipar] = new TH1F(hname.c_str(),htitle.c_str(),100,-pscale*sigmas[ipar],pscale*sigmas[ipar]);
      hname = string("p") + LHelix::paramName(tpar);
      htitle = string("Pull ") + LHelix::paramTitle(tpar);
      dpullgenh[ipar] = new TH1F(hname.c_str(),htitle.c_str(),100,-nsig,nsig);
      hname = string("e") + LHelix::paramName(tpar);
      htitle = string("Error ") + LHelix::paramTitle(tpar);
      fiterrh[ipar] = new TH1F(hname.c_str(),htitle.c_str(),100,0.0,pscale*sigmas[ipar]);
      xax->SetBinLabel(ipar+1,LHelix::paramName(tpar).c_str());
      yax->SetBinLabel(ipar+1,LHelix::paramName(tpar).c_str());
    }
    for(unsigned itry=0;itry<ntries;itry++){
      // randomize the helix
      Vec4 torigin(TR->Gaus(0.0,3.0), TR->Gaus(0.0,3.0), TR->Gaus(0.0,3.0),TR->Gaus(0.0,3.0));
      double tphi = TR->Uniform(-M_PI,M_PI);
      double tcost = TR->Uniform(0.5,0.8);
      double tsint = sqrt(1.0-tcost*tcost);
      Mom4 tmomv(mom*tsint*cos(tphi),mom*tsint*sin(tphi),mom*tcost,pmass);
      LHelix tlhel(torigin,tmomv,icharge,context);
      PLHelix tplhel(tlhel);
      Vec3 vel; tplhel.velocity(0.0,vel);
      tplhel.setRange(TRange(-0.5*zrange/vel.Z()-tbuff,0.5*zrange/vel.Z()+tbuff));
      shits.clear();
      createHits(tplhel,smat, shits,addmat);
      auto seedhel = tplhel.nearestPiece(0.0);
      createSeed(seedhel);
      thits.clear();
      for(auto& shit : shits) thits.push_back(&shit);
      KKTRK kktrk(seedhel,thits,config);
      kktrk.fit();
      // compare parameters at the first traj of both true and fit
      auto const& tpars = tplhel.front().params();
      auto const& fpars = kktrk.fitTraj().front().params();
      // accumulate parameter difference and pull
      vector<double> cerr(6,0.0);
      for(size_t ipar=0;ipar< LHelix::NParams(); ipar++){
	cerr[ipar] = sqrt(fpars.covariance()[ipar][ipar]);
	dpgenh[ipar]->Fill(fpars.parameters()[ipar]-tpars.parameters()[ipar]);
	dpullgenh[ipar]->Fill((fpars.parameters()[ipar]-tpars.parameters()[ipar])/cerr[ipar]);
	fiterrh[ipar]->Fill(cerr[ipar]);
      }
      // accumulate average correlation matrix
      auto const& cov = fpars.covariance();
      //    auto cormat = cov;
      for(unsigned ipar=0; ipar <LHelix::NParams();ipar++){
	for(unsigned jpar=ipar;jpar < LHelix::NParams(); jpar++){
	  double corr = cov[ipar][jpar]/(cerr[ipar]*cerr[jpar]);
	  //	cormat[ipar][jpar] = corr;
	  corravg->Fill(ipar,jpar,fabs(corr));
	}
      }
      // accumulate chisquared info
      double chiprob = TMath::Prob(kktrk.status().chisq_,kktrk.status().ndof_);
      niter->Fill(kktrk.status().niter_);
      ndof->Fill(kktrk.status().ndof_);
      chisq->Fill(kktrk.status().chisq_);
      chisqndof->Fill(kktrk.status().chisq_/kktrk.status().ndof_);
      chisqprob->Fill(chiprob);
      logchisqprob->Fill(log10(chiprob));
    }
    TCanvas* dpcan = new TCanvas("dpcan","dpcan",800,600);
    dpcan->Divide(3,2);
    for(size_t ipar=0;ipar<LHelix::NParams();++ipar){
      dpcan->cd(ipar+1);
      dpgenh[ipar]->Fit("gaus");
    }
    dpcan->Write();
    TCanvas* pullcan = new TCanvas("pullcan","pullcan",800,600);
    pullcan->Divide(3,2);
    for(size_t ipar=0;ipar<LHelix::NParams();++ipar){
      pullcan->cd(ipar+1);
      dpullgenh[ipar]->Fit("gaus");
    }
    pullcan->Write();
    TCanvas* perrcan = new TCanvas("perrcan","perrcan",800,600);
    perrcan->Divide(3,2);
    for(size_t ipar=0;ipar<LHelix::NParams();++ipar){
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
