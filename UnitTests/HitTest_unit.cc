//
// ToyMC test of hits on a LHelix-based
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
#include "CLHEP/Units/PhysicalConstants.h"

#include <iostream>
#include <stdio.h>
#include <iostream>
#include <getopt.h>

#include "TH1F.h"
#include "TSystem.h"
#include "THelix.h"
#include "TPolyLine3D.h"
#include "TFile.h"
#include "TAxis3D.h"
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

using namespace KinKal;
using namespace std;
// avoid confusion with root
using KinKal::TLine;
typedef KinKal::PKTraj<LHelix> PLHelix;

// ugly global variables
double zrange(3000.0), rmax(800.0); // tracker dimension
double sprop(0.8*CLHEP::c_light), sdrift(0.065), rstraw(2.5);
double sigt(3); // drift time resolution in ns

void print_usage() {
  printf("Usage: FitTest  --momentum f --costheta f --azimuth f --particle i --charge i --zrange f --nhits i --hres f --seed i --ambigdoca f --ddoca f\n");
}

// helper function
KinKal::TLine GenerateStraw(PLHelix const& traj, double htime, TRandom* TR) {
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

int main(int argc, char **argv) {
  int opt;
  double mom(105.0), cost(0.7), phi(0.5);
  double ambigdoca(0.5); // minimum DOCA to consider sign accurate
  double masses[5]={0.511,105.66,139.57, 493.68, 938.0};
  int imass(0), icharge(-1);
  double pmass;
  unsigned nhits(40);
  int iseed(124223);
  float rwire(0.025), wthick(0.015);
  float ddoca(0.1);
  float scatfrac(0.995);

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
    {"ambigdoca",     required_argument, 0, 'd'  },
    {"ddoca",     required_argument, 0, 'x'  },
  };

  int long_index =0;
  while ((opt = getopt_long_only(argc, argv,"",
	  long_options, &long_index )) != -1) {
    cout << "found option " << opt << " argument " <<optarg << endl;
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
      case 'd' : ambigdoca = atof(optarg);
		 break;
      case 'x' : ddoca = atof(optarg);
		 break;
      default: print_usage();
	       exit(EXIT_FAILURE);
    }
  }

  TRandom* TR = new TRandom3(iseed);
  pmass = masses[imass];
  TFile htfile("HitTest.root","RECREATE");
// define the context
  UniformBField BF(1.0); // 1 Tesla
  Context context(BF);
  CVD2T d2t(sdrift,sigt*sigt);
  Vec4 origin(0.0,0.0,0.0,0.0);
  float sint = sqrt(1.0-cost*cost);
  Mom4 momv(mom*sint*cos(phi),mom*sint*sin(phi),mom*cost,pmass);
  LHelix lhel(origin,momv,icharge,context);
  PLHelix plhel(lhel);
  // truncate the range according to the Z range
  Vec3 vel; plhel.velocity(0.0,vel);
  plhel.range() = TRange(-0.5*zrange/vel.Z(),0.5*zrange/vel.Z());
// create Canvase
  TCanvas* hcan = new TCanvas("hcan","Hits",1000,1000);
  TPolyLine3D* hel = new TPolyLine3D(100);
  Vec4 hpos;
  double tstep = plhel.range().range()/100.0;
  for(int istep=0;istep<101;++istep){
  // compute the position from the time
    hpos.SetE(plhel.range().low() + tstep*istep);
    plhel.position(hpos);
    // add these positions to the TPolyLine3D
    hel->SetPoint(istep, hpos.X(), hpos.Y(), hpos.Z());
  }
  // draw the helix
  hel->SetLineColor(kBlue);
  hel->Draw();
// divide time range
  double dt = plhel.range().range()/(nhits-1);
  cout << "True " << plhel << endl;
  // create straw material
  MatDBInfo matdbinfo;
  const DetMaterial* wallmat = matdbinfo.findDetMaterial("straw-wall");
  const DetMaterial* gasmat = matdbinfo.findDetMaterial("straw-gas");
  const DetMaterial* wiremat = matdbinfo.findDetMaterial("straw-wire");

  const_cast<DetMaterial*>(wallmat)->setScatterFraction(scatfrac);
  const_cast<DetMaterial*>(gasmat)->setScatterFraction(scatfrac);
  const_cast<DetMaterial*>(wiremat)->setScatterFraction(scatfrac);

  StrawMat smat(rstraw,wthick,rwire, *wallmat, *gasmat, *wiremat);
  TGraph* ggplen = new TGraph(nhits); ggplen->SetTitle("Gas Pathlength;DOCA (mm);Pathlength (mm)");
  TGraph* gwplen = new TGraph(nhits); gwplen->SetTitle("Wall Pathlength;DOCA (mm);Pathlength (mm)");
  TGraph* ggeloss = new TGraph(nhits); ggeloss->SetTitle("Gas Energy Loss;DOCA (mm);Energy Loss (MeV)");
  TGraph* gweloss = new TGraph(nhits); gweloss->SetTitle("Wall Energy Loss;DOCA (mm);Energy Loss (MeV)");
  TGraph* ggscat = new TGraph(nhits); ggscat->SetTitle("Gas Scattering;DOCA (mm);Scattering (radians)");
  TGraph* gwscat = new TGraph(nhits); gwscat->SetTitle("Wall Scattering;DOCA (mm);Scattering (radians)");

  // generate hits
  std::vector<StrawHit> hits;
  std::vector<TPolyLine3D*> tpl;
  cout << "ambigdoca = " << ambigdoca << endl;
  for(size_t ihit=0; ihit<nhits; ihit++){
    double htime = plhel.range().low() + ihit*dt;
    auto tline = GenerateStraw(plhel,htime,TR);
//    cout << "TLine " << tline << endl;
// try to create TPoca for a hit
    TDPoca<PLHelix,TLine> tp(plhel,tline);
//    cout << "TPoca status " << tp.statusName() << " doca " << tp.doca() << " dt " << tp.deltaT() << endl;
// only set ambiguity if DOCA is above this value
    WireHit::LRAmbig ambig(WireHit::null);
    if(fabs(tp.doca())> ambigdoca)
      ambig = tp.doca() < 0 ? WireHit::left : WireHit::right;
    // construct the hit from this trajectory
    StrawHit sh(tline,context,d2t,smat,ambigdoca,ambig);
    hits.push_back(sh);
   // compute residual
    Residual res;
    sh.resid(tp,res);
    cout << res << " ambig " << ambig << endl;
    TPolyLine3D* line = new TPolyLine3D(2);
    Vec3 plow, phigh;
    tline.position(tline.range().low(),plow);
    tline.position(tline.range().high(),phigh);
    line->SetPoint(0,plow.X(),plow.Y(), plow.Z());
    line->SetPoint(1,phigh.X(),phigh.Y(), phigh.Z());
    line->SetLineColor(kRed);
    line->Draw();
    tpl.push_back(line);
    // compute material effects
//    double adot = tp.dirDot();
    double adot =0.0; // transverse
    double gpath = smat.gasPath(tp.doca(),ddoca,adot);
    double wpath = smat.wallPath(tp.doca(),ddoca,adot);
    ggplen->SetPoint(ihit,fabs(tp.doca()),gpath );
    gwplen->SetPoint(ihit,fabs(tp.doca()),wpath);
    cout << "doca " << tp.doca() << " gas path " << smat.gasPath(tp.doca(),ddoca,adot)
    << " wall path " << smat.wallPath(tp.doca(),ddoca,adot) << endl;

    // compute material effects
    double geloss = gasmat->energyLoss(mom,gpath,pmass);
    double weloss = wallmat->energyLoss(mom,wpath,pmass);
    double gscat = gasmat->scatterAngleRMS(mom,gpath,pmass);
    double wscat = wallmat->scatterAngleRMS(mom,wpath,pmass);
    ggeloss->SetPoint(ihit,fabs(tp.doca()),geloss);
    gweloss->SetPoint(ihit,fabs(tp.doca()),weloss);
    ggscat->SetPoint(ihit,fabs(tp.doca()),gscat);
    gwscat->SetPoint(ihit,fabs(tp.doca()),wscat);
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
  hcan->Write();
// test updating the hit residual and derivatives with different trajectories
  vector<double> delpars { 0.5, 0.1, 0.5, 0.5, 0.005, 5.0}; // small parameter changes for derivative calcs
  unsigned nsteps(10);
  vector<TGraph*> hderivg(LHelix::NParams());
  for(size_t ipar=0;ipar < LHelix::NParams();ipar++){
    auto tpar = static_cast<LHelix::ParamIndex>(ipar);
    hderivg[ipar] = new TGraph(hits.size()*nsteps);
    std::string title = LHelix::paramTitle(tpar) + " Residual Derivative Test;"
    + LHelix::paramName(tpar) + " Exact #Delta (mm);"
    + LHelix::paramName(tpar) + " Algebraic #Delta (mm)";
    hderivg[ipar]->SetTitle(title.c_str());
  }
  unsigned ipt(0);
  for(size_t istep=0;istep<nsteps; istep++){
    for(size_t ipar=0;ipar < LHelix::NParams();ipar++){
      double dpar = delpars[ipar]*(-0.5 + float(istep)/float(nsteps));
      // modify the helix
      LHelix modlhel =lhel;
      modlhel.params().parameters()[ipar] += dpar;
      PLHelix modplhel(modlhel);
      ROOT::Math::SVector<double,6> dpvec;
      dpvec[ipar] += dpar;
      // update the hits
      for(auto& hit : hits) {
	KKHit kkhit(hit,plhel);
	Residual ores = kkhit.refResid(); // original residual
	kkhit.update(modplhel);// refer to moded helix
	Residual mres = kkhit.refResid();
	double dr = ores.resid()-mres.resid(); // this sign is confusing.  I think
	// it means the fit needs to know how much to change the ref parameters, which is
	// opposite from how much the ref parameters are different from the measurement
	// compare the change with the expected from the derivatives
	kkhit.update(plhel);// refer back to original
	auto pder = kkhit.dRdP();
	double ddr = ROOT::Math::Dot(pder,dpvec);
	hderivg[ipar]->SetPoint(ipt++,dr,ddr);
      }
    }
  }
  TCanvas* hderivgc = new TCanvas("hderiv","hderiv",800,600);
  hderivgc->Divide(3,2);
  for(size_t ipar=0;ipar<6;++ipar){
    hderivgc->cd(ipar+1);
    hderivg[ipar]->Draw("A*");
  };
  hderivgc->Write();

  TCanvas* mateff = new TCanvas("mateff","mateff",800,600);
  mateff->Divide(3,2);
  mateff->cd(1);
  ggplen->Draw("A*");
  mateff->cd(2);
  ggeloss->Draw("A*");
  mateff->cd(3);
  ggscat->Draw("A*");
  mateff->cd(4);
  gwplen->Draw("A*");
  mateff->cd(5);
  gweloss->Draw("A*");
  mateff->cd(6);
  gwscat->Draw("A*");

  mateff->Write();

  htfile.Write();
  htfile.Close();
  exit(EXIT_SUCCESS);
}
