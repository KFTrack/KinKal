//
// ToyMC test of hits
//
#include "MatEnv/MatDBInfo.hh"
#include "MatEnv/DetMaterial.hh"
#include "KinKal/PKTraj.hh"
#include "KinKal/LHelix.hh"
#include "KinKal/TLine.hh"
#include "KinKal/TPoca.hh"
#include "KinKal/StrawHit.hh"
#include "KinKal/ScintHit.hh"
#include "KinKal/StrawMat.hh"
#include "KinKal/KKMHit.hh"
#include "KinKal/Residual.hh"
#include "KinKal/BField.hh"
#include "KinKal/Vectors.hh"
#include "KinKal/KKHit.hh"
#include "CLHEP/Units/PhysicalConstants.h"
#include "UnitTests/ToyMC.hh"

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
#include "TH2F.h"
#include "TDirectory.h"
#include "TProfile.h"
#include "TProfile2D.h"

using namespace MatEnv;
using namespace KinKal;
using namespace std;
// avoid confusion with root
using KinKal::TLine;

void print_usage() {
  printf("Usage: HitTest  --momentum f --costheta f --azimuth f --particle i --charge i --lighthit i --zrange f --nhits i --hres f --seed i --ambigdoca f --ddoca f --By f --Bgrad f --simmat i\n");
}

template <class KTRAJ>
int HitTest(int argc, char **argv, const vector<double>& delpars) {
  typedef KinKal::PKTraj<KTRAJ> PKTRAJ;
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
  typedef vector<DXINGPTR> DXINGCOL;
  typedef Residual<KTRAJ::NParams()> RESIDUAL;

  int opt;
  double mom(105.0);
  double masses[5]={0.511,105.66,139.57, 493.68, 938.0};
  int imass(0), icharge(-1);
  double pmass;
  unsigned nhits(40);
  int iseed(124223);
  double ddoca(0.1);
  double Bgrad(0.0), By(0.0);
  bool simmat(true), lighthit(true);
  double zrange(3000.0); // tracker dimension

  static struct option long_options[] = {
    {"momentum",     required_argument, 0, 'm' },
    {"simmat",     required_argument, 0, 'b'  },
    {"particle",     required_argument, 0, 'p'  },
    {"charge",     required_argument, 0, 'q'  },
    {"zrange",     required_argument, 0, 'z'  },
    {"seed",     required_argument, 0, 's'  },
    {"hres",     required_argument, 0, 'h'  },
    {"lighthit",     required_argument, 0, 'l'  },
    {"nhits",     required_argument, 0, 'n'  },
    {"ddoca",     required_argument, 0, 'x'  },
    {"By",     required_argument, 0, 'y'  },
    {"Bgrad",     required_argument, 0, 'g'  },
    {NULL, 0,0,0}
  };

  int long_index =0;
  while ((opt = getopt_long_only(argc, argv,"",
	  long_options, &long_index )) != -1) {
    switch (opt) {
      case 'm' : mom = atof(optarg);
		 break;
      case 'p' : imass = atoi(optarg);
		 break;
      case 'q' : icharge = atoi(optarg);
		 break;
      case 'z' : zrange = atof(optarg);
		 break;
      case 'n' : nhits = atoi(optarg);
		 break;
      case 'l' : lighthit = atoi(optarg);
		 break;
      case 'b' : simmat = atoi(optarg);
		 break;
      case 's' : iseed = atoi(optarg);
		 break;
      case 'x' : ddoca = atof(optarg);
		 break;
      case 'y' : By = atof(optarg);
		 break;
      case 'g' : Bgrad = atof(optarg);
		 break;
      default: print_usage();
	       exit(EXIT_FAILURE);
    }
  }

  pmass = masses[imass];
  TFile htfile((KTRAJ::trajName()+"HitTest.root").c_str(),"RECREATE");
  // construct BField
  Vec3 bnom(0.0,By,1.0);
  BField* BF;
  if(Bgrad != 0){
    BF = new GradBField(1.0-0.5*Bgrad,1.0+0.5*Bgrad,-0.5*zrange,0.5*zrange);
    bnom = BF->fieldVect(Vec3(0.0,0.0,0.0));
  } else {
    BF = new UniformBField(bnom);
  }
  KKTest::ToyMC<KTRAJ> toy(*BF, mom, icharge, zrange, iseed, nhits, simmat, lighthit, -1.0, pmass );
  PKTRAJ tptraj;
//  cout << "True " << tptraj << endl;
  StrawMat const& smat = toy.strawMaterial();
  TGraph* ggplen = new TGraph(nhits); ggplen->SetTitle("Gas Pathlength;Doca (mm);Pathlength (mm)"); ggplen->SetMinimum(0.0);
  TGraph* gwplen = new TGraph(nhits); gwplen->SetTitle("Wall Pathlength;Doca (mm);Pathlength (mm)"); gwplen->SetMinimum(0.0);
  TGraph* ggeloss = new TGraph(nhits); ggeloss->SetTitle("Gas Energy Change;Doca (mm);Energy Change (MeV)"); ggeloss->SetMaximum(0.0);
  TGraph* gweloss = new TGraph(nhits); gweloss->SetTitle("Wall Energy Change;Doca (mm);Energy Change (MeV)"); gweloss->SetMaximum(0.0);
  TGraph* ggscat = new TGraph(nhits); ggscat->SetTitle("Gas Scattering;Doca (mm);Scattering (radians)"); ggscat->SetMinimum(0.0);
  TGraph* gwscat = new TGraph(nhits); gwscat->SetTitle("Wall Scattering;Doca (mm);Scattering (radians)"); gwscat->SetMinimum(0.0);
  std::vector<TPolyLine3D*> tpl;
  // generate hits
  THITCOL thits;
  DXINGCOL dxings; // this program shares det xing ownership with KKTrk
  toy.simulateParticle(tptraj, thits, dxings);
// create Canvas
  TCanvas* hcan = new TCanvas("hcan","Hits",1000,1000);
  TPolyLine3D* hel = new TPolyLine3D(100);
  Vec4 hpos;
  double tstep = tptraj.range().range()/100.0;
  for(int istep=0;istep<101;++istep){
  // compute the position from the time
    hpos.SetE(tptraj.range().low() + tstep*istep);
    tptraj.position(hpos);
    // add these positions to the TPolyLine3D
    hel->SetPoint(istep, hpos.X(), hpos.Y(), hpos.Z());
  }
  // draw the helix
  hel->SetLineColor(kBlue);
  hel->Draw();
  unsigned ihit(0);
  for(auto const& thit : thits) {
   // compute residual
    RESIDUAL  res;
    thit->resid(tptraj,res);
    TPolyLine3D* line = new TPolyLine3D(2);
    Vec3 plow, phigh;
    STRAWHITPTR shptr = std::dynamic_pointer_cast<STRAWHIT> (thit); 
    SCINTHITPTR lhptr = std::dynamic_pointer_cast<SCINTHIT> (thit);
    if(shptr.use_count() > 0){
      auto const& tline = shptr->wire();
      plow = tline.position(tline.range().low());
      phigh = tline.position(tline.range().high());
      line->SetLineColor(kRed);
    } else if (lhptr.use_count() > 0){
      auto const& tline = lhptr->sensorAxis();
      plow = tline.position(tline.range().low());
      phigh = tline.position(tline.range().high());
      line->SetLineColor(kCyan);
    }
    line->SetPoint(0,plow.X(),plow.Y(), plow.Z());
    line->SetPoint(1,phigh.X(),phigh.Y(), phigh.Z());
    line->Draw();
    tpl.push_back(line);
    // compute material effects
//    double adot = tp.dirDot();
    double adot =0.0; // transverse
//    double adot = tp.dirDot();
    double doca = fabs(res.tPoca().doca());
    double gpath = smat.gasPath(doca,ddoca,adot);
    double wpath = smat.wallPath(doca,ddoca,adot);
    ggplen->SetPoint(ihit,doca,gpath );
    gwplen->SetPoint(ihit,doca,wpath);
//    cout << "Residual " << res << " gas path " << smat.gasPath(doca,ddoca,adot)
//    << " wall path " << smat.wallPath(doca,ddoca,adot) << endl;
    // compute material effects
    double geloss = smat.gasMaterial().energyLoss(mom,gpath,pmass);
    double weloss = smat.wallMaterial().energyLoss(mom,wpath,pmass);
    double gscat = smat.gasMaterial().scatterAngleRMS(mom,gpath,pmass);
    double wscat = smat.wallMaterial().scatterAngleRMS(mom,wpath,pmass);
    ggeloss->SetPoint(ihit,doca,geloss);
    gweloss->SetPoint(ihit,doca,weloss);
    ggscat->SetPoint(ihit,doca,gscat);
    gwscat->SetPoint(ihit,doca,wscat);
    ihit++;
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
  unsigned nsteps(10);
  vector<TGraph*> hderivg(KTRAJ::NParams());
  for(size_t ipar=0;ipar < KTRAJ::NParams();ipar++){
    auto tpar = static_cast<typename KTRAJ::ParamIndex>(ipar);
    hderivg[ipar] = new TGraph(thits.size()*nsteps);
    std::string title = KTRAJ::paramTitle(tpar) + " Residual Derivative Test;"
    + KTRAJ::paramName(tpar) + " Exact #Delta (mm);"
    + KTRAJ::paramName(tpar) + " Algebraic #Delta (mm)";
    hderivg[ipar]->SetTitle(title.c_str());
  }
  unsigned ipt(0);
  for(size_t istep=0;istep<nsteps; istep++){
    for(size_t ipar=0;ipar < KTRAJ::NParams();ipar++){
      auto tpar = static_cast<typename KTRAJ::ParamIndex>(ipar);
      double dpar = delpars[ipar]*(-0.5 + double(istep)/double(nsteps));
      // update the hits
      for(auto& thit : thits) {
        KKHit kkhit(thit,tptraj);
        RESIDUAL ores = kkhit.refResid(); // original residual
            // modify the helix
        KTRAJ modktraj = tptraj.nearestPiece(kkhit.time());
        modktraj.params().parameters()[ipar] += dpar;
        PKTRAJ modtptraj(modktraj);
        ROOT::Math::SVector<double,6> dpvec;
        dpvec[ipar] += dpar;
        kkhit.update(modtptraj);// refer to moded helix
        RESIDUAL mres = kkhit.refResid();
        double dr = ores.value()-mres.value(); // this sign is confusing.  I think
        // it means the fit needs to know how much to change the ref parameters, which is
        // opposite from how much the ref parameters are different from the measurement
        // compare the change with the expected from the derivatives
        kkhit.update(tptraj);// refer back to original
        auto pder = ores.dRdP();
        double ddr = ROOT::Math::Dot(pder,dpvec);
        hderivg[ipar]->SetPoint(ipt++,dr,ddr);
        if(dr*ddr < 0.0)cout << "Sign error " << KTRAJ::paramName(tpar) << " hit " << *thit
        << " doca " << ores.tPoca().doca() << " DirDot " << ores.tPoca().dirDot() <<" Exact change " << dr << " deriv " << ddr << endl;
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
