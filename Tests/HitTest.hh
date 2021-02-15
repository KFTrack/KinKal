//
// ToyMC test of hits
//
#include "KinKal/MatEnv/MatDBInfo.hh"
#include "KinKal/MatEnv/DetMaterial.hh"
#include "KinKal/Trajectory/ParticleTrajectory.hh"
#include "KinKal/Trajectory/LoopHelix.hh"
#include "KinKal/Trajectory/Line.hh"
#include "KinKal/Trajectory/ClosestApproach.hh"
#include "KinKal/Tests/SimpleWireHit.hh"
#include "KinKal/Tests/ScintHit.hh"
#include "KinKal/Detector/StrawMaterial.hh"
#include "KinKal/Detector/Residual.hh"
#include "KinKal/Detector/BFieldMap.hh"
#include "KinKal/General/Vectors.hh"
#include "KinKal/Fit/HitConstraint.hh"
#include "KinKal/General/PhysicalConstants.h"
#include "KinKal/Tests/ToyMC.hh"

#include <iostream>
#include <stdio.h>
#include <iostream>
#include <getopt.h>

#include "TH1F.h"
#include "TSystem.h"
#include "THelix.h"
#include "TPolyLine3D.h"
#include "TFile.h"
#include "TF1.h"
#include "TFitResult.h"
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
using KinKal::Line;

void print_usage() {
  printf("Usage: HitTest  --momentum f --particle i --charge i --strawhit i --scinthit i --zrange f --nhits i --hres f --seed i --ambigdoca f --ddoca f --By f --Bgrad f --simmat_ i --prec f\n");
}

template <class KTRAJ>
int HitTest(int argc, char **argv, const vector<double>& delpars) {
  using PKTRAJ = ParticleTrajectory<KTRAJ>;
  using KKHIT = HitConstraint<KTRAJ>;
  using HIT = Hit<KTRAJ>;
  using HITPTR = std::shared_ptr<HIT>;
  using HITCOL = vector<HITPTR>;
  using EXING = ElementXing<KTRAJ>;
  using EXINGPTR = std::shared_ptr<EXING>;
  using EXINGCOL = std::vector<EXINGPTR>;
  using STRAWHIT = WireHit<KTRAJ>;
  using STRAWHITPTR = std::shared_ptr<STRAWHIT>;
  using SCINTHIT = ScintHit<KTRAJ>;
  using SCINTHITPTR = std::shared_ptr<SCINTHIT>;
  using STRAWXING = StrawXing<KTRAJ>;
  using STRAWXINGPTR = shared_ptr<STRAWXING>;

  int status = 0;

  int opt;
  double mom(105.0);
  double masses[5]={0.511,105.66,139.57, 493.68, 938.0};
  double ambigdoca(-1.0);// minimum doca to set ambiguity, default sets for all hits
  int imass(0), icharge(-1);
  double pmass;
  unsigned nhits(40);
  int iseed(124223);
  double ddoca(0.1);
  double Bgrad(0.0), By(0.0);
  bool simmat_(true), scinthit_(true), strawhit_(true);
  double precision(1e-8);
  double zrange(3000.0); // tracker dimension

  static struct option long_options[] = {
    {"momentum",     required_argument, 0, 'm' },
    {"simmat_",     required_argument, 0, 'b'  },
    {"particle",     required_argument, 0, 'p'  },
    {"charge",     required_argument, 0, 'q'  },
    {"zrange",     required_argument, 0, 'z'  },
    {"seed",     required_argument, 0, 's'  },
    {"hres",     required_argument, 0, 'h'  },
    {"scinthit",     required_argument, 0, 'l'  },
    {"strawhit",     required_argument, 0, 'S'  },
    {"ambigdoca",     required_argument, 0, 'd'  },
    {"nhits",     required_argument, 0, 'n'  },
    {"ddoca",     required_argument, 0, 'x'  },
    {"By",     required_argument, 0, 'y'  },
    {"Bgrad",     required_argument, 0, 'g'  },
    {"prec",     required_argument, 0, 'P'  },
    {NULL, 0,0,0}
  };

  int long_index =0;
  while ((opt = getopt_long_only(argc, argv,"",
	  long_options, &long_index )) != -1) {
    switch (opt) {
      case 'm' : mom = atof(optarg);
		 break;
      case 'P' : precision = atof(optarg);
		 break;
      case 'p' : imass = atoi(optarg);
		 break;
      case 'q' : icharge = atoi(optarg);
		 break;
      case 'z' : zrange = atof(optarg);
		 break;
      case 'n' : nhits = atoi(optarg);
		 break;
      case 'l' : scinthit_ = atoi(optarg);
		 break;
      case 'S' : strawhit_ = atoi(optarg);
		 break;
      case 'd' : ambigdoca = atof(optarg);
		 break;
      case 'b' : simmat_ = atoi(optarg);
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
  // construct BFieldMap
  VEC3 bnom(0.0,By,1.0);
  BFieldMap* BF;
  if(Bgrad != 0){
    BF = new GradientBFieldMap(1.0-0.5*Bgrad,1.0+0.5*Bgrad,-0.5*zrange,0.5*zrange);
    bnom = BF->fieldVect(VEC3(0.0,0.0,0.0));
  } else {
    BF = new UniformBFieldMap(bnom);
  }
  KKTest::ToyMC<KTRAJ> toy(*BF, mom, icharge, zrange, iseed, nhits, simmat_, scinthit_,false, ambigdoca, pmass );
  toy.setInefficiency(0.0);
  PKTRAJ tptraj;
//  cout << "True " << tptraj << endl;
  StrawMaterial const& smat = toy.strawMaterial();
  TGraph* ggplen = new TGraph(nhits); ggplen->SetTitle("Gas Pathlength;Doca (mm);Pathlength (mm)"); ggplen->SetMinimum(0.0);
  TGraph* gwplen = new TGraph(nhits); gwplen->SetTitle("Wall Pathlength;Doca (mm);Pathlength (mm)"); gwplen->SetMinimum(0.0);
  TGraph* ggeloss = new TGraph(nhits); ggeloss->SetTitle("Gas Energy Change;Doca (mm);Energy Change (MeV)"); ggeloss->SetMaximum(0.0);
  TGraph* gweloss = new TGraph(nhits); gweloss->SetTitle("Wall Energy Change;Doca (mm);Energy Change (MeV)"); gweloss->SetMaximum(0.0);
  TGraph* ggscat = new TGraph(nhits); ggscat->SetTitle("Gas Scattering;Doca (mm);Scattering (radians)"); ggscat->SetMinimum(0.0);
  TGraph* gwscat = new TGraph(nhits); gwscat->SetTitle("Wall Scattering;Doca (mm);Scattering (radians)"); gwscat->SetMinimum(0.0);
  std::vector<TPolyLine3D*> tpl;
  // generate hits
  HITCOL thits;
  EXINGCOL dxings; // this program shares det xing ownership with Track
  toy.simulateParticle(tptraj, thits, dxings);
// create Canvas
  TCanvas* hcan = new TCanvas("hcan","Hits",1000,1000);
  TPolyLine3D* hel = new TPolyLine3D(100);
  double tstep = tptraj.range().range()/100.0;
  for(int istep=0;istep<101;++istep){
  // compute the position from the time
    VEC4 hpos = tptraj.position4(tptraj.range().begin() + tstep*istep);
    // add these positions to the TPolyLine3D
    hel->SetPoint(istep, hpos.X(), hpos.Y(), hpos.Z());
  }
  // draw the helix
  hel->SetLineColor(kBlue);
  hel->Draw();
  unsigned ihit(0);
  for(auto& thit : thits) {
    Residual res;
    ClosestApproachData tpdata;
    STRAWHIT* strawhit = dynamic_cast<STRAWHIT*>(thit.get());
    SCINTHIT* scinthit = dynamic_cast<SCINTHIT*>(thit.get());
    if(strawhit && strawhit_){
      strawhit->update(tptraj);
      res = strawhit->residual(0);
      tpdata = strawhit->closestApproach();
    } else if(scinthit && scinthit_){
      scinthit->update(tptraj);
      res = scinthit->residual(0);
      tpdata = scinthit->closestApproach();
    } else
      continue;
    TPolyLine3D* line = new TPolyLine3D(2);
    VEC3 plow, phigh;
    STRAWHITPTR shptr = std::dynamic_pointer_cast<STRAWHIT> (thit); 
    SCINTHITPTR lhptr = std::dynamic_pointer_cast<SCINTHIT> (thit);
    if((bool)shptr){
      auto const& tline = shptr->wire();
      plow = tline.position3(tline.range().begin());
      phigh = tline.position3(tline.range().end());
      line->SetLineColor(kRed);
    } else if ((bool)lhptr){
      auto const& tline = lhptr->sensorAxis();
      plow = tline.position3(tline.range().begin());
      phigh = tline.position3(tline.range().end());
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
    double doca = fabs(tpdata.doca());
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
  vector<TGraph*> hderivg(NParams());
  for(size_t ipar=0;ipar < NParams();ipar++){
    auto tpar = static_cast<typename KTRAJ::ParamIndex>(ipar);
    hderivg[ipar] = new TGraph(thits.size()*nsteps);
    std::string title = KTRAJ::paramName(tpar) + " Residual Derivative Test;Exact #Delta Residual (mm);Algebraic #Delta Residual (mm)";
    hderivg[ipar]->SetTitle(title.c_str());
  }
  unsigned ipt(0);
//  cout << tptraj << endl;
  for(auto& thit : thits) {
    KKHIT kkhit(thit,tptraj,precision);
    Residual ores;
    ClosestApproachData tpdata;
    STRAWHIT* strawhit = dynamic_cast<STRAWHIT*>(thit.get());
    SCINTHIT* scinthit = dynamic_cast<SCINTHIT*>(thit.get());
    if(strawhit && strawhit_){
      ores = strawhit->residual(0);
      tpdata = strawhit->closestApproach();
    } else if(scinthit && scinthit_){
      ores = scinthit->residual(0);
      tpdata = scinthit->closestApproach();
    } else
      continue;
    auto pder = ores.dRdP();
    for(size_t ipar=0;ipar < NParams();ipar++){
      if(fabs(pder[ipar])>1e-10){
	auto tpar = static_cast<typename KTRAJ::ParamIndex>(ipar);
	// update the hits
	for(size_t istep=0;istep<nsteps; istep++){
	  double dpar = delpars[ipar]*(-0.5 + double(istep)/double(nsteps));
	  // modify the helix
	  KTRAJ modktraj = tptraj.nearestPiece(kkhit.time());
	  modktraj.params().parameters()[ipar] += dpar;
	  PKTRAJ modtptraj(modktraj);
	  KinKal::DVEC dpvec;
	  dpvec[ipar] = dpar;
	  kkhit.update(modtptraj);// refer to moded helix
	  Residual mres;
	  if(strawhit){
	    mres = strawhit->residual(0);
	  } else if(scinthit) {
	    mres = scinthit->residual(0);
	  }
	  double dr = ores.value()-mres.value(); // this sign is confusing.  I think
	  // it means the fit needs to know how much to change the ref parameters, which is
	  // opposite from how much the ref parameters are different from the measurement
	  // compare the change with the expected from the derivatives
	  double ddr = ROOT::Math::Dot(pder,dpvec);
	  hderivg[ipar]->SetPoint(ipt++,dr,ddr);
	  if(fabs(dr - ddr) > 1.0 ){
	    cout << "Large ddiff " << KTRAJ::paramName(tpar) << " " << *thit << " delta " << dpar 
	      << " doca " << tpdata.doca() << " DirDot " << tpdata.dirDot() <<" Exact change " << dr << " deriv " << ddr << endl;
	    status = 2;
	  }
	}
      }
    }
  }
  // test
  TF1* pline = new TF1("pline","[0]+[1]*x");

  TCanvas* hderivgc = new TCanvas("hderiv","hderiv",800,600);
  hderivgc->Divide(3,2);
  for(size_t ipar=0;ipar<6;++ipar){
    if(fabs(hderivg[ipar]->GetMaximum()-hderivg[ipar]->GetMinimum())>1e-10){
      pline->SetParameters(0.0,1.0);
      hderivgc->cd(ipar+1);
      TFitResultPtr pfitr = hderivg[ipar]->Fit(pline,"SQ","AC*");
      hderivg[ipar]->Draw("AC*");
      if(fabs(pfitr->Parameter(0))> delpars[ipar] || fabs(pfitr->Parameter(1)-1.0) > 1e-2){
	cout << "Parameter " 
	  << KTRAJ::paramName(typename KTRAJ::ParamIndex(ipar))
	  << " Residual derivative Out of tolerance : Offset " << pfitr->Parameter(0) << " Slope " << pfitr->Parameter(1) << endl;
	status = 1;
      }
    }
  }
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
  cout << "Return status = " << status << endl;  
  exit(status);
}
