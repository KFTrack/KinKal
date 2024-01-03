//
// ToyMC test of hits
//
#include "KinKal/Trajectory/ParticleTrajectory.hh"
#include "KinKal/Trajectory/LoopHelix.hh"
#include "KinKal/Trajectory/ClosestApproach.hh"
#include "KinKal/Examples/SimpleWireHit.hh"
#include "KinKal/Examples/ScintHit.hh"
#include "KinKal/Detector/StrawMaterial.hh"
#include "KinKal/Detector/Residual.hh"
#include "KinKal/General/BFieldMap.hh"
#include "KinKal/General/Vectors.hh"
#include "KinKal/General/PhysicalConstants.h"
#include "KinKal/Tests/ToyMC.hh"

#include <iostream>
#include <cstdio>
#include <iostream>
#include <getopt.h>

#include "TH1F.h"
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

void print_usage() {
  printf("Usage: HitTest  --momentum f --particle i --charge i --strawhit i --scinthit i --zrange f --nhits i --hres f --seed i --ambigdoca f --By f --Bgrad f --simmat_ i --prec f\n");
}

template <class KTRAJ>
int HitTest(int argc, char **argv, const vector<double>& delpars) {
  using PTRAJ = ParticleTrajectory<KTRAJ>;
  using HIT = Hit<KTRAJ>;
  using HITPTR = std::shared_ptr<HIT>;
  using HITCOL = vector<HITPTR>;
  using EXING = ElementXing<KTRAJ>;
  using EXINGPTR = std::shared_ptr<EXING>;
  using EXINGCOL = std::vector<EXINGPTR>;
  using WIREHIT = SimpleWireHit<KTRAJ>;
  using WIREHITPTR = std::shared_ptr<WIREHIT>;
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
  double Bgrad(0.0), By(0.0);
  bool simmat_(true), scinthit_(true), strawhit_(true);
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
  KKTest::ToyMC<KTRAJ> toy(*BF, mom, icharge, zrange, iseed, nhits, simmat_, scinthit_,ambigdoca, pmass );
  toy.setInefficiency(0.0);
  PTRAJ tptraj;
  //  cout << "True " << tptraj << endl;
  StrawMaterial const& smat = toy.strawMaterial();
  TGraph* ggplen = new TGraph(nhits); ggplen->SetTitle("Gas Pathlength;Doca (mm);Pathlength (mm)"); ggplen->SetMinimum(0.0);
  TGraph* gwplen = new TGraph(nhits); gwplen->SetTitle("Wall Pathlength;Doca (mm);Pathlength (mm)"); gwplen->SetMinimum(0.0);
  TGraph* ggeloss = new TGraph(nhits); ggeloss->SetTitle("Gas Energy Change;Pathlength (mm);Energy Change (MeV)"); ggeloss->SetMaximum(0.0);
  TGraph* gweloss = new TGraph(nhits); gweloss->SetTitle("Wall Energy Change;Pathlength (mm);Energy Change (MeV)"); gweloss->SetMaximum(0.0);
  TGraph* ggscat = new TGraph(nhits); ggscat->SetTitle("Gas Scattering;Pathlength (mm);Scattering (radians)"); ggscat->SetMinimum(0.0);
  TGraph* gwscat = new TGraph(nhits); gwscat->SetTitle("Wall Scattering;Pathlength (mm);Scattering (radians)"); gwscat->SetMinimum(0.0);
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
  StrawXingConfig sxconfig(0.3,5.0,10.0,false);
  for(auto& thit : thits) {
    Residual res;
    ClosestApproachData tpdata;
    WIREHIT* strawhit = dynamic_cast<WIREHIT*>(thit.get());
    SCINTHIT* scinthit = dynamic_cast<SCINTHIT*>(thit.get());
    if(strawhit && strawhit_){
      res = strawhit->refResidual(0);
      tpdata = strawhit->closestApproach().tpData();
    } else if(scinthit && scinthit_){
      res = scinthit->refResidual(0);
      tpdata = scinthit->closestApproach().tpData();
    } else
      continue;
    TPolyLine3D* line = new TPolyLine3D(2);
    VEC3 plow, phigh;
    WIREHITPTR shptr = std::dynamic_pointer_cast<WIREHIT> (thit);
    SCINTHITPTR lhptr = std::dynamic_pointer_cast<SCINTHIT> (thit);
    if((bool)shptr){
      auto const& tline = shptr->wire();
      plow = tline.start();
      phigh = tline.end();
      line->SetLineColor(kRed);
    } else if ((bool)lhptr){
      auto const& tline = lhptr->sensorAxis();
      plow = tline.start();
      phigh = tline.end();
      line->SetLineColor(kCyan);
    }
    line->SetPoint(0,plow.X(),plow.Y(), plow.Z());
    line->SetPoint(1,phigh.X(),phigh.Y(), phigh.Z());
    line->Draw();
    tpl.push_back(line);
    // compute material effects
    if(strawhit && strawhit_){
      double doca = fabs(tpdata.doca());
      double gaspath, wirepath, wallpath;
      smat.pathLengths(tpdata,sxconfig,wallpath,gaspath,wirepath);
      ggplen->SetPoint(ihit,doca,gaspath);
      gwplen->SetPoint(ihit,doca,wallpath);
      // compute material effects
      double geloss = smat.gasMaterial().energyLoss(mom,gaspath,pmass);
      double weloss = smat.wallMaterial().energyLoss(mom,wallpath,pmass);
      double gscat = smat.gasMaterial().scatterAngleRMS(mom,gaspath,pmass);
      double wscat = smat.wallMaterial().scatterAngleRMS(mom,wallpath,pmass);
      ggeloss->SetPoint(ihit,gaspath,geloss);
      gweloss->SetPoint(ihit,wallpath,weloss);
      ggscat->SetPoint(ihit,gaspath,gscat);
      gwscat->SetPoint(ihit,wallpath,wscat);
    }
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
  MetaIterConfig miconfig;
  for(auto& thit : thits) {
    thit->updateState(miconfig,false);
    Residual ores;
    ClosestApproachData tpdata;
    WIREHIT* strawhit = dynamic_cast<WIREHIT*>(thit.get());
    SCINTHIT* scinthit = dynamic_cast<SCINTHIT*>(thit.get());
    if(strawhit && strawhit_){
      ores = strawhit->refResidual(0);
      tpdata = strawhit->closestApproach().tpData();
    } else if(scinthit && scinthit_){
      ores = scinthit->refResidual(0);
      tpdata = scinthit->closestApproach().tpData();
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
          KTRAJ modktraj = tptraj.nearestPiece(thit->time());
          modktraj.params().parameters()[ipar] += dpar;
          PTRAJ modtptraj(modktraj);
          KinKal::DVEC dpvec;
          dpvec[ipar] = dpar;
          thit->updateReference(modtptraj.backPtr());// refer to moded helix
          thit->updateState(miconfig,false);
          Residual mres;
          if(strawhit){
            mres = strawhit->refResidual(WIREHIT::dresid);
          } else if(scinthit) {
            mres = scinthit->refResidual(0);
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
    if(fabs(hderivg[ipar]->GetRMS(1)) > 1e-5 && fabs(hderivg[ipar]->GetRMS(2)) > 1e-5){
      pline->SetParameters(0.0,1.0);
      TFitResultPtr pfitr = hderivg[ipar]->Fit(pline,"SQ","AC*");
      if(fabs(pfitr->Parameter(0))> 100*delpars[ipar] || fabs(pfitr->Parameter(1)-1.0) > 1e-2){
        cout << "Parameter "
          << KTRAJ::paramName(typename KTRAJ::ParamIndex(ipar))
          << " Residual derivative Out of tolerance : Offset " << pfitr->Parameter(0) << " Slope " << pfitr->Parameter(1) << endl;
        status = 1;
      }
    } else {
      cout << "Zero derivatives for parameter " << ipar << endl;
    }
    hderivgc->cd(ipar+1);
    hderivg[ipar]->Draw("AC*");
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
