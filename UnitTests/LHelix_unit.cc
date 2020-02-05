// 
// test basic functions of the Loop Helix TTraj class
//
#include "BTrk/KinKal/LHelix.hh"
#include "BTrk/KinKal/Context.hh"

#include <iostream>
#include <stdio.h>
#include <iostream>
#include <getopt.h>

#include "TH1F.h"
#include "TSystem.h"
#include "THelix.h"
#include "TPolyLine3D.h"
#include "TArrow.h"
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

void print_usage() {
  printf("Usage: LHelix  --momentum f --costheta f --azimuth f --particle i --charge i --zorigin f --torigin --tmin ff--tmax f \n");
}

int main(int argc, char **argv) {
  int opt;
  float mom(105.0), cost(0.7), phi(0.5);
  float masses[5]={0.511,105.66,139.57, 493.68, 938.0};
  int imass(0), icharge(-1);
  float pmass, oz(0.0), ot(0.0);
  float tmin(-5.0), tmax(5.0);

  static struct option long_options[] = {
    {"momentum",     required_argument, 0, 'm' },
    {"costheta",     required_argument, 0, 'c'  },
    {"azimuth",     required_argument, 0, 'a'  },
    {"particle",     required_argument, 0, 'p'  },
    {"charge",     required_argument, 0, 'q'  },
    {"zorigin",     required_argument, 0, 'z'  },
    {"torigin",     required_argument, 0, 't'  },
    {"tmin",     required_argument, 0, 's'  },
    {"tmax",     required_argument, 0, 'e'  }


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
      case 'z' : oz = atof(optarg);
		 break;
      case 't' : ot = atof(optarg);
		 break;
      case 's' : tmin = atof(optarg);
		 break;
      case 'e' : tmax = atof(optarg);
		 break;
      default: print_usage(); 
	       exit(EXIT_FAILURE);
    }
  }

  pmass = masses[imass];

  printf("Testing LHelix with momentum = %f, costheta = %f, phi = %f, mass = %f, charge = %i, z = %f, t = %f \n",mom,cost,phi,pmass,icharge,oz,ot);
// define the context

  Context context;
  context.Bz_ = 1.0; // 1 Tesla
  Vec4 origin(0.0,0.0,oz,ot);
  float sint = sqrt(1.0-cost*cost);
  Mom4 momv(mom*sint*cos(phi),mom*sint*sin(phi),mom*cost,pmass);
  LHelix lhel(origin,momv,icharge,context);
  LHelix invlhel(origin,momv,icharge,context);
  // also inverse helix
  invlhel.invertCT();
  Mom4 testmom;
  lhel.momentum(0.0,testmom);
  cout << "Helix with momentum " << testmom << " position " << origin << " has parameters: " << endl;
  for(size_t ipar=0;ipar < LHelix::npars_;ipar++)
    cout << LHelix::paramName(static_cast<LHelix::paramIndex>(ipar) ) << " : " << lhel.param(ipar) << endl;
  Vec3 vel;
  lhel.velocity(ot,vel);
  double dot = vel.Dot(testmom)/KinKal::c_;
  cout << "velocity dot mom = " << dot << endl;
  cout << "momentum beta =" << momv.Beta() << " LHelix beta = " << lhel.beta() << endl;
  Vec3 mdir;
  lhel.direction(ot,mdir);

// now graph this as a polyline over the specified time range.
  double tstep = 0.1; // nanoseconds
  double trange = tmax-tmin;
  int nsteps = (int)rint(trange/tstep);
// create Canvase
  TCanvas* hcan = new TCanvas("hcan","Helix",1000,1000);
//TPolyLine to graph the result
  TPolyLine3D* hel = new TPolyLine3D(nsteps+1);
  TPolyLine3D* invhel = new TPolyLine3D(nsteps+1);
  Vec4 hpos;
  for(int istep=0;istep<nsteps+1;++istep){
  // compute the position from the time
    hpos.SetE(tmin + tstep*istep);
    lhel.position(hpos);
    // add these positions to the TPolyLine3D
    hel->SetPoint(istep, hpos.X(), hpos.Y(), hpos.Z());
    invlhel.position(hpos);
    invhel->SetPoint(istep, hpos.X(), hpos.Y(), hpos.Z());
  }
  // draw the helix
  if(icharge > 0)
    hel->SetLineColor(kBlue);
  else
    hel->SetLineColor(kRed);
  hel->Draw();
//  invhel->SetLineColor(kCyan);
//  invhel->Draw();
  // draw the origin and axes
  TAxis3D* rulers = new TAxis3D();
  rulers->GetXaxis()->SetAxisColor(kBlue);
  rulers->GetXaxis()->SetLabelColor(kBlue);
  rulers->GetYaxis()->SetAxisColor(kCyan);
  rulers->GetYaxis()->SetLabelColor(kCyan);
  rulers->GetZaxis()->SetAxisColor(kOrange);
  rulers->GetZaxis()->SetLabelColor(kOrange);
//   TAxis3D::ToggleRulers();
//   TAxis3D::ToggleRulers();
//  cout <<"axis = " << rulers << endl;
  rulers->Draw();

//  TEveArrow* start = new TEveArrow(pos._x,pos._y,pos._z,0.0,0.0,0.0);
//  start->Draw();

//  TVector3* start =new TVector3(pos._x,pos._y,pos._z);
//  start->Draw();
  TPolyLine3D* refmom = new TPolyLine3D(2);
  refmom->SetPoint(0,origin.X(),origin.Y(),origin.Z());
  Vec3 lend = origin.Vect() + lhel.rad()*mdir;
  refmom->SetPoint(1,lend.X(),lend.Y(),lend.Z());
  int momcolor;
  if(icharge>0)
    momcolor = kRed;
  else
    momcolor = kBlack;
  refmom->SetLineColor(momcolor);
  refmom->Draw();

  TPolyMarker3D* refp = new TPolyMarker3D(1,24);
  refp->SetPoint(0,origin.X(),origin.Y(),origin.Z());
  refp->Draw();
  TPolyMarker3D* refmomp = new TPolyMarker3D(1,22);
  refmomp->SetPoint(1,lend.X(),lend.Y(),lend.Z());
  refmomp->SetMarkerColor(momcolor);
  refmomp->Draw();

  TPolyMarker3D* helixp = new TPolyMarker3D(1,3);
  helixp->SetMarkerColor(kGreen);
  hpos.SetE(ot);
  lhel.position(hpos);
  helixp->SetPoint(0,hpos.X(),hpos.Y(),hpos.Z());
  helixp->Draw();

  TPolyMarker3D* testp = new TPolyMarker3D(1,25);
  testp->SetMarkerColor(kOrange);
  testp->SetPoint(0,hpos.X(),hpos.Y(),hpos.Z()+2*M_PI*lhel.lam());
  testp->Draw();

  TPolyMarker3D* startp = new TPolyMarker3D(1,21);
  startp->SetMarkerColor(kBlue);
  hpos.SetE(tmin);
  lhel.position(hpos);
  startp->SetPoint(0,hpos.X(),hpos.Y(),hpos.Z());
  startp->Draw();

  TPolyMarker3D* endp = new TPolyMarker3D(1,22);
  endp->SetMarkerColor(kBlue);
  hpos.SetE(tmax);
  lhel.position(hpos);
  endp->SetPoint(0,hpos.X(),hpos.Y(),hpos.Z());
  endp->Draw();

  TLegend* leg = new TLegend(0.8,0.8,1.0,1.0);
  char title[80];
  snprintf(title,80,"Helix, mass=%3.1g MeV/c^{2}, q=%1i",pmass,icharge);
  leg->AddEntry(hel,title,"L");
  snprintf(title,80,"Initial Momentum =%3.1g MeV/c",mom);
  leg->AddEntry(refmom,title,"L");
  leg->AddEntry(refp,"Initial Position","P");
  snprintf(title,80,"Helix, t=%4.2g ns",ot);
  leg->AddEntry(helixp,title,"P");
  snprintf(title,80,"Helix, t=%4.2g ns",ot+tmin);
  leg->AddEntry(startp,title,"P");
  snprintf(title,80,"Helix, t=%4.2g ns",ot+tmax);
  leg->AddEntry(endp,title,"P");
  leg->Draw();

  snprintf(title,80,"LHelix_m%3.1f_p%3.2f_q%i.root",pmass,mom,icharge);
  cout << "Saving canvas to " << title << endl;
  hcan->SaveAs(title); 

  exit(EXIT_SUCCESS);
}

