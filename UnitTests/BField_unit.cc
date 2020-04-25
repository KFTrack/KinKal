// 
// test basic functions of BField class
//
#include "KinKal/LHelix.hh"
#include "KinKal/BField.hh"

#include <iostream>
#include <stdio.h>
#include <iostream>
#include <getopt.h>

#include "TH1F.h"
#include "TSystem.h"
#include "THelix.h"
#include "TPolyLine3D.h"
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
#include "TDirectory.h"

using namespace KinKal;
using namespace std;

void print_usage() {
  printf("Usage: BField  --momentum f --costheta f --charge i --bz f --dbz f --dby f --bgrad f\n");
}

int main(int argc, char **argv) {
  typedef LHelix KTRAJ;

  return 0;
}


