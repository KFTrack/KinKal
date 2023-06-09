#include "KinKal/Trajectory/ConstantDistanceToTime.cc"
#include "KinKal/Trajectory/DistanceToTime.hh"
#include <iostream>
#include <functional>

#include "TH1F.h"
#include "TSystem.h"
#include "THelix.h"
#include "TFile.h"
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
#include "TProfile.h"
#include "TProfile2D.h"

using namespace std;

TGraph* graph(int numIter, double start, double stepSize, DistanceToTime* d, function<double(double, DistanceToTime*)> fn) {
    double x[numIter];
    double y[numIter];

    for (int i = 0; i < numIter; i++) {
        x[i] = i * stepSize + start;
        y[i] = fn(x[i], d);
    }
    return new TGraph(numIter, x, y);
}

double timeWrapper(double x, DistanceToTime* d) {
    return d->time(x);
}

double distanceWrapper(double x, DistanceToTime* d) {
    return d->distance(x);
}

double inverseSpeedWrapper(double x, DistanceToTime* d) {
    return d->inverseSpeed(x);
}

double speedWrapper(double x, DistanceToTime* d) {
    return d->speed(x);
}

int main(int argc, char **argv) {
    // cout << "Hello World" << endl;
    ConstantDistanceToTime* d = new ConstantDistanceToTime(-136);
    
    TGraph* g1 = graph(200, 0, 1, d, &timeWrapper);
    g1->SetTitle("Distance->deltaT;Distance (mm);deltaT (ns)");
    TGraph* g2 = graph(200, 0, 1, d, &inverseSpeedWrapper);
    g2->SetTitle("Inverse Speed (mm/ns);Distance (mm);deltaT (ns)");
    TGraph* g3 = graph(200, 0, 1, d, &speedWrapper);
    g3->SetTitle("Speed (ns/mm); Distance (mm); deltaT (ns)");
    TGraph* g4 = graph(200, 0, 0.0099, d, &distanceWrapper);
    g4->SetTitle("DeltaT->Distance; deltaT (ns); Distance (mm)");

    TFile mefile("ConstantDistanceToTime.root","RECREATE");

    TCanvas* can = new TCanvas("ConstantDistanceToTime", "ConstantDistanceToTime", 1000, 1000);
    can->Divide(2, 2);
    can->cd(1);
    g1->Draw("AC*");
    can->cd(2);
    g2->Draw("AC*");
    can->cd(3);
    g3->Draw("AC*");
    can->cd(4);
    g4->Draw("AC*");
    can->Write();
    mefile.Write();
    mefile.Close();
}