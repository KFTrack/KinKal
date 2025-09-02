#include "KinKal/Examples/CaloDistanceToTime.cc"
#include "KinKal/Trajectory/DistanceToTime.hh"
#include "KinKal/Trajectory/ConstantDistanceToTime.hh"
#include <iostream>
#include <functional>

#include "TH1F.h"
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

TGraph* graph(size_t numIter, double start, double stepSize, DistanceToTime* d, function<double(double, DistanceToTime*)> fn) {
    std::vector<double> x(numIter,0);
    std::vector<double> y(numIter,0);
    for (size_t i = 0; i < numIter; i++) {
        x[i] = i * stepSize + start;
        y[i] = fn(x[i], d);
    }
    return new TGraph(numIter, x.data(), y.data());
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
    double static const calorimeterLength = 200;
    CaloDistanceToTime* d = new CaloDistanceToTime(85.76, calorimeterLength-27.47, 0.003);
    //ConstantDistanceToTime* linear = new ConstantDistanceToTime(85.76);
    TGraph* g1 = graph(220, -10, 1, d, &timeWrapper);
    g1->SetTitle("Distance->deltaT;Distance (mm);deltaT (ns)");
    TGraph* g2 = graph(220, -10, 1, d, &inverseSpeedWrapper);
    g2->SetTitle("Inverse Speed (ns/mm);Distance (mm); dt/dx (ns/mm)");
    TGraph* g3 = graph(220, -10, 1, d, &speedWrapper);
    g3->SetTitle("Speed (mm/ns); Distance (mm); dx/dt (mm/ns)");
    TGraph* g4 = graph(361, -0.5, 0.00623, d, &distanceWrapper);
    g4->SetTitle("DeltaT->Distance; deltaT (ns); Distance (mm)");

    TFile mefile("CaloDistanceToTime.root","RECREATE");

    TCanvas* can = new TCanvas("CaloDistanceToTime", "CaloDistanceToTime", 1000, 1000);
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
