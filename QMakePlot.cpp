#include <string>
#include <TPad.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TLegend.h>
#include <TStyle.h>
#include "QMakePlot.h"

using std::string;

QMakePlot::QMakePlot (string filename, int MeshZ, int MeshP, int MeshR) : QAnalysis(filename) {
  zbin = MeshZ;
  pbin = MeshP;
  rbin = MeshR;
  spot = new int[3];
  spot[0] = 1;
  spot[1] = 1;
  spot[2] = 18;
}

QMakePlot::~QMakePlot() {
}

bool QMakePlot::ShowTemperatureDistribution () {
  bool   scale = false;
  double t0    = 1.;       // [sec]
  double tf    = 21.;      // [sec]
  int    nt    = 20;
  double dt    = (tf-t0)/nt;
  
  TCanvas* c0 = new TCanvas("tmpdis", "tmpdis", 1500, 1000);
  c0->Divide(4,5);
  
  gStyle->SetOptStat(0);

  for (int i=0; i<nt; i++) {
    c0->cd(i+1);
	gPad->SetTicks(1,1);
	gPad->SetRightMargin(0.14);
	PlotTemperatureDistribution(i*dt+t0, 1, scale);
  }
  c0->Print("tempdistribution.pdf");
  
  return true;
}

bool QMakePlot::ShowQuenchVelocity () {
  int nphi = 4;
  
  TCanvas* c0 = new TCanvas("velocity", "velocity", 700, 500);
  c0->Divide(2,2);

  for (int i=0; i<nphi; i++) {
    c0->cd(i+1);
    gPad->SetTicks(1,1);
	gPad->SetRightMargin(0.14);
	PlotQuenchVelocity(i+1);
  }
  c0->Print("velocity.pdf");

  return true;
}

bool QMakePlot::ShowTemperature() {
  TMultiGraph* mg = new TMultiGraph();
  TGraph* gr[5];
  for (int i=0; i<5; i++) gr[i] = new TGraph();
  int color[5] = {kRed, kAzure+3, kOrange+7, kGreen+3, kBlack};

  gr[0] = (TGraph*)GetTemperatureGraph(spot[0], spot[1], spot[2]);
  gr[1] = (TGraph*)GetTemperatureGraph(1, 1, 2);
  gr[2] = (TGraph*)GetTemperatureGraph(1, 1, 18);
  gr[3] = (TGraph*)GetTemperatureGraph(90, 1, 2);
  gr[4] = (TGraph*)GetTemperatureGraph(90, 1, 18);
  //gr[5] = (TGraph*)GetTemperatureGraph(90, 2, 2);
  
  for (int i=0; i<5; i++) {
    gr[i]->SetLineStyle(i+1);
	gr[i]->SetLineColor(color[i]);
	gr[i]->SetLineWidth(3);
	mg->Add(gr[i], "l");
  }
  
  TCanvas* c0 = new TCanvas("temperature", "temperature", 650, 500);
  c0->SetTicks(1,1);
  c0->SetLeftMargin(0.11);
  c0->SetBottomMargin(0.11);
  
  gStyle->SetTitleSize(0.05, "xy");
  gStyle->SetLabelSize(0.05, "xy");

  mg->SetTitle("; Time [sec]; Temperature [K]");
  mg->Draw("a");
  
  TLegend* lg = new TLegend(0.12, 0.68, 0.4, 0.89);
  lg->SetTextSize(0.045);
  lg->SetBorderSize(0);
  lg->SetFillStyle(0);
  lg->AddEntry(gr[0], "Quench Spot", "l");
  lg->AddEntry(gr[1], "(1, 1, 2)", "l");
  lg->AddEntry(gr[2], "(1, 1, 18)", "l");
  lg->AddEntry(gr[3], "(90, 1, 2)", "l");
  lg->AddEntry(gr[4], "(90, 1, 18)", "l");
  lg->Draw();

  c0->Print("temperature.pdf");

  return true;
}

bool QMakePlot::ShowElectric() {
  // ************************************** //
  // set margin
  // ************************************** //
  double leftMargin   = 0.15;
  double rightMargin  = 0.05;
  double headMargin   = 0.06;
  double bottomMargin = 0.10;
  double figHeight    = 0.27;
  double figVstep     = 0.01;
  // ************************************** //
  
  TCanvas* c0 = new TCanvas("power", "power", 650, 700);
  gStyle->SetTitleOffset(1.8, "y");
  TPad* pad[3];

  pad[0] = new TPad("current", "current", 0., 0., 1., 1.);
  pad[0]->SetMargin(leftMargin, rightMargin, bottomMargin + 2*figHeight + 2*figVstep, headMargin);
  //pad[0]->SetGrid();
  pad[0]->SetTicks(1,1);
  //pad[0]->SetAttLinePS(kBlack, 1, 2);
  pad[0]->SetFillColor(0);
  pad[0]->SetFillStyle(0);
  pad[0]->Draw();
  pad[0]->cd(0);
  PlotCurrent();

  pad[1] = new TPad("resist", "resist", 0., 0., 1., 1.);
  pad[1]->SetMargin(leftMargin, rightMargin, bottomMargin + figHeight + figVstep, headMargin + figHeight + figVstep);
  //pad[1]->SetGrid();
  pad[1]->SetTicks(1,1);
  //pad[1]->SetAttLinePS(kBlack, 1, 2);
  pad[1]->SetFillColor(0);
  pad[1]->SetFillStyle(0);
  pad[1]->Draw();
  pad[1]->cd(0);
  PlotVoltage();

  pad[2] = new TPad("power", "power", 0., 0., 1., 1.);
  pad[2]->SetMargin(leftMargin, rightMargin, bottomMargin, headMargin + 2*figHeight + 2*figVstep);
  //pad[2]->SetGrid();
  pad[2]->SetTicks(1,1);
  //pad[2]->SetAttLinePS(kBlack, 1, 2);
  pad[2]->SetFillColor(0);
  pad[2]->SetFillStyle(0);
  pad[2]->Draw();
  pad[2]->cd(0);
  PlotPower();
  
  c0->Update();
  c0->Print("electric.pdf");

  return true;
}

bool QMakePlot::ShowDeltaTemp() {
  TCanvas* c0 = new TCanvas("deltaT", "deltaT", 600, 500);
  c0->SetTicks(1,1);

  
  return false;
}

