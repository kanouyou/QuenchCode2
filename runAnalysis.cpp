#include <iostream>
#include <string>
#include <TApplication.h>
#include <TMultiGraph.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TStyle.h>
#include "QAnalysis.h"

using namespace std;

int main(int argc, char** argv) {
  
  TApplication* app = new TApplication("app", &argc, argv);
  
  TCanvas* c0 = new TCanvas("temp", "temp", 1000, 1000);
  c0->Divide(3,3);

  gStyle->SetOptStat(0);
  gStyle->SetNumberContours(99);
  //gStyle->SetTitleSize(0.05);
  gStyle->SetTitleOffset(1.7, "y");
  //gStyle->SetLabelSize(0.045);

  double time = 1.;

  QAnalysis* ana = new QAnalysis("output.root");

  for (int i=0; i<8; i++) {
    c0->cd(i+1);
    gPad->SetTicks(1,1);
	gPad->SetRightMargin(0.13);
	ana->PlotTemperatureDistribution(time*(i+1), 1, false);
  }
  
  // Plot time-temperature
  TCanvas* c1 = new TCanvas("timetemp", "timetemp", 1000, 1000);
  c1->Divide(2,2);
  
  TMultiGraph* mg = new TMultiGraph();
  TGraph* gr[3];
  
  gr[0] = ana->GetTemperatureGraph( 1, 1, 18);
  gr[1] = ana->GetTemperatureGraph(34, 1, 18);
  gr[2] = ana->GetTemperatureGraph( 1, 1,  2);
  
  for (int i=0; i<3; i++) mg->Add(gr[i], "l");
  
  c1->cd(1);
  gPad->SetTicks(1,1);
  mg->SetTitle(";Time [sec];Temperature [K]");
  mg->Draw("a");
  
  // Plot Quenched Cell
  c1->cd(2);
  gPad->SetTicks(1,1);
  
  ana->PlotResistance();
  
  // plot voltage
  c1->cd(3);
  gPad->SetTicks(1,1);
  ana->PlotVoltage();
  
  // plot queched time
  c1->cd(4);
  gPad->SetTicks(1,1);
  ana->PlotVelocityDistribution(1);
  
  // Plot current, resistance, power
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
  
  TCanvas* cv = new TCanvas("power", "power", 650, 700); 
  TPad* pad[3];

  pad[0] = new TPad("current", "current", 0., 0., 1., 1.);
  pad[0]->SetMargin(leftMargin, rightMargin, bottomMargin + 2*figHeight + 2*figVstep, headMargin);
  pad[0]->SetGrid();
  pad[0]->SetTicks(1,1);
  pad[0]->SetAttLinePS(kBlack, 1, 2);
  pad[0]->SetFillColor(0);
  pad[0]->SetFillStyle(0);
  pad[0]->Draw();
  pad[0]->cd(0);

  ana->PlotCurrent();

  pad[1] = new TPad("resist", "resist", 0., 0., 1., 1.);
  pad[1]->SetMargin(leftMargin, rightMargin, bottomMargin + figHeight + figVstep, headMargin + figHeight + figVstep);
  pad[1]->SetGrid();
  pad[1]->SetTicks(1,1);
  pad[1]->SetAttLinePS(kBlack, 1, 2);
  pad[1]->SetFillColor(0);
  pad[1]->SetFillStyle(0);
  pad[1]->Draw();
  pad[1]->cd(0);
  
  ana->PlotResistance();

  pad[2] = new TPad("power", "power", 0., 0., 1., 1.);
  pad[2]->SetMargin(leftMargin, rightMargin, bottomMargin, headMargin + 2*figHeight + 2*figVstep);
  pad[2]->SetGrid();
  pad[2]->SetTicks(1,1);
  pad[2]->SetAttLinePS(kBlack, 1, 2);
  pad[2]->SetFillColor(0);
  pad[2]->SetFillStyle(0);
  pad[2]->Draw();
  pad[2]->cd(0);
  
  ana->PlotPower();
  
  cv->Update();

  // Plot delta temperature
  /*TCanvas* c2 = new TCanvas("delta", "delta", 600, 500);
  c2->SetTicks(1,1);
  c2->SetAttLinePS(kBlack, 1, 2);

  ana->PlotDeltaTemperature();
*/
  app->Run();

  return 0;
}

