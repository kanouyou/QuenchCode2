#include <iostream>
#include <TH2F.h>
#include <TGraph.h>
#include "QAnalysis.h"

using namespace std;

QAnalysis::QAnalysis(string fileName) {
  file = new TFile(fileName.c_str());
  zmax = 90.;
  zbin = 90;
}

TTree* QAnalysis::ReadTree(TFile* file) {
  TTree* tree = (TTree*)file->Get("tree");
  
  posID = new int[3];
  fpos  = new double[3];
  tree->SetBranchAddress(         "id",    posID);
  tree->SetBranchAddress(   "position",     fpos);
  tree->SetBranchAddress(       "time",    &time);
  tree->SetBranchAddress(    "current", &current);
  tree->SetBranchAddress(      "field",   &field);
  tree->SetBranchAddress( "resistance",       &R);
  tree->SetBranchAddress( "quenchcell", &qchCell);
  tree->SetBranchAddress("temperature",    &temp);
  tree->SetBranchAddress( "quenchTime", &qchTime);

  return tree;
}

void QAnalysis::PlotTemperatureDistribution(double qtime, int numPhi, bool scale) {
  int    rbin = 19;
  double rmax = 19.;
  TH2F*  hist = new TH2F(Form("tempDis%.1f", qtime), Form("tempDis%.1f", qtime), zbin, 0, zmax, rbin, 0, rmax);
  TTree* tree = ReadTree(file);
  
  for (int i=0; i<tree->GetEntries(); i++) {
    tree->GetEntry(i);
	if ( time==qtime && posID[1]==numPhi ) 
      hist->Fill(posID[0], posID[2], temp);
  }
  
  if (scale==true) hist->GetZaxis()->SetRangeUser(0., GetMaximum("temp"));
  hist->SetTitle(Form("time = %.1f [sec]; Z; R; Temperature [K]", qtime));
  hist->GetXaxis()->SetTitleSize(0.05);
  hist->GetYaxis()->SetTitleSize(0.05);
  hist->GetZaxis()->SetTitleSize(0.05);
  hist->GetXaxis()->SetLabelSize(0.05);
  hist->GetYaxis()->SetLabelSize(0.05);
  hist->GetZaxis()->SetLabelSize(0.05);
  hist->Draw("cont4z");
}

void QAnalysis::PlotSpotTemperature(int z, int phi, int r) {
  TTree* tree = ReadTree(file);
  TGraph* gr  = new TGraph();
  
  int cnt = 0;
  for (int i=0; i<tree->GetEntries(); i++) {
    tree->GetEntry(i);
	if ( posID[0]==z && posID[1]==phi && posID[2]==r ) {
	  gr->SetPoint(cnt, time, temp);
	  cnt += 1;
    }
  }

  gr->SetTitle(Form("Position: [%i, %i, %i]; Time [sec]; Temperature [K]", z, phi, r));
  gr->SetLineColor(kRed);
  gr->SetLineWidth(2);
  gr->Draw("al");

  cout << "Finished PlotSpotTemperature()" << endl;
}

TGraph* QAnalysis::GetTemperatureGraph(int z, int phi, int r) {
  TTree* tree = ReadTree(file);
  TGraph* gr  = new TGraph();
  
  int cnt = 0;
  for (int i=0; i<tree->GetEntries(); i++) {
    tree->GetEntry(i);
	if ( posID[0]==z && posID[1]==phi && posID[2]==r ) {
	  gr->SetPoint(cnt, time, temp);
	  cnt += 1;
    }
  }

  gr->SetLineColor(kRed);
  gr->SetLineWidth(2);
  return gr;
}

void QAnalysis::PlotCurrent() {
  TTree* tree = ReadTree(file);
  TGraph* gr  = new TGraph();
  
  int z = 1; int phi = 1; int r = 18;
  int cnt = 0;
  for (int i=0; i<tree->GetEntries(); i++) {
    tree->GetEntry(i);
	if ( posID[0]==z && posID[1]==phi && posID[2]==r ) {
      gr->SetPoint(cnt, time, current);
	  cnt += 1;
	}
  }
  
  gr->GetXaxis()->SetLabelSize(0);
  gr->GetXaxis()->SetTickLength(0.01);
  gr->SetTitle("; ; Current [A]");
  gr->SetLineColor(kAzure+3);
  gr->SetLineWidth(2);
  gr->Draw("al"); 
}

double QAnalysis::GetMaximum(string name) {
  TTree* tree = ReadTree(file);
  double max  = tree->GetMaximum(name.c_str());
  return max;
}

double QAnalysis::GetMinimum(string name) {
  TTree* tree = ReadTree(file);
  double min  = tree->GetMinimum(name.c_str());
  return min;
}

void QAnalysis::PlotQuenchedNumber() {
  TTree* tree = ReadTree(file);
  TGraph* gr  = new TGraph();
  
  int meshNum = zbin * 19 * 4;
  int cnt = 0;
  for (int i=0; i<tree->GetEntries(); i++) {
    tree->GetEntry(i);
	if ( posID[0]==1 && posID[1]==1 && posID[2]==18 ) {
      gr->SetPoint(cnt, time, static_cast<double>(qchCell)/meshNum);
	  cnt += 1;
	}
  }

  gr->SetTitle("Quenched Cell; Time [sec]; Quenched Ratio [%]");
  gr->SetLineWidth(2);
  gr->SetLineColor(kBlue+1);
  gr->Draw("al");
}

void QAnalysis::PlotResistance() {
  TTree* tree = ReadTree(file);
  TGraph* gr  = new TGraph();
  
  int cnt = 0;
  for (int i=0; i<tree->GetEntries(); i++) {
    tree->GetEntry(i);
	if ( posID[0]==1 && posID[1]==1 && posID[2]==18 ) {
      gr->SetPoint(cnt, time, R);
	  cnt += 1;
	}
  }
  
  gr->GetXaxis()->SetTickLength(0.01);
  gr->GetXaxis()->SetLabelSize(0);
  gr->SetTitle("; ; Resistance [#Omega]");
  gr->SetLineWidth(2);
  gr->SetLineColor(kRed);
  gr->Draw("al");
}

void QAnalysis::PlotVoltage() {
  TTree* tree = ReadTree(file);
  TGraph* gr  = new TGraph();
  
  int cnt = 0;
  for (int i=0; i<tree->GetEntries(); i++) {
    tree->GetEntry(i);
	if ( posID[0]==1 && posID[1]==1 && posID[2]==18 ) {
      gr->SetPoint(cnt, time, R * current);
	  cnt += 1;
	}
  }
  
  gr->GetXaxis()->SetTickLength(0.01);
  gr->GetXaxis()->SetLabelSize(0);
  gr->SetTitle("; Time [sec]; Voltage [V]");
  gr->SetLineWidth(2);
  gr->SetLineColor(kPink-1);
  gr->Draw("al");
}

void QAnalysis::PlotPower() {
  TTree* tree = (TTree*)ReadTree(file);
  TGraph* gr  = new TGraph();

  int cnt = 0;
  for (int i=0; i<tree->GetEntries(); i++) {
    tree->GetEntry(i);
	if (posID[0]==1 && posID[1]==1 && posID[2]==18) {
      gr->SetPoint(cnt, time, R*current*current*1e-3);
	  cnt++;
	}
  }
  
  gr->GetXaxis()->SetTickLength(0.01);
  gr->SetTitle("; Time [sec]; Power [kW]");
  gr->SetLineWidth(2);
  gr->SetLineColor(kTeal+3);
  gr->Draw("al");
}

void QAnalysis::PlotVelocityDistribution(int phi) {
  TTree* tree = ReadTree(file);
  
  int rbin = 19;
  TH2F*  hist = new TH2F("quenchv", "quenchv", zbin, 0, zmax, rbin, 0, 19);
  
  for (int i=0; i<tree->GetEntries(); i++) {
    tree->GetEntry(i);
    if ( time==70. && posID[1]==phi )
	  hist->Fill(posID[0], posID[2], qchTime);
  }
  
  hist->SetTitle("; Z; R; Quenched Time [sec]");
  hist->Draw("colz");
}

double QAnalysis::GetMaxTemperature(double fTime) {
  TTree* tree = (TTree*)ReadTree(file);
  vector<double> fTemperature;

  for (int i=0; i<tree->GetEntries(); i++) {
    tree->GetEntry(i);
    if (time==fTime) 
      fTemperature.push_back(temp);
  }
  
  double fMaximum = fTemperature[0];

  for (vector<int>::size_type i=0; i<fTemperature.size(); i++) {
    if (fMaximum<fTemperature[i])
	  fMaximum = fTemperature[i];
  }

  return fMaximum;
}

double QAnalysis::GetMinTemperature(double fTime) {
  TTree* tree = (TTree*)ReadTree(file);
  vector<double> fTemperature;
  
  for (int i=0; i<tree->GetEntries(); i++) {
    tree->GetEntry(i);
    if (time==fTime) 
      fTemperature.push_back(temp);
  }
  
  double fMinimum = fTemperature[0];

  for (vector<int>::size_type i=0; i<fTemperature.size(); i++) {
    if (fMinimum>fTemperature[i])
	  fMinimum = fTemperature[i];
  }

  return fMinimum;
}

void QAnalysis::PlotDeltaTemperature() {
  TGraph* gr = new TGraph();
  int nt = 80;
  double fMin, fMax;
  
  for (int i=0; i<nt; i++) {
    fMin = GetMinTemperature(static_cast<double>(i+1));
	fMax = GetMaxTemperature(static_cast<double>(i+1));
    gr->SetPoint(i, static_cast<double>(i+1), fMax-fMin);
	cout << "Time: " << i << " [sec]; Maximum: " << fMax << " [K]; Minimum: " << fMin << " [K]" << endl;
  }

  gr->SetLineStyle(2);
  gr->SetLineWidth(2);
  gr->SetLineColor(kOrange);
  gr->SetMarkerStyle(20);
  gr->SetMarkerColor(kOrange);
  gr->SetTitle("; Time [sec]; #Delta T [K]");
  gr->Draw("al");
}

double* QAnalysis::GetQuenchPoint(int z, int phi, int r) {
  double* position = new double[3];
  TTree* tree = ReadTree(file);

  for (int i=0; i<tree->GetEntries(); i++) {
    tree->GetEntry(i);
	if (time==1. && posID[0]==z && posID[1]==phi && posID[2]==r) {
      position[0] = fpos[0];
	  position[1] = fpos[1];
	  position[2] = fpos[2];
	  break;
	}
  }

  return position;
}

void QAnalysis::PlotQuenchVelocity(int phi) {
  TTree* tree = ReadTree(file);
  
  int quenchedZ = 1;
  int quenchedP = 1;
  int quenchedR = 18;
  double* qPos = new double[3];
  qPos = GetQuenchPoint(quenchedZ, quenchedP, quenchedR);

  double cPos;
  int rbin = 9;
  TH2F* hist = new TH2F(Form("QuenchVelocity%i",phi), Form("QuenchVelocity%i",phi), zbin, 0, zmax, rbin, 0, 9);
  
  for (int i=0; i<tree->GetEntries(); i++) {
    tree->GetEntry(i);
    if ( time==70. && posID[1]==phi && posID[2]%2==0) {
	  cPos = sqrt(pow(abs(fpos[0]-qPos[0]),2) + pow(abs(fpos[1]-qPos[1]),2) + pow(abs(fpos[2]-qPos[2]),2));
	  hist->Fill(posID[0], posID[2]/2, cPos / qchTime);
    }
  }
  
  hist->GetXaxis()->SetTitleSize(0.05);
  hist->GetYaxis()->SetTitleSize(0.05);
  hist->GetZaxis()->SetTitleSize(0.05);
  hist->GetXaxis()->SetLabelSize(0.05);
  hist->GetYaxis()->SetLabelSize(0.05);
  hist->GetZaxis()->SetLabelSize(0.05);
  hist->SetTitle(Form("No. #phi = %i; Z; R; Quenched Velocity [m/sec]",phi));
  hist->Draw("cont4z");
}

