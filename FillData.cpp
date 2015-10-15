#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <TH2F.h>
#include <TString.h>
#include "FillData.hh"

using namespace std;

FillData::FillData(string fileName) {
  this->fileName = fileName;
  file = new TFile(fileName.c_str(), "recreate");
  tree = new TTree("tree", "Data File filled with FDM simulation result");
  
  id  = new int[3];
  fpos = new double[3];
  tree->Branch(         "id",    id,       "id[3]/I");
  tree->Branch(   "position",  fpos, "position[3]/D");
  tree->Branch(       "time",    &t,        "time/D");
  tree->Branch(    "current",    &I,     "current/D");
  tree->Branch(      "field",    &B,       "field/D");
  tree->Branch( "resistance",    &R,  "resistance/D");
  tree->Branch( "resist", &res, "resist/D");
  tree->Branch( "quenchcell",    &n,  "quenchcell/I");
  tree->Branch("temperature",    &T,        "temp/D");
  tree->Branch( "quenchTime", &qchT,  "quenchTime/D");
}


void FillData::FillTree(double time, int i, int j, int k, double x, double y, double z, 
                        double cI, double cB, double resist, double cR, int cn, double temperature, double quench) {
  t     = time;
  id[0] = i;
  id[1] = j;
  id[2] = k;
  fpos[0] = x;
  fpos[1] = y;
  fpos[2] = z;
  I = cI;
  B = cB;
  res = resist;
  R = cR;
  n = cn;
  T = temperature;
  qchT = quench;
  tree->Fill();
}

void FillData::Close() {
  file->Write();
  file->Close();
}


void FillData::Plot(string fName) {
  TFile* file = new TFile(fName.c_str());
  TTree* tree = (TTree*)file->Get("tree");
  
  double temp, time;
  double pos[3];
  int    vec[3];
  
  tree->SetBranchAddress("id", vec);
  tree->SetBranchAddress( "position",   pos);
  tree->SetBranchAddress("temperature", &temp);
  tree->SetBranchAddress("time", &time);
  
  double zmin = 0.;
  double zmax = 340.;
  int    zbin = 68;
  double rmin = 0.;
  double rmax = 25.;
  int    rbin = 8;

  TH2F* hist = new TH2F("temperature", "temp", zbin, zmin, zmax,
  											   rbin, rmin, rmax);
  for (int i=0; i<tree->GetEntries(); i++) {
    tree->GetEntry(i);
	if (vec[1]==1 && time==99.2)
	  hist->Fill(pos[0], pos[2], temp);
  }
  hist->Draw("colz");
}


double FillData::GetTheta(double xc, double yc) {
  double Theta;
  if (xc>=0)
    Theta = acos(yc/sqrt(pow(xc, 2) + pow(yc, 2)));
  else
    Theta = M_PI + acos(yc/sqrt(pow(xc, 2) + pow(yc, 2)));
  
  Theta = Theta * 180 / M_PI;
  return Theta;
}
