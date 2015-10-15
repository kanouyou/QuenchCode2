#include <iostream>
#include <TApplication.h>
#include "QMakePlot.h"

using std::string;

int main (int argc, char** argv) {
  
  TApplication* app = new TApplication("app", &argc, argv);
  
  string filename = "output.root";
  int MeshR = 19;
  int MeshP = 4;
  int MeshZ = 90;

  QMakePlot* plot = new QMakePlot(filename, MeshZ, MeshP, MeshR);
  //plot->ShowTemperatureDistribution();
  //plot->ShowQuenchVelocity();
  //plot->ShowTemperature();
  plot->ShowElectric();

  app->Run();
  return 0;
}
