#include <iostream>
#include <string>
#include <cmath>
#include "QSuperconduct.h"

using namespace std;

void QSuperconduct::SetStrand(double radius, int no){
  // unit: mm
  double r = radius;
  Asc = M_PI * pow(r, 2) * no * R_nbti / (R_al+R_cu+R_nbti);
}


void QSuperconduct::SetRatio(double R_al, double R_cu, double R_nbti){
  this->R_al = R_al;
  this->R_cu = R_cu;
  this->R_nbti = R_nbti;
}


void QSuperconduct::SetConductor(double width, double thickness) {
  w_cdt = width;
  t_cdt = thickness;
  Asc = width * thickness * R_nbti / (R_al + R_cu + R_nbti);
}


void QSuperconduct::SetParameter(double B, double T) {
  this->B = B;
  this->T = T;
}


void QSuperconduct::SetOperatingCurrent(double Iop) {
  this->Iop = Iop;
}


void QSuperconduct::SetTemperature(double temp) {
  T = temp;
}


void QSuperconduct::SetField(double field) {
  B = field;
}


double QSuperconduct::GetCriticalCurrent(){
  double Jc = CriticalCurrent(B, T) * Asc;
  return Jc;
}


double QSuperconduct::GetCriticalTemperature() {
  double Tc = CriticalTemperature(B);
  return Tc;
}


double QSuperconduct::GetCurrentSharingTemp(double field, double iop){
  double T0 = 4.2;
  double Tc = CriticalTemperature(field);
  double Ic = CriticalCurrent(field, T0) * Asc;

  Tcs = T0 + (Tc - T0)*(1 - iop/Ic);
  return Tcs;
}


/*double QSuperconduct::GetOperatingCurrent(double tcs) {
  double Iopt = CriticalCurrent(tcs) * Asc;
  return Iopt;
}*/


double QSuperconduct::CriticalCurrent(double field, double temp) {
// L. Bottura's fitting equation
// LHC project report 358
// Jc(B,T): critical current density
// Tc(B): critical temperature
// Bc2(T): upper critical field
// Tc0(B=0): maxmium critical temperature
// Bc20(T=0): maxmium critical field
// Tcs(B,Jop): current sharing temperature
  const int m = 5;
  float t, b;
  float Bc2, jc;
  float n = 1.7;
  float Tc0[m]   = { 9.2,  8.5,  8.9,  9.2,  9.35};
  float Bc20[m]  = {14.5, 14.2, 14.4, 14.4, 14.25};
  float C0[m]    = {23.8, 28.6, 28.5, 37.7, 28.40};
  float alpha[m] = {0.57, 0.76, 0.64, 0.89,  0.80};
  float beta[m]  = {0.90, 0.85, 0.75, 1.10,  0.89};
  float gamma[m] = {1.90, 1.76, 2.30, 2.09,  1.87};

  int data = 0;
  
  t = temp / Tc0[data];
  Bc2 = Bc20[data] * (1 - pow(t, n));
  b = field / Bc2;
  
  // J0: J(T=4.2K, B=5T)
  float J0 = 3000.0;   // A/mm2
  jc = C0[data] * pow(b, alpha[data]) * pow((1-b), beta[data]) * pow((1-pow(t, n)), gamma[data]) / field;
  
  jc = jc * J0;
  return jc;
}


double QSuperconduct::CriticalTemperature(double field) {
  double tc;
  const int m = 5;
  float Tc0[m]   = { 9.2,  8.5,  8.9,  9.2,  9.35};
  float Bc20[m]  = {14.5, 14.2, 14.4, 14.4, 14.25};
 
  int data = 0;
  float n = 1.7;
  tc = Tc0[data] * pow((1-field/Bc20[data]), 1/n);
  return tc;
}


double QSuperconduct::GetRadiationEffect(string name, int time, double flux) {
  double timeFactor = time / 365.;
  double neuFactor  = 2.7e-22;
  double rho = 0.;
  double RRR = 0.;

  if (name=="Al" || name=="Aluminium") {
    rho = 0.0675;
    RRR = 27. / (rho + flux * timeFactor * neuFactor);
  }
  else if (name=="Cdt" || name=="Conductor") { 
    rho = 0.0135;
    RRR = 27. / (rho + flux * timeFactor * neuFactor);
  }
  else
    cout << "Error: no such kind of material!" << endl;
  
  return RRR;
}


