#ifndef QuenchMain_HH
#define QuenchMain_HH

#include <string>
#include <TCanvas.h>
#include <TGraph2D.h>
#include <fstream>
#include "QSuperconduct.h"
#include "MParameter.hh"
#include "QMaterial.h"

// meshing
namespace Mesh {
  const int    Mz   = 90 + 1;
  const int    Mphi = 4 + 1;
  const int    Mr   = 19 + 1;
  // time mesh:
  const double t0   = 0.;
  const double tf   = 80.;
  const int    Mt   = 800000;
  const double dt   = (tf - t0) / Mt;
  const int    pt   = 500;
}

// unit: [m]
namespace Kapton {
  const double l_Resin   = 0.5 * 1e-3;
  const double k_Resin   = 0.05;
  const double thinkness = 0.135 * 1e-3;
  const double lr = thinkness * 2 + l_Resin;
  const double lz = thinkness * 2;
  QMaterial* mat = new QMaterial("Kapton");
}

namespace Aluminium {
  const double density = 2700.;
  const double strip   = 1. * 1e-3;
  const double inner   = 3. * 1e-3;
  QMaterial* mat = new QMaterial("Al");
}

namespace Shell {
  const double lr  = 50. * 1e-3;
  const double lHe = 50. * 1e-3;
  const double RRR = 40.;
}

namespace Solenoid {
  const int    turn     = 270;
  const int    layer    = 9;
  const double current  = 2700.;    // [A]
  const double resistor = 0.185;    // [Ohm]
  const double induct   = 12.69;    // [H]
  const double field    = 5.1;      // [T]
  const int    operTime = 90;       // [days]
}

namespace Conductor {
  const double density = 4000.;
  const double lz      = 4.73 * 1e-3;
  const double lr      = 15. * 1e-3;
  const double rAl     = 7.3 / 9.3;
  const double rNbTi   = 1. / 9.3;
  const double rCu     = 1. / 9.3;
  const double rArea   = static_cast<double>(Solenoid::turn) / (Mesh::Mz - 1);
  const double aAl     = lz * lr * rAl;                                              // cross section of Al stabilizer
  const double areaAl  = lz * lr * rAl * rArea;                                      // cross section of cell of Al
  const double cdtvol  = lz * lr * (2*M_PI*0.8) / (Mesh::Mphi-1);                    // volume of Al stabilizer
  const double volume  = rArea * lz * lr * (2*M_PI*0.8) / (Mesh::Mphi-1);            // volume of cell of Al
  const double current = Solenoid::current;
  const double RRR_Cu   = 50.;
  QSuperconduct* sp = new QSuperconduct();
}

namespace Dimension {
  const double lz = Solenoid::turn  * (Conductor::lz + Kapton::lz * 2);
  const double lr = Solenoid::layer * (Conductor::lr + Kapton::lr * 2) + (Solenoid::layer-1) * Aluminium::strip
                    + Aluminium::inner;
  const double lp = 2 * M_PI * 0.80;
  const double dz = lz / (Mesh::Mz   - 1);
  const double dr = lr / (Mesh::Mr   - 1);
  const double dp = lp / (Mesh::Mphi - 1);
}

namespace Field {
  const double factor = 1.;
  const int    n    = 5;
  //const double B[n] = {5.4, 4.0, 2.9, 1.1, 0.2};
  const double B[n] = {0.2, 1.1, 2.9, 4.0, 5.4};
}

namespace Radiation {
  const int nz = 5;
  const int np = 4;
  const int nr = 9;
  // PHITS data (flux[n/m2/sec]):
  double neutron[5][4][9] = {{{2.7e+20, 4.0e+20, 5.3e+20, 6.8e+20, 8.5e+20, 1.0e+21, 1.3e+21, 1.5e+21, 1.9e+21},
	   						  {2.8e+20, 4.1e+20, 5.4e+20, 6.9e+20, 8.4e+20, 1.0e+21, 1.3e+21, 1.5e+21, 1.8e+21},
		   					  {2.6e+20, 3.8e+20, 5.0e+20, 6.3e+20, 7.9e+20, 9.6e+20, 1.2e+21, 1.4e+21, 1.7e+21},
		   					  {2.7e+20, 4.0e+20, 5.3e+20, 6.7e+20, 8.3e+20, 1.0e+21, 1.2e+21, 1.5e+21, 1.8e+21}},
		   					 {{5.0e+20, 7.4e+20, 9.8e+20, 1.3e+21, 1.5e+21, 1.9e+21, 2.3e+21, 2.8e+21, 3.3e+21},
		   					  {5.0e+20, 7.4e+20, 9.7e+20, 1.2e+21, 1.5e+21, 1.9e+21, 2.2e+21, 2.7e+21, 3.2e+21},
		   					  {4.4e+20, 6.4e+20, 8.5e+20, 1.1e+21, 1.3e+21, 1.6e+21, 1.9e+21, 2.3e+21, 2.8e+21},
		   					  {4.9e+20, 7.2e+20, 9.8e+20, 1.2e+21, 1.5e+21, 1.9e+21, 2.2e+21, 2.7e+21, 3.2e+21}},
		   					 {{6.3e+20, 9.4e+20, 1.2e+21, 1.6e+21, 2.0e+21, 2.4e+21, 2.8e+21, 3.4e+21, 4.0e+21},
		   					  {6.6e+20, 9.8e+20, 1.3e+21, 1.6e+21, 2.0e+21, 2.5e+21, 3.0e+21, 3.6e+21, 4.2e+21},
		   					  {6.7e+20, 9.8e+20, 1.3e+21, 1.6e+21, 2.0e+21, 2.4e+21, 2.9e+21, 3.5e+21, 4.2e+21},
		   					  {6.7e+20, 9.9e+20, 1.3e+21, 1.7e+21, 2.0e+21, 2.5e+21, 3.0e+21, 3.6e+21, 4.3e+21}},
		   					 {{5.6e+20, 8.2e+20, 1.1e+21, 1.4e+21, 1.7e+21, 2.0e+21, 2.5e+21, 3.0e+21, 3.5e+21},
		   					  {6.5e+20, 9.4e+20, 1.2e+21, 1.5e+21, 1.9e+21, 2.3e+21, 2.8e+21, 3.3e+21, 3.9e+21},
		   					  {7.4e+20, 1.1e+21, 1.4e+21, 1.8e+21, 2.2e+21, 2.6e+21, 3.2e+21, 3.8e+21, 4.5e+21},
		   					  {6.5e+20, 9.5e+20, 1.2e+21, 1.6e+21, 1.9e+21, 2.3e+21, 2.8e+21, 3.4e+21, 4.0e+21}},
		   					 {{4.3e+20, 6.2e+20, 8.0e+20, 1.0e+21, 1.3e+21, 1.5e+21, 1.8e+21, 2.2e+21, 2.7e+21},
		   					  {4.7e+20, 6.8e+20, 8.9e+20, 1.1e+21, 1.4e+21, 1.7e+21, 2.0e+21, 2.4e+21, 2.9e+21},
		   					  {6.3e+20, 9.0e+20, 1.2e+21, 1.5e+21, 1.8e+21, 2.2e+21, 2.6e+21, 3.2e+21, 3.8e+21},
		   					  {4.9e+20, 7.1e+20, 9.4e+20, 1.2e+21, 1.4e+21, 1.7e+21, 2.1e+21, 2.5e+21, 3.0e+21}}};
  double Heat[nz][np][nr] = {{{0.0039839,0.0052346,0.0066784,0.0079471,0.0087378,0.010078 ,0.010896,0.011942,0.013295},
							 {0.0042479,0.0053834,0.0064554,0.0078091,0.0088091,0.010061 ,0.011338,0.012771,0.013936},
						     {0.0039344,0.0050329,0.006656 ,0.0073079,0.0086217,0.0099013,0.011124,0.012205,0.013177},
						     {0.0042332,0.0052832,0.0068928,0.007934 ,0.0089066,0.010333 ,0.011082,0.012561,0.013923}},
						    {{0.0078342,0.010151 ,0.012021 ,0.013344 ,0.015653 ,0.018076 ,0.020017,0.022166,0.023062},
						     {0.0083035,0.0098356,0.011727 ,0.013569 ,0.016335 ,0.017761 ,0.02019 ,0.021779,0.022267},
						     {0.007858 ,0.0097305,0.011221 ,0.012861 ,0.013835 ,0.016509 ,0.018195,0.019405,0.020307},
						     {0.0071509,0.009249 ,0.011935 ,0.014357 ,0.015637 ,0.017828 ,0.019267,0.021666,0.022812}},
						    {{0.009906 ,0.012878 ,0.015556 ,0.018677 ,0.020531 ,0.022587 ,0.025725,0.027969,0.030156},
						     {0.0108   ,0.013215 ,0.015806 ,0.019539 ,0.022341 ,0.023553 ,0.026644,0.02862 ,0.032423},
						     {0.013019 ,0.014827 ,0.016966 ,0.019033 ,0.021219 ,0.024848 ,0.026114,0.029677,0.032087},
						     {0.011332 ,0.013547 ,0.016293 ,0.019425 ,0.022645 ,0.025015 ,0.026941,0.029077,0.031504}},
						    {{0.009647 ,0.011685 ,0.014956 ,0.016235 ,0.0199   ,0.021879 ,0.023387,0.025336,0.025921},
						     {0.010918 ,0.014488 ,0.017292 ,0.019154 ,0.021663 ,0.02413  ,0.026612,0.027808,0.030146},
						     {0.013048 ,0.016393 ,0.01887  ,0.022813 ,0.024618 ,0.028264 ,0.02946 ,0.032857,0.035303},
						     {0.011586 ,0.014782 ,0.018667 ,0.020779 ,0.022632 ,0.023467 ,0.02553 ,0.028957,0.030982}},
						    {{0.0081708,0.0093208,0.010351 ,0.011625 ,0.013956 ,0.014655 ,0.016447,0.017737,0.019267},
						     {0.0092167,0.010401 ,0.012    ,0.013363 ,0.015914 ,0.017457 ,0.018349,0.020093,0.020839},
						     {0.012266 ,0.015218 ,0.017407 ,0.018611 ,0.021234 ,0.023326 ,0.02317 ,0.027249,0.0281  },
						     {0.0092928,0.010844 ,0.012343 ,0.01396  ,0.015126 ,0.017406 ,0.019553,0.022371,0.023277}}};
}

static double z      [Mesh::Mz+1][Mesh::Mphi+1][Mesh::Mr+1];
static double phi    [Mesh::Mz+1][Mesh::Mphi+1][Mesh::Mr+1];
static double r      [Mesh::Mz+1][Mesh::Mphi+1][Mesh::Mr+1];
static double T      [Mesh::Mz+1][Mesh::Mphi+1][Mesh::Mr+1];
static double preT   [Mesh::Mz+1][Mesh::Mphi+1][Mesh::Mr+1];
static double RRR_Al [Mesh::Mz+1][Mesh::Mphi+1][Mesh::Mr+1];
static double RRR_Cdt[Mesh::Mz+1][Mesh::Mphi+1][Mesh::Mr+1];
static double kr     [Mesh::Mz+1][Mesh::Mphi+1][Mesh::Mr+1];
static double kz     [Mesh::Mz+1][Mesh::Mphi+1][Mesh::Mr+1];
static double kp     [Mesh::Mz+1][Mesh::Mphi+1][Mesh::Mr+1];
static double C      [Mesh::Mz+1][Mesh::Mphi+1][Mesh::Mr+1];
static double Heat   [Mesh::Mz+1][Mesh::Mphi+1][Mesh::Mr+1];
static double dz     [Mesh::Mz+1][Mesh::Mphi+1][Mesh::Mr+1];
static double dzz    [Mesh::Mz+1][Mesh::Mphi+1][Mesh::Mr+1];
static double dr     [Mesh::Mz+1][Mesh::Mphi+1][Mesh::Mr+1];
static double drr    [Mesh::Mz+1][Mesh::Mphi+1][Mesh::Mr+1];
static double dp     [Mesh::Mz+1][Mesh::Mphi+1][Mesh::Mr+1];
static double dpp    [Mesh::Mz+1][Mesh::Mphi+1][Mesh::Mr+1];
static double rho    [Mesh::Mz+1][Mesh::Mphi+1][Mesh::Mr+1];
static double resist [Mesh::Mz+1][Mesh::Mphi+1][Mesh::Mr+1];
static double qchTime[Mesh::Mz+1][Mesh::Mphi+1][Mesh::Mr+1];
static double csT    [Mesh::Mz+1][Mesh::Mphi+1][Mesh::Mr+1];
static double fieldB [Mesh::Mz+1][Mesh::Mphi+1][Mesh::Mr+1];
static bool   trigger[Mesh::Mz+1][Mesh::Mphi+1][Mesh::Mr+1];
static bool   cstrg  [Mesh::Mz+1][Mesh::Mphi+1][Mesh::Mr+1];

static double totRes  = 0.;
static int    qchCell = 0;
static double Time;
static double field   = Solenoid::field;
static double current = Solenoid::current;

double RTok(double rho){
	double T   = 4.5;
	double Lwf = 2.44e-8;

	// Resistivity of conduct at room temperature
	double k;
	k = Lwf * T / rho;
	return k;
}

double CriticalTemperature(double B){
	double T;
	T = 9.35*pow(1-B/14.25,1/1.7);
	return T;
}

void FillMagneticField() {
  for (int i=1; i<Mesh::Mz; i++) {
    for (int j=1; j<Mesh::Mphi; j++) {
      for (int k=1; k<Mesh::Mr; k++) {
        fieldB[i][j][k] = Field::factor * Field::B[k/(Field::n-1)] * exp(-(Solenoid::resistor + totRes) * Time / Solenoid::induct);
	    //fieldB[i][j][k] = Solenoid::field;
	  }
	}
  }
}

void GetConductorResistance(double I) {
  double res = 0.;
  int    n   = 0;
  double cs;

  for (int i=1; i<Mesh::Mz; i++) {
    for (int j=1; j<Mesh::Mphi; j++) {
      for (int k=1; k<Mesh::Mr; k++) {
	    resist[i][j][k] = 0.;
	    csT[i][j][k] = Conductor::sp->GetCurrentSharingTemp(fieldB[i][j][k],I);
        if (k%2==0) {
		  if ( T[i][j][k]>=CriticalTemperature(fieldB[i][j][k])) {
            n++;
		    res += AlResistivity(preT[i][j][k], RRR_Al[i][j][k], fieldB[i][j][k]) * Dimension::dp / Conductor::areaAl;
			resist[i][j][k] += AlResistivity(preT[i][j][k], RRR_Al[i][j][k], fieldB[i][j][k]) * Dimension::dp / Conductor::areaAl;
		  }
		  else if ( T[i][j][k]<CriticalTemperature(fieldB[i][j][k]) && T[i][j][k]>=csT[i][j][k] ) {
            cs = (T[i][j][k] - csT[i][j][k]) / (CriticalTemperature(fieldB[i][j][k]) - T[i][j][k]);
			res += cs * AlResistivity(preT[i][j][k], RRR_Al[i][j][k], fieldB[i][j][k]) * Dimension::dp / Conductor::areaAl;
		    resist[i][j][k] += cs * AlResistivity(preT[i][j][k], RRR_Al[i][j][k], fieldB[i][j][k]) * Dimension::dp / Conductor::areaAl;
		  }
		}
	  }
	}
  }
  
  qchCell = n;
  totRes  = res;
}


double GetRRRFromData(std::string mat, int ix, int iy, int iz) {
  QSuperconduct* rad = new QSuperconduct();
  int zfactor = (Mesh::Mz  -1) / Radiation::nz;
  int pfactor = (Mesh::Mphi-1) / Radiation::np;
  int rfactor = (Mesh::Mr  -1) / Radiation::nr;
  
  int cz = (ix-1) / zfactor;
  int cp = (iy-1) / pfactor;
  int cr = (iz-1) / rfactor;

  if (cr>=Radiation::nr)
    cr = cr - 1;
  
  // flip the z and r direction data
  //cz = Radiation::nz - 1 - cz;
  //cr = Radiation::nr - 1 - cr;
  double RRR = rad->GetRadiationEffect(mat, Solenoid::operTime, Radiation::neutron[cz][cp][cr]);

  //std::cout << cz << " " << cp << " " << cr << " " << " " << ix << " " << iy << " " << iz << " " << RRR << " " << Radiation::neutron[cz][cp][cr] << std::endl;
  return RRR;
}

void TemperatureInitialization() {
  Conductor::sp->SetStrand();
  for (int i=0; i<=Mesh::Mz; i++) {
    for (int j=0; j<=Mesh::Mphi; j++) {
      for (int k=0; k<=Mesh::Mr; k++)
	    preT[i][j][k] = T[i][j][k];
	}
  }
}

bool SetNode() {
  bool done = false;

  for (int i=1; i<Mesh::Mz; i++) {
    for (int j=1; j<Mesh::Mphi; j++) {
      for (int k=1; k<Mesh::Mr; k++) {
        switch (k) {
		  // shell
          case 19:
		    r  [i][j][k] = Aluminium::inner + (k-2)*Aluminium::strip/2 + (k-1)*(Conductor::lr/2+Kapton::lz) + Shell::lr/2;
			phi[i][j][k] = Dimension::dp/2 + (j-1) * Dimension::dp;
			z  [i][j][k] = Dimension::dz/2 + (i-1) * Dimension::dz; 
			break;
		  // aluminium inner strip
          case 1:
		    r  [i][j][k] = Aluminium::inner/2;
			phi[i][j][k] = Dimension::dp/2 + (j-1) * Dimension::dp;
			z  [i][j][k] = Dimension::dz/2 + (i-1) * Dimension::dz;
			break;
		  case 2:
		  case 3:
		  case 4:
		  case 5:
		  case 6:
		  case 7:
		  case 8:
		  case 9:
		  case 10:
		  case 11:
          case 12:
		  case 13:
		  case 14:
		  case 15:
		  case 16:
		  case 17:
		  case 18:
		    r  [i][j][k] = Aluminium::inner + (k-2)*Aluminium::strip/2 + (k-1)*(Conductor::lr/2 + Kapton::lz);
			phi[i][j][k] = Dimension::dp/2 + (j-1) * Dimension::dp;
			z  [i][j][k] = Dimension::dz/2 + (i-1) * Dimension::dz;
			break;
		  default:
		    r  [i][j][k] = 0.;
			phi[i][j][k] = 0.;
			z  [i][j][k] = 0.;
			break;
		}
	  }
	}
  }
  
  done = true;

  return done;
}

bool PrintNode() {
  TGraph2D* gr2 = new TGraph2D();
  
  ofstream file;
  file.open("position.dat");
  file << "nz    nphi    nr    z [mm]    phi [mm]    r [mm] \n";

  int cc = 0;
  for (int i=1; i<Mesh::Mz; i++) {
    for (int j=1; j<Mesh::Mphi; j++) {
      for (int k=1; k<Mesh::Mr; k++) {
        gr2->SetPoint(cc, z[i][j][k], phi[i][j][k], r[i][j][k]);
		file << i << "  " << j << "  " << k << "  " << z[i][j][k]*1e+3 << "  " << phi[i][j][k]*1e+3 << "  " << r[i][j][k]*1e+3 << "\n";
		cc++;
	  }
	} 
  }
  
  TCanvas* c0 = new TCanvas("point", "point", 550, 500); 
  c0->SetTicks(1,1);
  
  gr2->SetMarkerSize(0.5);
  gr2->SetMarkerStyle(24);
  gr2->SetTitle("Node; Z [m]; #Phi [m]; R [m]");
  gr2->Draw("ap");
  c0->Print("point.eps"); 
  
  file.close();
  return true;
}

void SetGeometryParameter() {
  double k_Al, k_Cdt, k_tape;
  double term1, term2;
  //double C_Al;
  
  //std::cout << preT[1][1][8] << RRR_Al [1][1][8] << std::endl;
  for (int i=1; i<Mesh::Mz; i++) {
    for (int j=1; j<Mesh::Mphi; j++) {
      for (int k=1; k<Mesh::Mr; k++) {
		k_Al   = Aluminium::mat->ThermalConduct(preT[i][j][k], RRR_Al [i][j][k], field);
		k_Cdt  = Aluminium::mat->ThermalConduct(preT[i][j][k], RRR_Cdt[i][j][k], field);
		k_tape = Kapton::mat->ThermalConduct(preT[i][j][k], 0., 0.);
		
		C [i][j][k] = Aluminium::mat->SpecificHeat(preT[i][j][k]);
		
		// fill the thermal conductivity into array:
		switch (k) {
          // 1: layer with shell
          case  1:
			term1 = Kapton::lr + Conductor::lr/2 + Shell::lr/2;
			term2 = Kapton::lr/k_tape + Conductor::lr/2/k_Cdt + Shell::lr/2/k_Al;
			kr[i][j][k] = term1 / term2;
			kz[i][j][k] = k_Al;
			break;
		  // 19: inner aluminium strip:
		  case 19:
		    term1 = Aluminium::inner/2 + Kapton::lr + Conductor::lr/2;
			term2 = Aluminium::inner/2/k_Al + Kapton::lr/k_tape + Conductor::lr/2/k_Cdt;
		    kr[i][j][k] = term1 / term2;
			kz[i][j][k] = k_Al;
			break;
		  // conductor:
		  case  2:
		  case  4:
		  case  6:
		  case  8:
		  case 10:
		  case 12:
		  case 14:
		  case 16:
		  case 18:
	        term1 = Aluminium::strip/2 + Kapton::lr + Conductor::lr/2;
			term2 = Aluminium::strip/2/k_Al + Kapton::lr/k_tape + Conductor::lr/2/k_Cdt;
			kr[i][j][k] = term1 / term2;
			kz[i][j][k] = (Kapton::lz + Conductor::lz) / (Kapton::lz/k_tape + Conductor::lz/k_Cdt);
		    break;
		  // aluminium:
		  case  3:
		  case  5:
		  case  7:
		  case  9:
		  case 11:
		  case 13:
		  case 15:
		  case 17:
			term1 = Aluminium::strip/2 + Kapton::lr + Conductor::lr/2;
			term2 = Aluminium::strip/2/k_Al + Kapton::lr/k_tape + Conductor::lr/2/k_Cdt;
			kr[i][j][k] = term1 / term2;
            kz[i][j][k] = k_Al;
			break;
		  default: break;
		}
	  }
	}
  }
  
  //std::cout << k_Al << " " << k_Cdt << std::endl;
}

void SetQuenchSpot(int i, int j, int k, int qz, int qphi, int qr) {
  double heat;
  double csheat;
  
  if ( i==1 && j==1 && k==18 ) {
	heat = pow(current,2) * AlResistivity(preT[i][j][k], RRR_Al[i][j][k], field) * Dimension::dp / Conductor::aAl;
	Heat[i][j][k] = heat * Conductor::rArea / (Conductor::density * C[i][j][k] * Conductor::volume);
  }
  
  if ( preT[i][j][k]>=CriticalTemperature(fieldB[i][j][k]) ) {
    heat = pow(current,2) * AlResistivity(preT[i][j][k], RRR_Al[i][j][k], field) * Dimension::dp / Conductor::aAl;
	Heat[i][j][k] = heat * Conductor::rArea / (Conductor::density * C[i][j][k] * Conductor::volume);
  }
  else if ( preT[i][j][k]<CriticalTemperature(fieldB[i][j][k]) && preT[i][j][k]>=csT[i][j][k] ) {
    heat = pow(current,2) * AlResistivity(preT[i][j][k], RRR_Al[i][j][k], field) * Dimension::dp / Conductor::aAl;
    csheat = (preT[i][j][k] - csT[i][j][k]) / (CriticalTemperature(field) - csT[i][j][k]);
	Heat[i][j][k] = csheat * heat * Conductor::rArea / (Conductor::density * C[i][j][k] * Conductor::volume);
  }
  else if ( preT[i][j][k]<csT[i][j][k] )
    Heat[i][j][k] = 0.0;
  
}

void SetBoundary() {
  // r direction:
  for (int i=0;  i<=Mesh::Mz; i++){
    for (int j=0;  j<=Mesh::Mphi; j++){
	  T[i][j][0]        = T[i][j][1];
	  T[i][j][Mesh::Mr] = T[i][j][Mesh::Mr-1];
	}
  }
  // phi direction:
  for (int i=0;  i<=Mesh::Mz; i++){
	for (int k=0;  k<=Mesh::Mr; k++){
	  T[i][0][k]          = T[i][Mesh::Mphi-1][k];
	  T[i][Mesh::Mphi][k] = T[i][1][k];
	}
  }
  // z direction:	
  for (int j=0;  j<=Mesh::Mphi; j++){
	for (int k=0;  k<=Mesh::Mr; k++){
	  T[Mesh::Mz][j][k] = T[Mesh::Mz-1][j][k];
	  T[0][j][k]        = T[1][j][k];
	}
  }
}

void WriteQuenchVelocity() {
  // write time when temperature exceed to the critical temperature
  for (int i=1; i<Mesh::Mz; i++) {
    for (int j=1; j<Mesh::Mphi; j++) {
      for (int k=1; k<Mesh::Mr; k++) {
        if ( trigger[i][j][k]==false && T[i][j][k]>=CriticalTemperature(fieldB[i][j][k]) ) {
          qchTime[i][j][k] = Time;
		  trigger[i][j][k] = true;
		}
	  }
	}
  }
}


#endif
