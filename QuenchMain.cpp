#include <iostream>
#include <cmath>
#include <TApplication.h>
#include <TSystem.h>
#include <TROOT.h>
#include "QuenchMain.hh"
#include "FillData.hh"

using namespace std;

int main(int argc, char** argv){
  
  bool setnode = SetNode();
  
  if ( setnode )
    PrintNode();
  
  for (int i=1; i<Mesh::Mz; i++) {
	for (int j=1; j<Mesh::Mphi; j++) {
	  for (int k=1; k<Mesh::Mr; k++) {
	    //RRR_Cdt[i][j][k] = 50;
		//RRR_Al [i][j][k] = 50;
		// RRR data from PHITS approach
		RRR_Cdt[i][j][k] = GetRRRFromData("Conductor", i, j, k);
	    RRR_Al [i][j][k] = GetRRRFromData("Aluminium", i, j, k);
	  }
	}
  }
  
  //ofstream vfile;
  //vfile.open("velocity.dat");

  for (int i=0;  i<=Mesh::Mz; i++){
	for (int j=0; j<=Mesh::Mphi; j++){
	  for (int k=0; k<=Mesh::Mr; k++){
		T      [i][j][k] = 4.5;
		preT   [i][j][k] = 4.5;
	    Heat   [i][j][k] = 0.;
	    trigger[i][j][k] = false;
	    // field
		//fieldB [i][j][k] = Field::factor * Field::B[k/(Field::n-1)];
	  }
	}
  }
	
  double Qr, Qrr, Qp, Qpp, Qz, Qzz, Q;
  double kAl, kCdt, ktape, kShell, term1, term2, term3;
  
  FillData* dataflow = new FillData();

  for (int t=1; t<=Mesh::Mt; t++){
    // current time:
    Time = Mesh::t0 + t * Mesh::dt;
    // initialization:
    TemperatureInitialization();
    // calculate the resistance of quenched cell:
	GetConductorResistance(current);
    // current and field decay:
    current = Solenoid::current * exp(-(Solenoid::resistor + totRes) * Time / Solenoid::induct);
	//field   = Solenoid::field   * exp(-(Solenoid::resistor + totRes) * Time / Solenoid::induct);
	// magnetic field
    FillMagneticField();

	if (t%Mesh::pt==0) printf("numOfQuench= %i ", qchCell);
    
	// Geometry and Parameter Initialization
	for (int i=1; i<Mesh::Mz; i++) {
      for (int j=1; j<Mesh::Mphi; j++) {
        for (int k=1; k<Mesh::Mr; k++) {
	      kAl   = Aluminium::mat->ThermalConduct(preT[i][j][k], RRR_Al [i][j][k], fieldB[i][j][k]);
		  kShell= Aluminium::mat->ThermalConduct(preT[i][j][k], Shell::RRR, fieldB[i][j][k]);
		  kCdt  = Aluminium::mat->ThermalConduct(preT[i][j][k], RRR_Cdt[i][j][k], fieldB[i][j][k]);
		  ktape = Kapton::mat->ThermalConduct(preT[i][j][k], 0., 0.);
		  
		  C [i][j][k] = Aluminium::mat->SpecificHeat(preT[i][j][k]);
		  dz [i][j][k] = pow(Dimension::dz, 2);
		  dzz[i][j][k] = pow(Dimension::dz, 2);
          dp [i][j][k] = pow(Dimension::dp, 2);
		  dpp[i][j][k] = pow(Dimension::dp, 2);
          rho[i][j][k] = Aluminium::density;

		  // fill the thermal conductivity into array:
		  switch (k) {
            // 1: layer with shell
            case  1:
			  term1 = Kapton::lr + Conductor::lr/2 + Shell::lr/2;    // = dl from Shell cell to conductor cell
			  //term2 = Kapton::lr/ktape + Conductor::lr/2/kCdt + Shell::lr/2/kAl;
			  term2 = Kapton::lr/ktape + Shell::lr/2/kShell;
			  term3 = Kapton::lr + Shell::lr/2;
			  kr [i][j][k] = term3 / term2;
		   	  kz [i][j][k] = kShell;
			  kp [i][j][k] = kShell;
			  dr [i][j][k] = term1 * Shell::lr;
			  drr[i][j][k] = Shell::lr * Shell::lr;
			  break;
		    // 19: inner aluminium strip:
		    case 19:
		      term1 = Aluminium::inner/2 + Kapton::lr + Conductor::lr/2;
			  //term2 = Aluminium::inner/2/kAl + Kapton::lr/ktape + Conductor::lr/2/kCdt;
			  term2 = Aluminium::inner/2/kAl + Kapton::lr/ktape;
			  term3 = Aluminium::inner/2 + Kapton::lr;
		      kr [i][j][k] = term3 / term2;
			  kz [i][j][k] = kAl;
			  kp [i][j][k] = kAl;
			  dr [i][j][k] = term1 * Aluminium::inner;
			  drr[i][j][k] = Aluminium::inner * Aluminium::inner;
			  if (i==Mesh::Mz-1) {
                dz [i][j][k] = (Shell::lHe + Dimension::dz/2) * Dimension::dz;
				dzz[i][j][k] = Dimension::dz * Dimension::dz;
			  }
			  else if (i==1) {
                dz [i][j][k] = Dimension::dz * Dimension::dz;
				dzz[i][j][k] = (Shell::lHe + Dimension::dz/2) * Dimension::dz;
			  }
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
			  //term2 = Aluminium::strip/2/kAl + Kapton::lr/ktape + Conductor::lr/2/kCdt;
			  term2 = Kapton::lr/ktape + Conductor::lr/2/kCdt;
			  term3 = Kapton::lr + Conductor::lr/2;
			  kr [i][j][k] = term3 / term2;
			  kz [i][j][k] = (Kapton::lz + Conductor::lz) / (Kapton::lz/ktape + Conductor::lz/kCdt);
		      kp [i][j][k] = kCdt;
			  dr [i][j][k] = term1 * (Conductor::lr + Kapton::lr);
			  drr[i][j][k] = term1 * (Conductor::lr + Kapton::lr);
			  rho[i][j][k] = Conductor::density;
			  SetQuenchSpot(i, j, k, 1, 1, 10);
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
			  //term2 = Aluminium::strip/2/kAl + Kapton::lr/ktape + Conductor::lr/2/kCdt;
			  term2 = Aluminium::strip/2/kAl + Kapton::lr/ktape;
			  term3 = Aluminium::strip/2 + Kapton::lr;
			  kr [i][j][k] = term3 / term2;
              kz [i][j][k] = kAl;
			  kp [i][j][k] = ktape;
			  dr [i][j][k] = term1 * Aluminium::strip;
			  drr[i][j][k] = term1 * Aluminium::strip;
			  if (i==Mesh::Mz-1) {
                dz [i][j][k] = (Shell::lHe + Dimension::dz/2) * Dimension::dz;
				dzz[i][j][k] = Dimension::dz * Dimension::dz;
			  }
			  else if (i==1) {
                dz [i][j][k] = Dimension::dz * Dimension::dz;
				dzz[i][j][k] = (Shell::lHe + Dimension::dz/2) * Dimension::dz;
			  }
			  break;
		    default: break;
		  }
	    }
	  }
    }

	//SetGeometryParameter();

	for (int i=1;  i<Mesh::Mz; i++){
	  for (int j=1; j<Mesh::Mphi; j++){
		for (int k=1; k<Mesh::Mr; k++){	
		  
          Qz  = kz[i][j][k] * (preT  [i][j][k] - preT[i+1][j][k]) / (dz [i][j][k] * rho[i][j][k] * C[i][j][k]);
	      Qzz = kz[i][j][k] * (preT[i-1][j][k] - preT  [i][j][k]) / (dzz[i][j][k] * rho[i][j][k] * C[i][j][k]);
		  Qp  = kp[i][j][k] * (preT  [i][j][k] - preT[i][j+1][k]) / (dp [i][j][k] * rho[i][j][k] * C[i][j][k]);
		  Qpp = kp[i][j][k] * (preT[i][j-1][k] - preT  [i][j][k]) / (dpp[i][j][k] * rho[i][j][k] * C[i][j][k]);
		  Qr  = kr[i][j][k] * (preT  [i][j][k] - preT[i][j][k+1]) / (dr [i][j][k] * rho[i][j][k] * C[i][j][k]);
		  Qrr = kr[i][j][k] * (preT[i][j][k-1] - preT  [i][j][k]) / (drr[i][j][k] * rho[i][j][k] * C[i][j][k]);

		  Q = Qzz + Qpp + Qrr - Qz - Qp - Qr + Heat[i][j][k];
		  T[i][j][k] = preT[i][j][k] + Mesh::dt * Q;
        }
	  }
	}
	
	// setup boundary:
	SetBoundary();
	
	// check quench velocity:
    WriteQuenchVelocity();

	double T1 = T[1][1][18];
	double T2 = T[Mesh::Mz-1][1][18];
	double T3 = T[1][(Mesh::Mphi-1)/2][18];
	double T4 = T[1][1][2];
		
	if (t%Mesh::pt==0) {
	  printf("time= %.2f [sec]  T1= %.2f  T2= %.2f  T3= %.2f  T4= %.2f [K]  I= %.1f R= %.2e P= %.2e Pcell=%.2e B= %.2f\n",
	         Time, T1, T2, T3, T4, current, totRes, current*totRes, Heat[1][1][18], fieldB[1][1][18]);
	  for (int i=1; i<Mesh::Mz; i++) {
        for (int j=1; j<Mesh::Mphi; j++) {
          for (int k=1; k<Mesh::Mr; k++)
            dataflow->FillTree(Time, i, j, k, r[i][j][Mesh::Mr-k], phi[i][j][Mesh::Mr-k], z[i][j][Mesh::Mr-k],
			                   current, fieldB[i][j][k], resist[i][j][k], totRes, qchCell, T[i][j][k], qchTime[i][j][k]);
		}
	  }
	}
   
    /*if ( t%Mesh::pt==0 ) {
     //gSystem->ProcessEvents();
	 int ccc = 0;
     for (int i=1; i<Mesh::Mz; i++) {
       for (int j=1; j<Mesh::Mr; j++) {
         gr->SetPoint(ccc, z[i][1][j], r[i][1][j], T[i][1][j]);
         ccc++;
	   }
	 }
	 gr->Draw("colz");
	 c0->Modified();
	 c0->Update();
	}*/
  } 
    
  dataflow->Close();
  
  //app->Run();

  return 0;
}




