#include <iostream>
#include <string>
#include <cmath>
#include "QMaterial.h"

using namespace std;

// constructor
QMaterial::QMaterial(string materialName) {
	name = materialName;
}


void QMaterial::SetMetal(double RRR, double magfld) {
	this->RRR = RRR;
	field = magfld;
}


double QMaterial::Resistivity(double temp, double materialRRR, double magfld) {
	// give the input parameter to private
	temperature = temp;
	RRR = materialRRR;
	field = magfld;

	// separate the material
	static double res;
	if (name=="Al" || name=="Aluminium") res = AlResistivity(temp, RRR, field);
	else if (name=="Cu" || name=="Copper") res = AlResistivity(temp, RRR, field);
	else cout << "Error: there is no this kind of material!" << endl;
	
	resistivity = res;
	return res;
}


double QMaterial::AlResistivity(double T, double rrr, double B) {
	double rho0, rhoi, rhoi0, rho, rhoB;
	// fitting parameter for Al resistivity
	const int n = 7;
	double p[n] = {1.671e-17, 4.36, 2.841e+10, 1.18, 64.0, 4.428, 1.2031};
	double rhoRT = 2.75e-8;			// Aluminium resistivity at room temperature [Ohm*m]
	// fitting parameter for Al magnetoresistance
	const int nb = 5;
	double pb[nb] = {3.62857, 2.90419e-5, 3.79649e+6, 10975.9, 0.761609};

	rho0 = rhoRT / rrr;
	rhoi = p[0] * pow(T, p[1]) / (1 + p[0] * p [2] * pow(T, p[1]-p[3]) * exp(-pow(p[4]/T, p[5])));
	rhoi0 = p[6] * rhoi * rho0 / (rho0 + rhoi);
	rho = rho0 + rhoi + rhoi0;
	if (B>0.0) {
		double h = B * 10 * rhoRT / rho;
		rhoB = h*h*(pb[0] - pb[1]*h) * rho / (pb[2] + pb[3]*h + pb[4]*h*h) + rho;
		return rhoB;
	}
	else{
		return rho;
	}
}


double QMaterial::CuResistivity(double T, double rrr, double B) {
	double rhoRT, rho0, rhoi, rhoi0, rho;
// fitting parameter for copper:
	const int nf = 7;
	double p[nf] = {1.171e-17, 4.49, 3.841e+10, 1.14, 50.0, 6.428, 0.4531};
// Copper resistivity at room temperature [Ohm*m]:
	rhoRT = 1.553e-8;
	rho0 = rhoRT/rrr;		// resistivity at cryogenic temperature
	rhoi = p[0] * pow(T, p[1]) / (1 + p[0] * p [2] * pow(T, p[1]-p[3]) * exp(-pow(p[4]/T, p[5])));
	rhoi0 = p[6] * rhoi * rho0 / (rho0 + rhoi);
	rho = rho0 + rhoi + rhoi0;
// Magnetoresistivity:
	if (B>0.0) {
		const int nb = 5;
		double pb[nb] = {-2.662, 0.3168, 0.6229, -0.1839, 0.01827};
		double ax = pb[0]*pow(log10(RRR*B),0) + pb[1]*pow(log10(RRR*B),1) + pb[2]*pow(log10(RRR*B),2) + pb[3]*pow(log10(RRR*B),3)+pb[4]*pow(log10(RRR*B),4);
		double rhoB = rho * (1 + pow(10,ax));
		return rhoB;
	}
	else {
		return rho;
	}
}


double QMaterial::ThermalConduct(double temp, double materialRRR, double magfld) {
	double k = 0.;;
	RRR = materialRRR;
	temperature = temp;
	field = magfld;
	if (name=="Kapton") k = KaptonThermalConduct(temp);
	else if (name=="Al" || name=="Aluminium") k = AlThermalConduct(temp, RRR, field);
	else if (name=="Cu" || name=="Copper") k = CuThermalConduct(temp, RRR, field);
	else cout << "Error!" << endl;
	thermal = k;
	return k;
}


double QMaterial::AlThermalConduct(double T, double rrr, double B) {
	double k;
	double Lwf = 2.44e-8;
	double rho = AlResistivity(T, rrr, B);
	k = Lwf * T / rho;
	return k;
}


double QMaterial::CuThermalConduct(double T, double rrr, double B) {
// from NIST fitting equation:
	double k, kB, W0, Wi, Wi0, rho;
	double beta = 0.634/rrr;
	double betar = beta/0.0003;
	const int n = 7;
	double p[n] = {1.754e-8, 2.763, 1102.0, -0.165, 70.0, 1.756, 0.838/pow(betar, 0.1661)};
	W0 = beta / T;
	Wi = p[0]*pow(T, p[1]) / (1 + p[0]*p[2]*pow(T, (p[1]+p[3]))*exp(-pow(p[4]/T, p[5])));
	Wi0 = p[6] * Wi * W0 / (Wi + W0);
	k = 1 / (W0 + Wi + Wi0);
	// magnetic effect:
	if (B>0.0) {
	// k(B=0), rho(B)
		rho = CuResistivity(T, rrr, B);
		kB = CuResistivity(T, rrr, 0.0) * k / rho;
		return kB;
	}
	else {
		return k;
	}
}


double QMaterial::KaptonThermalConduct(double T) {
	double k, ax;
	ax = 0.0;
	const int n = 8;
	double a[n] = {5.73101, -39.5199, 79.9313, -83.8572, 50.9157, -17.9835, 3.42413, -0.27133};
// NIST equation:
	if (T<4.3) k = 0.00378 + 0.00161*T;
	else {
		for (int i=0; i<n; i++) {
			ax = ax + a[i]*pow(log10(T), i);
		}
		k = pow(10, ax);
	}
	/*const int n = 4;
	double a[n] = {15.549709e-3, -17.417525e-3, 82.690158e-4, -81.952869e-5};
	if (T<4.3) k = 0.00378 + 0.00161*T;
	else {
		k = 0.0;
		for (int i=0; i<n; i++) {
			k = k + a[i]*pow(log(T+1), i);
		}
	}*/
	return k;
}


double QMaterial::SpecificHeat(double temp) {
	double SH = 0.;
	if (name=="Al" || name=="Aluminium") SH = AlSpecificHeat(temp);
	else if (name=="Cu" || name=="Copper") SH = CuSpecificHeat(temp);
	else if (name=="NbTi") SH = NbTiSpecificHeat(temp, field);
	else if (name=="Kapton") SH = KaptonSpecificHeat(temp);
	else cout << "Error!: no such kind of material..." << endl;
	specificHeat = SH;
	return specificHeat;
}


double QMaterial::AlSpecificHeat(double T) {
	double C;
	double point1 = 22.67; double point2 = 46.0;
	const int n = 4;
	double a1[n] = {-0.207489, 0.165759, -0.0142572, 0.00146459};
	double a2[n] = {7.88e-13, 6.93201, -0.07139, 46.4363};
	double a3[n] = {6.273517, -0.5469, 0.000925, -156.932};
	if (T<point1) C = a1[0] + a1[1]*T + a1[2]*pow(T, 2) + a1[3]*pow(T, 3);
	else if (T>=point1 && T<point2) C = a2[0] * pow(T, a2[1]) * exp(a2[2]*T) * exp(a2[3]/T) * 4.186e+3;
	else C = a3[0] * pow(T, a3[1]) * exp(a3[2]*T) * exp(a3[3]/T) * 4.186e+3;
	return C;
}


double QMaterial::CuSpecificHeat(double T) {
	double C;
// Fitted by YeYang:
	/*double point1 = 23.; double point2 = 55.; double point3 = 250.;
	const int n = 4;
	double p[n];
	if (T<point1) {
		p[0] = -0.0104251;
		p[1] = 0.0832507;
		p[2] = -0.0118811;
		p[3] = 0.00131157;
	}
	else if (T>=point1 && T<point2) {
		p[0] = 25.2218;
		p[1] = -3.41019;
		p[2] = 0.144053;
		p[3] = -0.000939968;
	}
	else if (T>=point2 && T<point3) {
		p[0] = -158.094;
		p[1] = 6.40481;
		p[2] = -0.0273429;
		p[3] = 4.08432e-5;
	}
	else if (T>=point3) {
		p[0] = 147.891;
		p[1] = 1.78704;
		p[2] = -0.00452361;
		p[3] = 4.21855e-6;
	}
	else cout << "Error: out of range" << endl;
	C = p[0] + p[1]*T + p[2]*pow(T, 2) + p[3]*pow(T, 3);*/
// Fitted by NIST (unit: J/K/kg)
	const int n = 8;
	double a[n] = {-1.91844, -0.15973, 8.61013, -18.996, 21.9661, -12.7328, 3.54322, -0.3797};
	double ax = 0.0;
	for (int i=0; i<n; i++) ax = ax + a[i] * pow(log10(T), i);
	C = pow(10, ax);
	return C;
}


double QMaterial::NbTiSpecificHeat(double T, double B) {
	const int n = 5;
	double C = 0.0; double a[n];
	double Tc = 9.4;		// the Tc should change with magnetic field
	double rho = 6538.;
	if (T>0 && T<Tc) {
		a[0] = 0.0;
		a[1] = 64. * B;
		a[2] = 0.0;
		a[3] = 49.1;
		a[4] = 0.0;
	}
	else if (T>=Tc && T<28.358) {
		a[0] = 0.0;
		a[1] = 928.;
		a[2] = 0.0;
		a[3] = 16.24;
		a[4] = 0.0;
	}
	else if (T>=28.358 && T<50.99) {
		a[0] = 41383;
		a[1] = -7846.1;
		a[2] = 553.71;
		a[3] = 11.9838;
		a[4] = -0.2177;
	}
	else if (T>=50.99 && T<165.8) {
		a[0] = -1.53e+6;
		a[1] = 83022.;
		a[2] = -716.3;
		a[3] = 2.976;
		a[4] = -0.00482;
	}
	else if (T>=165.8 && T<496.54) {
		a[0] = 1.24e+6;
		a[1] = 13706;
		a[2] = -51.66;
		a[3] = 0.09296;
		a[4] = -6.29e-5;
	}
	else if (T>=496.54) {
		a[0] = 2.45e+6;
		a[1] = 955.5;
		a[2] = -0.257;
		a[3] = 0.;
		a[4] = 0.;
	}
	else cout << "Error: out of the range" << endl;
	for (int i=0; i<n; i++) C = C + a[i]*pow(T,i);
	return C/rho;
}


double QMaterial::KaptonSpecificHeat(double T) {
	const int n = 8;
	double C; double ax = 0.0;
	double a[n] = {-1.3684, 0.65892, 2.8719, 0.42651, -3.0088, 1.9558, -0.51998, 0.051574};
	for (int i=0; i<n; i++) ax = ax + a[i]*pow(log10(T), i);
	C = pow(10, ax);
	return C;
}


void QMaterial::Print() {
	cout << " * Material: " << name << endl;
	cout << "   -> temperature: " << temperature << "[K]; " << "magnetic field: " << field << "[T]" << endl;
	cout << "   -> RRR: " << RRR << "; " << "resistivity: " << resistivity << "[Ohm*m]" << endl;
}


