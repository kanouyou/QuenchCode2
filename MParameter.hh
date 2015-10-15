double AlResistivity(double& T, double& RRR, double& B){
	double P1 = 1.671e-17;
	double P2 = 4.36;
	double P3 = 2.841e+10;
	double P4 = 1.18;
	double P5 = 64.;
	double P6 = 4.428;
	double P7 = 1.2031;
	double rho0 = 2.45e-8 / RRR;
	double rhoi = P1 * pow(T,P2) / (1 + P1 * P3 * pow(T, P2-P4) * exp(-pow(P5/T, P6)));
	double rhoi0 = P7 * rhoi * rho0 / (rho0 + rhoi);
	double rho = rho0 + rhoi + rhoi0;
	double rhoB;
	if (B==0.0){
		rhoB = rho;
	}
	else{
		double h = B*10*2.75e-8/rho;
		rhoB = h*h*(3.62857-2.90419e-5*h)*rho/(3.79649e+6+10975.9*h+0.761609*h*h)+rho;
	}
	return rhoB;
}

double AlSpecificHeat(double T){
	float point1 = 22.67;
	float point2 = 46.;
	float c0 = -0.207489;
	float c1 = 0.165759;
	float c2 = -0.0142572;
	float c3 = 0.00146459;
	if (T < point1){
		double C = c0 + c1*T + c2*pow(T,2) + c3*pow(T,3);
		return C;
	}
	else if (T >= point1 && T < point2){
		double a = 7.88e-13;
		double b = 6.93201;
		double c = -0.07139;
		double d = 46.4363;
		double C = a * pow(T,b) * exp(c*T) * exp(d/T) * 4.186 * 1e+3;
		return C;
	}
	else{
		double a = 6.273517;
		double b = -0.5469;
		double c = 0.000925;
		double d = -156.932;
		double C = a * pow(T,b) * exp(c*T) * exp(d/T) * 4.186 * 1e+3;
		return C;
	}
}

double CuResistivity(double T, double RRR, double B){
	double P1 = 1.171e-17;
	double P2 = 4.49;
	double P3 = 3.841e+10;
	double P4 = 1.14;
	double P5 = 50.0;
	double P6 = 6.428;
	double P7 = 0.4531;
	double a0 = -2.662;
	double a1 = 0.3168;
	double a2 = -0.1839;
	double a3 = 0.01827;
	double rho0 = 1.553e-8 / RRR;
	double rhoi = P1 * pow(T,P2) / (1 + P1 * P3 * pow(T,P2-P4) * exp(-pow(P5/T,P6)));
	double rhoi0 = P7 * rhoi * rho0 / (rho0 + rhoi);
	double rho = rho0 + rhoi + rhoi0;
	double rhoB;
	if (B==0.){
		rhoB = rho;
	}
	else{
		double x = 1.553e-8 * B / rho;
		double ax = a0 + a1*log10(x) + a2*pow(log10(x),2) + a3*pow(log10(x),3);
		rhoB = rho*(1+pow(10,ax));
	}
	return rhoB;
}

double CuSpecificHeat(double T){
	double C;
	/*float point1 = 23.;
	float point2 = 55.;
	float point3 = 250.;
	if (T < point1){
		double p0 = -0.104251;
		double p1 = 0.0832507;
		double p2 = -0.0118811;
		double p3 = 0.00131157;
		C = p0 + p1*T + p2*pow(T,2) + p3*pow(T,3);
	}
	else if (T >= point1 && T < point2){
		double p0 = 25.2218;
		double p1 = -3.41019;
		double p2 = 0.144053;
		double p3 = -0.000939968;
		C = p0 + p1*T + p2*pow(T,2) + p3*pow(T,3);
	}
	else if (T >= point2 && T < point3){
		double p0 = -158.094;
		double p1 = 6.40481;
		double p2 = -0.0273429;
		double p3 = 4.08432e-5;
		C = p0 + p1*T + p2*pow(T,2) + p3*pow(T,3);
	}
	else if (T >= point3){
		double p0 = 147.891;
		double p1 = 1.76704;
		double p2 = -0.00452361;
		double p3 = 4.21855e-6;
		C = p0 + p1*T + p2*pow(T,2) + p3*pow(T,3);
	}*/
	const int n = 8;
    double a[n] = {-1.91844, -0.15973, 8.61013, -18.996, 21.9661, -12.7328, 3.54322, -0.3797};
	double ax = 0.;
	for (int i=0; i<n; i++) ax += a[i] * pow(log10(T), i);
	C = pow(10, ax);
	return C;
}

double NbTiSpecificHeat(double T, double B){
	double C  = 0.;
	double Tc = 6.5;
	double rho = 6538;
	if (T>0 && T<Tc){
		double a0 = 0;
		double a1 = 64 * B;
		double a2 = 0;
		double a3 = 49.1;
		double a4 = 0;
		C = (a0*pow(T,0) + a1*pow(T,1) + a2*pow(T,2) + a3*pow(T,3) + a4*pow(T,4))/rho;
	}
	else if (T>=Tc && T<28.358){
		double a0 = 0;
		double a1 = 928;
		double a2 = 0;
		double a3 = 16.24;
		double a4 = 0;
		C = (a0*pow(T,0) + a1*pow(T,1) + a2*pow(T,2) + a3*pow(T,3) + a4*pow(T,4))/rho;
	}
	else if (T>=28.358 && T<50.99){
		double a0 = 41383;
		double a1 = -7846.1;
		double a2 = 553.71;
		double a3 = 11.9838;
		double a4 = -0.00482;
		C = (a0*pow(T,0) + a1*pow(T,1) + a2*pow(T,2) + a3*pow(T,3) + a4*pow(T,4))/rho;
	}
	else if (T>=50.99 && T<165.8){
		double a0 = -1.53e+6;
		double a1 = 83022;
		double a2 = -716.3;
		double a3 = 2.976;
		double a4 = -0.00482;
		C = (a0*pow(T,0) + a1*pow(T,1) + a2*pow(T,2) + a3*pow(T,3) + a4*pow(T,4))/rho;
	}
	else if (T>=165.8 && T<496.54){
		double a0 = 1.24e+6;
		double a1 = 13706;
		double a2 = -51.66;
		double a3 = 0.09296;
		double a4 = -6.29e-5;
		C = (a0*pow(T,0) + a1*pow(T,1) + a2*pow(T,2) + a3*pow(T,3) + a4*pow(T,4))/rho;
	}
	else if (T>=496.54){
		double a0 = 2.45e+6;
		double a1 = 955.5;
		double a2 = -0.257;
		double a3 = 0;
		double a4 = 0;
		C = (a0*pow(T,0) + a1*pow(T,1) + a2*pow(T,2) + a3*pow(T,3) + a4*pow(T,4))/rho;
	}
	return C;
}

double ConductorSpecificHeat(double C_Al, double C_Cu, double C_NbTi){
	double R_Al = 7.3/9.3;
	double R_Cu = 1.0/9.3;
	double R_NbTi = 1.0/9.3;
	double C = C_Al * R_Al + C_Cu * R_Cu + C_NbTi * R_NbTi;
	return C;
}

double ConductorResistivity(double rho_Al, double rho_SC){
	double R_Al = 7.3/9.3;
	double R_SC = 1.0/9.3;
	double reverse_rho = R_Al/rho_Al + R_SC/rho_SC;
	return 1/reverse_rho;
}

