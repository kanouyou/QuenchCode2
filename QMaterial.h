//#ifdef MATERIAL_H_INCLUDED
//#define MATERIAL_H_INCLUDED
#include <string>

class QMaterial {
	
	public:
		QMaterial(std::string name);
		~QMaterial() { }
		void SetMetal(double RRR, double magfld=0.0);
		double Resistivity(double temp, double materialRRR, double magfld=0.0);
		double ThermalConduct(double temp, double materialRRR, double magfld=0.0);
		double SpecificHeat(double temp);
		void Print();
		// Get Parameter from Material:
		double GetResistivity() const { return resistivity; }
		std::string GetMaterial() const { return name; }
		double GetRRR() const { return RRR; }
		double GetField() const { return field; }

	private:
		std::string name;
		double temperature;
		double RRR;
		double field;
		double thermal;
		double resistivity;
		double specificHeat;
		// Electric resistivity:
		double AlResistivity(double temp, double rrr, double field);
		double CuResistivity(double temp, double rrr, double field);
		// Thermal conductivity:
		double AlThermalConduct(double T, double rrr, double B);
		double CuThermalConduct(double T, double rrr, double B);
		double KaptonThermalConduct(double T);
		// Specific heat:
		double AlSpecificHeat(double T);
		double CuSpecificHeat(double T);
		double NbTiSpecificHeat(double T, double B=0.0);
		double KaptonSpecificHeat(double T);
};

//#endif
