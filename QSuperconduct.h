#ifndef ___QSuperconduct
#define ___QSuperconduct

#include <string>

class QSuperconduct {
// this class is to give the critical temperature and current for NbTi
  public:
    QSuperconduct() {}
	~QSuperconduct() {}
	void SetStrand(double radius=1.15/2, int no=14);
	void SetConductor(double width=15.00, double thickness=4.73);
	void SetRatio(double R_al=7.3, double R_cu=1.0, double R_nbti=0.9);
	void SetParameter(double B, double T);
	void SetOperatingCurrent(double Iop);
	void SetTemperature(double temp);
	void SetField(double field);
	double GetCriticalCurrent();
	double GetCriticalTemperature();
	double GetCurrentSharingTemp(double field, double iop);
	double GetRadiationEffect(std::string mat, int time, double flux);
	//double GetOperatingCurrent(double tcs);

  protected:
    double Iop;
    double B, T;
	double Tcs;
	double Asc, R_al, R_cu, R_nbti;
	double w_cdt, t_cdt;

  private:
    double CriticalCurrent(double field, double temp);
	double CriticalTemperature(double field);
};

#endif
