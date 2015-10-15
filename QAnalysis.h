#ifndef QAnalysis_HH
#define QAnalysis_HH

#include <string>
#include <TGraph.h>
#include <TFile.h>
#include <TTree.h>

class QAnalysis {
  public:
    QAnalysis(std::string fileName);
	~QAnalysis() {}
    void PlotTemperatureDistribution(double qtime, int numPhi, bool scale=false);
    void PlotSpotTemperature(int z, int phi, int r);
    void PlotCurrent();
	void PlotQuenchedNumber();
	void PlotResistance();
	void PlotVoltage();
	void PlotPower();
	void PlotVelocityDistribution(int phi);
	TGraph* GetTemperatureGraph(int z, int phi, int r);
    double  GetMinimum(std::string name);
	double  GetMaximum(std::string name);
    double  GetMinTemperature(double fTime);
	double  GetMaxTemperature(double fTime);
	double* GetQuenchPoint(int z, int phi, int r);
    void    PlotDeltaTemperature();
	void    PlotQuenchVelocity(int phi);

  private:
    TTree* ReadTree(TFile* file);

  protected:
    TFile* file;
    int*   posID;
    double* fpos;
	double time;
	double current;
	double field;
	double R;
	int    qchCell;
	double temp;
	double qchTime;
    double zbin;
	double zmax;
};

#endif
