#ifndef QMakePlot_HH
#define QMakePlot_HH

#include "QAnalysis.h"

class QMakePlot : public QAnalysis {
  
  public:
    QMakePlot(std::string, int, int, int);
	~QMakePlot();
    bool ShowTemperatureDistribution();
    bool ShowQuenchVelocity();
    bool ShowTemperature();
    bool ShowElectric();
    bool ShowDeltaTemp();

  private:
    int zbin;
	int pbin;
	int rbin;
    int* spot;
};

#endif
