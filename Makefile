TARGET   = SimQuench.exe
TARGET2  = RunAnalysis.exe
CXX      = g++
OBJECTS  = QuenchMain.o QSuperconduct.o QMaterial.o FillData.o
OBJECTS2 = runPlot.o QAnalysis.o QMakePlot.o
CXXLIBS  = 
CXXFLAGS = -Wall -O3

ROOTFLAGS = `root-config --cflags`
ROOTLIBS  = `root-config --evelibs`

CXXLIBS  += $(ROOTLIBS)
CXXFLAGS += $(ROOTFLAGS)

INSTALLDIR = ./bin

.PHONY: all
all: $(TARGET) $(TARGET2)

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXLIBS) $^ -o $@

$(TARGET2): $(OBJECTS2)
	$(CXX) $(CXXLIBS) $^ -o $@

.o:.cpp
	$(CXX) $(CXXFLAGS) -c $<

.PHONY: clean
clean:
	$(RM) $(TARGET) $(TARGET2) $(OBJECTS) $(OBJECTS2)

.PHONY: install
install:
	mkdir -p $(INSTALLDIR)
	cp -p $(TARGET) $(TARGET2) $(INSTALLDIR)
