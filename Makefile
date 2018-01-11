CPP=g++ 
CPPFLAGS=-Wall -O2 -fopenmp

LDFLAGS=-lnetcdf_c++ -fopenmp
RM=rm -rf

DESTDIR=$(HOME)/bin

common_src= Constants.cpp  EqDate.cpp  GridConstruction.cpp LagrangianEngine.cpp  ParseParameters.cpp  VectorXYZ.cpp  VFlow.cpp  VTKdump.cpp mt19937ar.cpp mt19937arParallel.cpp
common_obj=$(common_src:.cpp=.o) 
common_dep=$(common_obj:.o=.d)  # one dependency file for each source

#InFlux
InFlux_src=InFlux.cpp
InFlux_obj=$(InFlux_src:.cpp=.o) 
InFlux_dep=$(InFlux_obj:.o=.d)  # one dependency file for each source

#InFluxCG
InFluxCG_src=InFluxCG.cpp
InFluxCG_obj=$(InFluxCG_src:.cpp=.o) 
InFluxCG_dep=$(InFluxCG_obj:.o=.d)  # one dependency file for each source

#InFluxNum
InFluxNum_src=InFluxNum.cpp
InFluxNum_obj=$(InFluxNum_src:.cpp=.o) 
InFluxNum_dep=$(InFluxNum_obj:.o=.d)  # one dependency file for each source

#Dproj
Dproj_src=Dproj.cpp
Dproj_obj=$(Dproj_src:.cpp=.o) 
Dproj_dep=$(Dproj_obj:.o=.d)  # one dependency file for each source

#Dproj2
Dproj2_src=Dproj2.cpp
Dproj2_obj=$(Dproj2_src:.cpp=.o) 
Dproj2_dep=$(Dproj2_obj:.o=.d)  # one dependency file for each source

#DprojCG
DprojCG_src=DprojCG.cpp
DprojCG_obj=$(DprojCG_src:.cpp=.o) 
DprojCG_dep=$(DprojCG_obj:.o=.d)  # one dependency file for each source

#Dproj
chist_src=chist.cpp
chist_obj=$(chist_src:.cpp=.o) 
chist_dep=$(chist_obj:.o=.d)  # one dependency file for each source

.PHONY: all InFlux InFluxCG InFluxNum Dproj DprojCG Dproj2 chist

all: all InFlux InFluxNum InFluxCG Dproj Dproj2 DprojCG chist


InFlux: $(common_obj) $(InFlux_obj)
	$(CPP) -o $@ $^ $(LDFLAGS)

InFluxCG: $(common_obj) $(InFluxCG_obj)
	$(CPP) -o $@ $^ $(LDFLAGS)

InFluxNum: $(common_obj) $(InFluxNum_obj)
	$(CPP) -o $@ $^ $(LDFLAGS)

Dproj: $(common_obj) $(Dproj_obj)
	$(CPP) -o $@ $^ $(LDFLAGS)

Dproj2: $(common_obj) $(Dproj2_obj)
	$(CPP) -o $@ $^ $(LDFLAGS)

DprojCG: $(common_obj) $(DprojCG_obj)
	$(CPP) -o $@ $^ $(LDFLAGS)

chist: $(common_obj) $(chist_obj)
	$(CPP) -o $@ $^ $(LDFLAGS)

-include $(common_dep) # include all dep files in makefile
-include $(InFlux_dep)
-include $(InFluxCG_dep)
-include $(InFluxNum_dep)
-include $(Dproj_dep) 
-include $(Dproj2_dep) 
-include $(DprojCG_dep) 
-include $(chist_dep) 

%.d: %.cpp 	# rule to generate a dep file by using the g++ prepocesor
	$(CPP) $(CPPFLAGS) -MM -MT $(@:.d=.o) $< -MF $@

%.o: %.c
	$(CPP) $(CPPFLAGS) -o $@ -c $<

.PHONY: debug InFlux InFluxNum Dproj chist
debug: CPPFLAGS+= -DDEBUG -ggdb # debug with gdb
debug: InFlux InFluxNum Dproj Dproj2 DprojCG chist InFluxCG

.PHONY: clean
clean: 
	$(RM) $(common_obj) $(InFlux_obj) $(InFluxNum_obj) $(Dproj_obj) $(chist_obj) *.d *~ *# InFlux InFLuxCG InfluxNum Dproj Dproj2 DprojCG chist

.PHONY: install
exec = InFlux InFluxNum Dproj Dproj2 DprojCG chist InFluxCG
install: $(exec)
	install $(exec) $(DESTDIR)
