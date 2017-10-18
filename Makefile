CPP=g++ 
CPPFLAGS=-Wall -O2 -fopenmp

LDFLAGS=-lnetcdf_c++ -fopenmp
RM=rm -rf

DESTDIR=$(HOME)/bin

common_src= Constants.cpp  EqDate.cpp  GridConstruction.cpp LagrangianEngine.cpp  ParseParameters.cpp  VectorXYZ.cpp  VFlow.cpp  VTKdump.cpp
common_obj=$(common_src:.cpp=.o) 
common_dep=$(common_obj:.o=.d)  # one dependency file for each source

#InFlux
InFlux_src=InFlux.cpp
InFlux_obj=$(InFlux_src:.cpp=.o) 
InFlux_dep=$(InFlux_obj:.o=.d)  # one dependency file for each source

#InFluxNum
InFluxNum_src=InFluxNum.cpp
InFluxNum_obj=$(InFluxNum_src:.cpp=.o) 
InFluxNum_dep=$(InFluxNum_obj:.o=.d)  # one dependency file for each source

#Dproj
Dproj_src=Dproj.cpp
Dproj_obj=$(Dproj_src:.cpp=.o) 
Dproj_dep=$(Dproj_obj:.o=.d)  # one dependency file for each source

.PHONY: all InFlux InFluxNum Dproj

all: all InFlux InFluxNum Dproj


InFlux: $(common_obj) $(InFlux_obj)
	$(CPP) -o $@ $^ $(LDFLAGS)

InFluxNum: $(common_obj) $(InFluxNum_obj)
	$(CPP) -o $@ $^ $(LDFLAGS)

Dproj: $(common_obj) $(Dproj_obj)
	$(CPP) -o $@ $^ $(LDFLAGS)

-include $(common_dep) # include all dep files in makefile
-include $(InFlux_dep)
-include $(InFluxNum_dep) 
-include $(Dproj_dep) 

%.d: %.cpp 	# rule to generate a dep file by using the g++ prepocesor
	$(CPP) $(CPPFLAGS) -MM -MT $(@:.d=.o) $< -MF $@

%.o: %.c
	$(CPP) $(CPPFLAGS) -o $@ -c $<

.PHONY: debug InFlux InFluxNum Dproj
debug: CPPFLAGS+= -DDEBUG -ggdb # debug with gdb
debug: InFlux InFluxNum Dproj

.PHONY: clean
clean: 
	$(RM) $(common_obj) $(InFlux_obj) $(InFluxNum_obj) $(Dproj_obj) *.d *~ *# InFlux InfluxNum Dproj

.PHONY: install
install: InFlux InFluxNum Dproj
	install $< $(DESTDIR)
