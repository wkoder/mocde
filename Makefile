CXXFLAGS = -O3 -g -Wall -fmessage-length=0
#CXXFLAGS = -g -Wall -fmessage-length=0
CXX = g++ $(CXXFLAGS)
OBJS = randomlib.o benchmark.o rand.o util.o mocde.o problemdef.o
#DEPS = $(OBJS) algorithm.h.gch common.h.gch global.h.gch individual.h.gch random.h.gch recomb.h.gch cec09.h.gch
DEPS = $(OBJS)
LIBS =
EXE = mocderunner

randomlib.o:
	$(CXX) -c randomlib.c $(LIBS)
rand.o:
	$(CXX) -c rand.c $(LIBS)
util.o:
	$(CXX) -c util.cpp $(LIBS)
problemdef.o:
	$(CXX) -c problemdef.cpp $(LIBS)
benchmark.o:
	$(CXX) -c benchmark.cpp $(LIBS)
mocde.o:
	$(CXX) -c mocde.cpp $(LIBS)
	
#---> MOEA/D
#algorithm.h.gch:
#	$(CXX) -c moead/algorithm.h $(LIBS)
#common.h.gch:
#	$(CXX) -c moead/common.h $(LIBS)
#global.h.gch:
#	$(CXX) -c moead/global.h $(LIBS)
#individual.h.gch:
#	$(CXX) -c moead/individual.h $(LIBS)
#random.h.gch:
#	$(CXX) -c moead/random.h $(LIBS)
#recomb.h.gch:
#	$(CXX) -c moead/recomb.h $(LIBS)
#cec09.h.gch:
#	$(CXX) -c moead/cec09.h $(LIBS)
	
algorithm.h.gch:
	$(CXX) -c moead/algorithm.h -o algorithm.h.gch $(LIBS)
common.h.gch:
	$(CXX) -c moead/common.h -o common.h.gch $(LIBS)
global.h.gch:
	$(CXX) -c moead/global.h -o global.h.gch $(LIBS)
individual.h.gch:
	$(CXX) -c moead/individual.h -o individual.h.gch $(LIBS)
random.h.gch:
	$(CXX) -c moead/random.h -o random.h.gch $(LIBS)
recomb.h.gch:
	$(CXX) -c moead/recomb.h -o recomb.h.gch $(LIBS)
cec09.h.gch:
	$(CXX) -c moead/cec09.h -o cec09.h.gch $(LIBS)
# <---

$(EXE).o:
	$(CXX) -c $(EXE).cpp
$(EXE): $(EXE).o $(DEPS)
	$(CXX) -o $(EXE) $(EXE).o $(OBJS) $(LIBS)
	
all:	clean $(EXE)

clean:
	rm -f *.o *.gch moead/*.gch $(EXE)

### MOEAD Makefile ###
