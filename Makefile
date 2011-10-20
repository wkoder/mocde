#CXXFLAGS =	-O3 -g -Wall -fmessage-length=0
CXXFLAGS = -g -Wall -fmessage-length=0
CXX = g++ $(CXXFLAGS)
OBJS =	randomlib.o benchmark.o rand.o util.o mocde.o problemdef.o
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
	$(CXX) -c mocde.cpp

$(EXE).o:
	$(CXX) -c $(EXE).cpp
$(EXE): $(EXE).o $(OBJS)
	$(CXX) -o $(EXE) $(EXE).o $(OBJS) $(LIBS)
	
all:	clean $(EXE)

clean:
	rm -f *.o $(EXE)
