CXXFLAGS =	-O3 -g -Wall -fmessage-length=0
#CXXFLAGS = -g -Wall -fmessage-length=0
CXX = g++ $(CXXFLAGS)
OBJS =	randomlib.o benchmark.o compactdifferentialevolution.o rand.o util.o
LIBS =
IMPL1 = mocde

randomlib.o:
	$(CXX) -c randomlib.c $(LIBS)
rand.o:
	$(CXX) -c rand.c $(LIBS)
util.o:
	$(CXX) -c util.cpp $(LIBS)
benchmark.o:
	$(CXX) -c benchmark.cpp $(LIBS)
compactdifferentialevolution.o:
	$(CXX) -c compactdifferentialevolution.cpp

$(IMPL1).o:
	$(CXX) -c $(IMPL1).cpp
$(IMPL1): $(IMPL1).o $(OBJS)
	$(CXX) -o $(IMPL1) $(IMPL1).o $(OBJS) $(LIBS)
	
all:	clean $(IMPL1)

clean:
	rm -f *.o $(IMPL1)
