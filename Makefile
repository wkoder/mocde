CXXFLAGS =	-O3 -g -Wall -fmessage-length=0
#CXXFLAGS = -g -Wall -fmessage-length=0
CXX = g++ $(CXXFLAGS)

OBJS =	randomlib.o benchmark.o compactdifferentialevolution.o rand.o util.o

LIBS =

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

cde.o:
	$(CXX) -c cde.cpp
cde: cde.o $(OBJS)
	$(CXX) -o cde cde.o $(OBJS) $(LIBS)
	
all:	clean cde

clean:
	rm -f *.o cde
