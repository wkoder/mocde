CXXFLAGS = -O3 -g -Wall -fmessage-length=0
#CXXFLAGS = -g -Wall -fmessage-length=0
CXX = g++ $(CXXFLAGS)
SRC = src/
BIN = bin/
OBJS = $(BIN)randomlib.o $(BIN)benchmark.o $(BIN)rand.o $(BIN)util.o $(BIN)mocde.o $(BIN)problemdef.o $(BIN)paes.o
#DEPS = $(OBJS)
DEPS = $(OBJS) $(BIN)algorithm.h.gch $(BIN)common.h.gch $(BIN)global.h.gch $(BIN)individual.h.gch \
				$(BIN)random.h.gch $(BIN)recomb.h.gch $(BIN)cec09.h.gch
LIBS =
EXE = mocderunner

init:
	mkdir -p $(BIN)

$(BIN)randomlib.o:
	$(CXX) -c $(SRC)randomlib.c -o $(BIN)randomlib.o $(LIBS)
$(BIN)rand.o:
	$(CXX) -c $(SRC)rand.c -o $(BIN)rand.o $(LIBS)
$(BIN)util.o:
	$(CXX) -c $(SRC)util.cpp -o $(BIN)util.o $(LIBS)
$(BIN)problemdef.o:
	$(CXX) -c $(SRC)problemdef.cpp -o $(BIN)problemdef.o $(LIBS)
$(BIN)benchmark.o:
	$(CXX) -c $(SRC)benchmark.cpp -o $(BIN)benchmark.o $(LIBS)
$(BIN)mocde.o:
	$(CXX) -c $(SRC)mocde.cpp -o $(BIN)mocde.o $(LIBS)
$(BIN)paes.o:
	$(CXX) -c $(SRC)paes/paes.cpp -o $(BIN)paes.o $(LIBS)
		
#---> MOEA/D
$(BIN)algorithm.h.gch:
	$(CXX) -c $(SRC)moead/algorithm.h -o $(BIN)algorithm.h.gch $(LIBS)
$(BIN)common.h.gch:
	$(CXX) -c $(SRC)moead/common.h -o $(BIN)common.h.gch $(LIBS)
$(BIN)global.h.gch:
	$(CXX) -c $(SRC)moead/global.h -o $(BIN)global.h.gch $(LIBS)
$(BIN)individual.h.gch:
	$(CXX) -c $(SRC)moead/individual.h -o $(BIN)individual.h.gch $(LIBS)
$(BIN)random.h.gch:
	$(CXX) -c $(SRC)moead/random.h -o $(BIN)random.h.gch $(LIBS)
$(BIN)recomb.h.gch:
	$(CXX) -c $(SRC)moead/recomb.h -o $(BIN)recomb.h.gch $(LIBS)
$(BIN)cec09.h.gch:
	$(CXX) -c $(SRC)moead/cec09.h -o $(BIN)cec09.h.gch $(LIBS)
# <---

$(EXE).o:
	$(CXX) -c $(SRC)$(EXE).cpp -o $(BIN)$(EXE).o
$(EXE): init $(EXE).o $(DEPS)
	$(CXX) -o $(BIN)$(EXE) $(BIN)$(EXE).o $(OBJS) $(LIBS)
	
all:	clean $(EXE)

clean:
	rm -f $(BIN)*.o $(BIN)*.gch $(BIN)$(EXE)
