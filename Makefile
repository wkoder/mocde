CXX = g++
CXXFLAGS = -O3 -g -Wall -fmessage-length=0
#CXXFLAGS = -g -Wall -fmessage-length=0
SRC = src/
BIN = bin/
VPATH = $(SRC):$(SRC)moead:$(SRC)paes:$(SRC)nsga2:$(SRC)jmetal
OBJS = $(patsubst $(SRC)%.cpp, $(BIN)%.o, $(wildcard $(SRC)*.cpp)) \
		$(patsubst $(SRC)%.c, $(BIN)%.o, $(wildcard $(SRC)*.c)) \
		$(patsubst $(SRC)paes/%.cpp, $(BIN)%.o, $(wildcard $(SRC)paes/*.cpp)) \
		$(patsubst $(SRC)nsga2/%.cpp, $(BIN)%.o, $(wildcard $(SRC)nsga2/*.cpp)) \
		$(patsubst $(SRC)nsga2/%.c, $(BIN)%.o, $(wildcard $(SRC)nsga2/*.c)) \
		$(patsubst $(SRC)jmetal/%.cpp, $(BIN)%.o, $(wildcard $(SRC)jmetal/*.cpp))
DEPS = $(OBJS) $(patsubst $(SRC)moead/%.h, $(BIN)%.h.gch, $(wildcard $(SRC)moead/*.h))
EXE = mocderunner

init:
	mkdir -p $(BIN)
compile: $(DEPS)
	

$(BIN)%.h.gch: $(SRC)moead/%.h
	$(CXX) $(CXXFLAGS) -c $< -o $@
$(BIN)%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@
$(BIN)%.o: %.c
	$(CXX) $(CXXFLAGS) -c $< -o $@
$(EXE): init
	make -j compile
	$(CXX) $(CXXFLAGS) $(OBJS) -o $(BIN)$(EXE)
	
all:	clean $(EXE)

clean:
	rm -f $(BIN)*.o $(BIN)*.h.gch $(BIN)$(EXE)
