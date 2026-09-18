CXX := g++
CXXFLAGS := -std=c++17 -O2 -Wall
SRC := src/rcwa_convergence.cpp
BIN := rcwa_convergence

.PHONY: all run clean

all: $(BIN)

$(BIN): $(SRC)
	$(CXX) $(CXXFLAGS) -o $(BIN) $(SRC)

run: $(BIN)
	./$(BIN)

clean:
	rm -f $(BIN)
