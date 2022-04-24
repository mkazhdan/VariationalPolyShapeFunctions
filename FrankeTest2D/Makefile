TARGET= FrankeTest2D
SOURCE= FrankeTest2D.cpp

COMPILER ?= gcc
#COMPILER ?= clang

ifeq ($(COMPILER),gcc)
	CC=gcc
	CXX=g++
	CFLAGS=-fopenmp
else
	CC=clang
	CXX=clang++
endif
MD=mkdir

CFLAGS += -fpermissive -Wno-deprecated -Wno-unused-result -Wno-format -msse2 -std=c++17 -O3 -funroll-loops -ffast-math
LFLAGS = -lgomp -O3 

SRC = ./
BIN = ./../Bin/Linux/
BIN = ./../Bin/Linux/
BIN_O = ./
INCLUDE = /usr/include/ -I./../include/


OBJECTS=$(addprefix $(BIN_O), $(addsuffix .o, $(basename $(SOURCE))))

all: $(BIN)
all: $(BIN)$(TARGET)

clean:
	rm -f $(BIN)$(TARGET)
	rm -f $(OBJECTS)

$(BIN):
	$(MD) -p $(BIN)

$(BIN)$(TARGET): $(OBJECTS)
	$(CXX) -o $@ $(OBJECTS) $(LFLAGS)

$(BIN_O)%.o: $(SRC)%.c
	$(CC) -c -o $@ $(CFLAGS) -I$(INCLUDE) $<

$(BIN_O)%.o: $(SRC)%.cpp
	$(CXX) -c -o $@ $(CFLAGS) -I$(INCLUDE) $<


