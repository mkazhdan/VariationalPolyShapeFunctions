TARGET= DeformableSolids3D
SOURCE= DeformableSolids3D.cpp

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

ifeq ($(OS),Windows_NT)
    detected_OS := Windows
else
    detected_OS := $(shell uname)
endif

CFLAGS += -fpermissive -Wno-deprecated -Wno-unused-result -Wno-format -msse2 -std=c++17 -O3 -funroll-loops -ffast-math

ifeq ($(detected_OS),Darwin)
		LFLAGS = -lgomp -lz -framework GLUT -framework OpenGL -lGLEW -lpng -ljpeg -O3 
else
        LFLAGS = -lgomp -lz -lGL -lGLU -lglut -lGLEW -lpng -ljpeg -O3 
endif

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


