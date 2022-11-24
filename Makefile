VERSION = 1


CC       = g++
CCFLAGS  = -ggdb -O2 -D_GLIBCXX_USE_CXX11_ABI=0 -std=c++11 -Wall
INCFLAGS = -Icommon/ -I.
LDFLAGS = -lm -pthread

SOURCES = Burer.cc BurerSolution.cc

OBJECTS = $(SOURCES:.cc=.o)

BIN = burer

all: burer

burer: main.cc $(SOURCES) $(HEADERS)
	$(CC) $(CCFLAGS) -o $(BIN) main.cc  $(SOURCES) $(INCFLAGS) $(LDFLAGS)

%.o: %.cc
	$(CC) $(CCFLAGS) -c $

