CC=g++
CFLAGS=-W -Wall -ansi -pedantic -std=c++11 -O3
LDFLAGS=-lboost_program_options -lpthread
EXEC=singleFileContinuousTime

all: $(EXEC)

prof: CFLAGS += -pg
prof: LDFLAGS += -pg
prof: singleFileContinuousTime

singleFileContinuousTime: main.o simul.o parseArguments.o parameters.o
		$(CC) -o $@ $^ $(LDFLAGS)

simul_test: simul_test.cpp simul.o parameters.o
		$(CC) -o $@ $^ $(LDFLAGS)

main.o: main.cpp parseArguments.h parameters.h
		$(CC) -o $@ -c $< $(CFLAGS)

simul.o: simul.cpp simul.h parameters.h
		$(CC) -o $@ -c $< $(CFLAGS)

parseArguments.o: parseArguments.cpp parseArguments.h parameters.h
		$(CC) -o $@ -c $< $(CFLAGS)

parameters.o: parameters.cpp parameters.h
		$(CC) -o $@ -c $< $(CFLAGS)

.PHONY: clean mrproper

clean:
		rm -f *.o

mrproper: clean
		rm -rf $(EXEC)
