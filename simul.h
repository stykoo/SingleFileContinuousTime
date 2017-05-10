#ifndef SIMUL_H
#define SIMUL_H

#include <vector>
#include <random>
#include "parameters.h"

//typedef std::minstd_rand MyRng;
typedef std::mt19937 MyRng;

struct State {
	std::vector<long> positions;  // Positions of the particles
	std::vector<bool> occupations;  // Occupations of the sites
};

struct Observables {
	std::vector< std::vector<long long> > moments;  // All the moments
};

int runSimulations(const Parameters &p);
void runMultipleSimulations(const Parameters &p, const long nbSimuls,
		                   std::vector<Observables> &sumObs,
						   const unsigned int seed);
void runOneSimulation(const Parameters &p, std::vector<Observables> &obs,
   					  MyRng &rndGen);
void initState(State &state, const Parameters &p, MyRng &rndGen);
void updateState(State &state, const Parameters &p, MyRng &rndGen);
void initObservables(std::vector<Observables> &obs, const Parameters &p);
void computeObservables(const State &state, const Parameters &p,
		                Observables &o);
void addObservables(std::vector<Observables> &obs1,
		            const std::vector<Observables> &obs2, const Parameters &p);
int exportObservables(const std::vector<Observables> &sumObs,
		              const Parameters &p);
void visualize(const State &state, const Parameters &p, const double t);

long periodicBC(const long x, const long L);
long periodicBCsym(const long x, const long L);
long periodicAdd1(const long x, const long L);
long periodicSubs1(const long x, const long L);
void periodicInc(long &x, const long L);
void periodicDec(long &x, const long L);

// Compute a^b with b a positive integer
template<typename T, typename U>
T mypow(const T a, const U b) {
    if (b <= 0) {
        return 1;
	} else if (b % 2 == 1) {
        return a * mypow(a, b - 1);
	} else {
		T c = mypow(a, b / 2);
		return c * c;
	}
}

#endif
