#ifndef SIMUL_H
#define SIMUL_H

#include <vector>
#include <random>
#include "parameters.h"

#define NB_OBS_OCCUPATIONS 4

struct State {
	std::vector<long> positions;  // Positions of the particles
	std::vector<long> occupations;  // Occupations of the sites

	std::vector<long> freeParticles;  // Particles that are free to move
	std::vector<bool> freeOnLeft;  // Is particle free to move to the left?
	std::vector<bool> freeOnRight;  // Is particle free to move to the right?
	long nbFreeParticles;
};

struct Observables {
	std::vector< std::vector<long long> > moments;  // All the moments
	std::vector<long long> moments1TP;  // Moments of TP 1

	// Observables related to the occupations
	// 0 -> occupations (eta_l)
	// 1 -> occupation * occupation of site X+1 (eta_l * eta_1)
	// 2 -> X * eta_l
	// 3 -> X * eta_l * eta_1 
	std::vector< std::vector<long> > occObs;
};

int runSimulations(const Parameters &p);
void runMultipleSimulations(const Parameters &p, const long nbSimuls,
		                   std::vector<Observables> &sumObs,
						   const unsigned int seed);
void runOneSimulation(const Parameters &p, std::vector<Observables> &obs,
   					  std::mt19937 &rndGen);
void initState(State &state, const Parameters &p, std::mt19937 &rndGen);
void updateState(State &state, const Parameters &p, std::mt19937 &rndGen);
void initObservables(std::vector<Observables> &obs, const Parameters &p);
void computeObservables(const State &state, const Parameters &p,
		                Observables &o);
void computeOccObs(const State &state, const Parameters &p,
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
