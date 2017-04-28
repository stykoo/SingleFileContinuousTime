#ifndef PARAMETERS_H
#define PARAMETERS_H

#define DEFAULT_PROBA_RIGHT 0.5
#define DEFAULT_NB_MOMENTS 4
#define DEFAULT_OUTPUT_PRECISION 15
#define DEFAULT_OUTPUT_FILE "observables.dat"
#define DEFAULT_THREADS 1
#define DEFAULT_VISU_SLEEP 200

#include <iostream>
#include <string>

struct Parameters {
	long nbSites;  // Number of sites
	long nbParticles;  // Number of vacancies
	long nbIters;  // Number of iterations
	long nbSimuls;  // Number of simulations
	int nbMoments;  // Number of moments to compute
	int nbThreads;  // Number of threads
	std::string output;
	
	bool alt;  // Use alternative algorithm for high density
	bool visu;  // Visualization
	int sleep;  // Number of milliseconds between two visualizations
	bool verbose;  // Verbose mode
};

int checkParameters(const Parameters &p);
void printParameters(const Parameters &p, std::ostream &stream = std::cout);

#endif
