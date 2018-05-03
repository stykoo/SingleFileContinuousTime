#ifndef PARAMETERS_H
#define PARAMETERS_H

#define DEFAULT_PROBA_RIGHT 0.5
#define DEFAULT_NB_MOMENTS 2
#define DEFAULT_OUTPUT_PRECISION 15
#define DEFAULT_OUTPUT_FILE "observables.dat"
#define DEFAULT_THREADS 1
#define DEFAULT_VISU_SLEEP 200

#include <iostream>
#include <string>
#include <vector>

struct Parameters {
	long nbSites;  // Number of sites
	long nbParticles;  // Number of vacancies

	double duration;  // Duration of the simulation
	double dt;  // Timestep for export of observables
	long nbSteps;  // This is duration / dt

	long nbSimuls;  // Number of simulations
	int nbMoments;  // Number of moments to compute
	int nbThreads;  // Number of threads
	std::string output;  // Name of the output file

	long nbTracers;  // Number of tracers
	std::vector<long> initPos;  // Initial positions of the tracers
	std::vector<double> probas;  // Probabilities to jump to the right
	
	bool alt;  // Use alternative algorithm for high density
	bool computeObs1TP;  // Compute observables for a single TP
	bool computeOcc;  // Compute occupations
	bool computeOccPosNeg;  // Compute occupations for sites +1 and -1
	bool computeVars;  // Compute variance of TPs
	bool computeW;  // Compute W
	bool visu;  // Visualization
	int sleep;  // Number of milliseconds between two visualizations
	bool verbose;  // Verbose mode
};

int checkParameters(const Parameters &p);
void printParameters(const Parameters &p, std::ostream &stream = std::cout);

#endif
