#include <algorithm>
#include <fstream>
#include <map>
#include <chrono>
#include <thread>
#include <iomanip>
#include "simul.h"

// Run all the simulations and store the moments.
int runSimulations(const Parameters &p) {
	std::random_device rd;
	std::vector<std::thread> threads;
	std::vector< std::vector<Observables> > allSumsObs(p.nbThreads);

	long nbSimulsPerThread = p.nbSimuls / p.nbThreads;
	
	if (p.nbSimuls % p.nbThreads != 0) {
		std::cerr << "Warning: nbSimuls is not a multiple of nbThreads."
			<< nbSimulsPerThread * p.nbThreads << " simulations will be done."
			<< std::endl;
	}

	// Threads
	for (int i = 0 ; i < p.nbThreads ; ++i) {
		threads.push_back(std::thread(runMultipleSimulations, p,
	        nbSimulsPerThread, std::ref(allSumsObs[i]), rd())); 
	}
	
	// Wait for everyone
	for (auto &th : threads) {
		th.join();
	}

	// Initialize the total sum
	std::vector<Observables> sumObs;
	initObservables(sumObs, p);

	// Add the observables to the total sum
	for (int k = 0 ; k < p.nbThreads ; ++k) {
		addObservables(sumObs, allSumsObs[k], p);
	}

	int status = exportObservables(sumObs, p);

	return status;
}

// Run several simulation and store the sum of the observables.
// This function is usually called as a thread.
void runMultipleSimulations(const Parameters &p, const long nbSimuls,
		                   std::vector<Observables> &sumObs,
						   const unsigned int seed) {
    // Random generator
	std::mt19937 rndGen(seed);

	// Initialize sumObs
	initObservables(sumObs, p);

	for (long s = 0 ; s < nbSimuls ; ++s) {
		if (p.verbose || p.visu) {
			std::cout << "Simulation " << s + 1 
				<< " on thread " << std::this_thread::get_id() << std::endl;
		}

		// Run a single simulation
		std::vector<Observables> obs;
		runOneSimulation(p, obs, rndGen);

		// Add the observables to the sum
		addObservables(sumObs, obs, p);

		if (p.visu) {
			std::cout << std::endl;
		}
	}
}

// Run a single simulation and compute the moments.
void runOneSimulation(const Parameters &p, std::vector<Observables> &obs,
					 std::mt19937 &rndGen) {
	// To generate the time
	double t = 0., tLast = 0.;
	long tDiscrete = 0;
	std::exponential_distribution<double> rndTime((double) p.nbParticles);

	// Initialize obs
	initObservables(obs, p);

	// Positions of the tracers
	State state;
	initState(state, p, rndGen);

	// Compute the initial observables
	computeObservables(state, p, obs[tDiscrete]);
	++tDiscrete;

	if (p.visu) {
		visualize(state, p, 0.);
	}

	// Loop over time
	while (tDiscrete < p.nbSteps) {
		// Evolve
		updateState(state, p, rndGen);
		t += rndTime(rndGen);

		while (t - tLast > p.dt) {
			computeObservables(state, p, obs[tDiscrete]);
			++tDiscrete;
			tLast += p.dt;
		}
		if (p.visu) {
			visualize(state, p, t);
		}
	}

	// Observables for occupations at final time.
	if (p.computeOcc) {
		computeOccObs(state, p, obs[p.nbSteps-1]);
	}
}

// Generate an initial state.
void initState(State &state, const Parameters &p, std::mt19937 &rndGen) {
	state.positions.resize(p.nbParticles);
	state.occupations.assign(p.nbSites, -1);

	// Assign the sites of the tracers
	for (long i = 0 ; i < p.nbTracers ; ++i) {
		state.positions[i] = p.initPos[i];
		state.occupations[p.initPos[i]] = i;
	}

	// Generate a random permutation of the remaining sites
	// and assign the first sites to the particles.
	std::vector<long> seq(p.nbSites - p.nbTracers);
	long i = 0, j = 0;
	for (long k = 0 ; k < p.nbTracers ; ++k) {
		while (i < p.initPos[k]) {
			seq[j++] = i++;
		}
		++i;
	}
	while (i < p.nbSites) {
		seq[j++] = i++;
	}

	std::shuffle(seq.begin(), seq.end(), rndGen);
	for (long i=p.nbTracers ; i < p.nbParticles ; ++i) {
		state.positions[i] = seq[i - p.nbTracers];
		state.occupations[seq[i - p.nbTracers]] = i;
	}

	// Alternative algorithm at high density
	if (p.alt) {
		// There are at most twice as many free particles as vacancies
		state.freeParticles.assign(2 * (p.nbSites - p.nbParticles), -1);
		state.freeOnLeft.assign(p.nbParticles, false);
		state.freeOnRight.assign(p.nbParticles, false);

		long k = 0;
		// Check which particles are free on left / right
		for (long i=0 ; i < p.nbParticles ; ++i) {
			long pos = periodicSubs1(state.positions[i], p.nbSites);
			if(state.occupations[pos] == -1){
				state.freeOnLeft[i] = true;
			}
			pos = periodicAdd1(state.positions[i], p.nbSites);
			if(state.occupations[pos] == -1){
				state.freeOnRight[i] = true;
			}

			// Add the particle to the 'list' of free particles
			if(state.freeOnLeft[i] || state.freeOnRight[i]) {
				state.freeParticles[k] = i;
				++k;
			}
		}
		state.nbFreeParticles = k;
	}
}

// Implement one step of the time evolution of the system.
void updateState(State &state, const Parameters &p, std::mt19937 &rndGen) {
	std::uniform_int_distribution<long> distInt(0, p.nbParticles - 1);
	std::uniform_real_distribution<double> distReal(0.0, 1.0);

	long part = distInt(rndGen);  // Random particle
	long pos = state.positions[part]; 
	double u = distReal(rndGen);

	// Probability to jump to the right
	double pr = (part < p.nbTracers) ? p.probas[part] : DEFAULT_PROBA_RIGHT;

	if (u < pr) {
		long posR = periodicAdd1(pos, p.nbSites);
		if(state.occupations[posR] == -1){
			state.occupations[pos] = -1;
			state.occupations[posR] = part;
			state.positions[part] = posR;
		}
	} else {
		long posL = periodicSubs1(pos, p.nbSites);
		if(state.occupations[posL] == -1){
			state.occupations[pos] = -1;
			state.occupations[posL] = part;
			state.positions[part] = posL;
		}
	}
}

// Initialize a vector of observables
void initObservables(std::vector<Observables> &obs, const Parameters &p) {
	obs.resize(p.nbSteps);
	for (long t = 0 ; t < p.nbSteps ; ++t) {
		obs[t].moments.resize(p.nbTracers);
		for (long i = 0 ; i < p.nbTracers ; ++i) {
			obs[t].moments[i].assign(p.nbTracers - i, 0);
		}
	}

	// If we store the observables for a single TP
	if (p.computeObs1TP) {
		for (long t = 0 ; t < p.nbSteps ; ++t) {
			obs[t].moments1TP.assign(p.nbMoments, 0);
		}
	}
		
	// We allocate memory for occupations only for the final iteration
	if (p.computeOcc) {
		obs[p.nbSteps-1].occObs.resize(NB_OBS_OCCUPATIONS);
		for (int i = 0 ; i < NB_OBS_OCCUPATIONS ; ++i) {
			obs[p.nbSteps-1].occObs[i].assign(p.nbSites, 0);
		}
	}
}

// Compute the observables (except those for the occupations).
void computeObservables(const State &state, const Parameters &p,
		                Observables &o) {
	std::vector<long> xsPer(p.nbTracers);
	for (long i = 0 ; i < p.nbTracers ; ++i) {
		xsPer[i] = periodicBCsym(state.positions[i] - p.initPos[i], p.nbSites);
	}
	for (long i = 0 ; i < p.nbTracers ; ++i) {
		for (long j = 0 ; j < p.nbTracers - i ; ++j) {
			o.moments[i][j] = 1;
			for (long k = j ; k < j + i + 1 ; ++k) {
				o.moments[i][j] *= xsPer[k];
			}
		}
	}

	if (p.computeObs1TP) {
		for (long i = 0 ; i < p.nbMoments ; ++i) {
			o.moments1TP[i] = mypow(xsPer[0], i + 1); 
		}
	}
}

// Compute the observables related to the occupations.
// We assume the memory is already allocated.
void computeOccObs(const State &state, const Parameters &p,
		           Observables &o) {
	long xPer = periodicBCsym(state.positions[0] - p.initPos[0], p.nbSites);

	for (long i = 0 ; i < p.nbSites ; ++i) {
		long k = periodicBC(state.positions[0] + i, p.nbSites);
		long l = periodicAdd1(state.positions[0], p.nbSites);
		long eta = (long) (state.occupations[k] >= 0);
		long eta1 = (long) (state.occupations[l] >= 0);

		o.occObs[0][i] = eta;
		o.occObs[1][i] = eta * eta1;
		o.occObs[2][i] = xPer * eta;
		o.occObs[3][i] = xPer * eta * eta1;
	}
}

// Add observables o2 to observables o1.
void addObservables(std::vector<Observables> &obs1,
		            const std::vector<Observables> &obs2, const Parameters &p)
{
	for (long t = 0 ; t < p.nbSteps ; ++t) {
		for (long i = 0 ; i < p.nbTracers ; ++i) {
			for (long j = 0 ; j < p.nbTracers - i ; ++j) {
				obs1[t].moments[i][j] += obs2[t].moments[i][j];
			}
		}
	}

	if (p.computeObs1TP) {
		for (long t = 0 ; t < p.nbSteps ; ++t) {
			for (long i = 0 ; i < p.nbMoments ; ++i) {
				obs1[t].moments1TP[i] += obs2[t].moments1TP[i];
			}
		}
	}

	if (p.computeOcc) {
		for (long i = 0 ; i < p.nbSites ; ++i) {
			for (int j = 0 ; j < NB_OBS_OCCUPATIONS ; ++j) {
				obs1[p.nbSteps-1].occObs[j][i] +=
					obs2[p.nbSteps-1].occObs[j][i];
			}
		}
	}
}

// Export the observables to a file.
int exportObservables(const std::vector<Observables> &sumObs,
		              const Parameters &p) {
	std::ofstream file(p.output);
	if (!file.is_open()) {
		return 1;
	}

	// Header
	file << "# SingleFileContinuousTime (" << __DATE__ <<  ", " << __TIME__
		<< "): ";
	printParameters(p, file);
	file << "\n# t";
	for (long i = 0 ; i < p.nbTracers ; ++i) {
		for (long j = 0 ; j < p.nbTracers - i ; ++j) {
			file << " ";
			for (long k = j ; k < j + i + 1 ; ++k) {
				file << "x" << k+1;
			}
		}
	}
	file << "\n";

	file << std::scientific << std::setprecision(DEFAULT_OUTPUT_PRECISION);

	// Data (we write the average and not the sum)
	for (long k = 0 ; k < p.nbSteps ; ++k) {
		file << k * p.dt;
		for (long i = 0 ; i < p.nbTracers ; ++i) {
			for (long j = 0 ; j < p.nbTracers - i ; ++j) {
				file << " " << ((double) sumObs[k].moments[i][j]) / p.nbSimuls;
			}
		}
		file << "\n";
	}

	file.close();

	if (p.computeObs1TP) {
		size_t pos = p.output.find_last_of(".");
		std::string rawname = p.output.substr(0, pos); 
		std::string ext = p.output.substr(pos); 
		std::string name = rawname + "_1TP" + ext;

		file.open(name);
		if (!file.is_open()) {
			return 1;
		}

		// Header
		file << "# SingleFileContinuousTime (" << __DATE__ <<  ", " << __TIME__
			<< "): ";
		printParameters(p, file);
		file << "\n# t";
		for (long i = 0 ; i < p.nbMoments ; ++i) {
			file << " x1^" << i + 1;
		}
		file << "\n";

		file << std::scientific << std::setprecision(DEFAULT_OUTPUT_PRECISION);

		// Data (we write the average and not the sum)
		for (long k = 0 ; k < p.nbSteps ; ++k) {
			file << k * p.dt;
			for (long i = 0 ; i < p.nbMoments ; ++i) {
				file << " " << ((double) sumObs[k].moments1TP[i]) / p.nbSimuls;
			}
			file << "\n";
		}

		file.close();
	}

	// Export occupations (at final time) if needed
	if (p.computeOcc) {
		size_t pos = p.output.find_last_of(".");
		std::string rawname = p.output.substr(0, pos); 
		std::string ext = p.output.substr(pos); 
		std::string name = rawname + "_occ_" + std::to_string(p.duration)
			               + ext;

		file.open(name);
		if (!file.is_open()) {
			return 1;
		}

		file << "# SingleFileContinuousTime (" << __DATE__ <<  ", " << __TIME__
		   	<< "): ";
		printParameters(p, file);
		file << ", t=" + std::to_string(p.nbSteps*p.dt) << "\n";
		file << "# i eta_i eta_i*eta_1 X*eta_i X*eta_i*eta_1 X\n";

		double xfin = ((double) sumObs[p.nbSteps-1].moments[0][0])
			          / p.nbSimuls;

		for (long j = 0 ; j < p.nbSites ; ++j) {
			file << j;
			for (int k = 0 ; k < NB_OBS_OCCUPATIONS ; ++k) {
		   		file << " " << (((double) sumObs[p.nbSteps-1].occObs[k][j])
					            / p.nbSimuls);
			}
			file << " " << xfin << "\n";
		}

		file.close();
	}
	return 0;
}

// Print the system into terminal.
void visualize(const State &state, const Parameters &p, const double t) {
	std::map<long,char> mymap;
	std::map<long,char>::iterator it;
	for (long i = 0 ; i < p.nbTracers ; ++i) {
		mymap[state.positions[i]] = 'X';
	}
	for (long i = p.nbTracers ; i < p.nbParticles ; ++i) {
		mymap[state.positions[i]] = 'O';
	}

	std::cout << "\r";
	for (long i = 0 ; i < p.nbSites ; ++i) {
		it = mymap.find(i);
		if (it != mymap.end()) {
			std::cout << it->second;
		} else {
			std::cout << '.';
		}
	}
	std::cout << std::fixed << std::setprecision(5) <<  "  (t = " << t << ")"
		<< std::flush;
	
	// Sleep
	std::this_thread::sleep_for(std::chrono::milliseconds(p.sleep));
}


// Wrap x into [0,L).
long periodicBC(const long x, const long L) {
	return (x < 0) ? (x % L + L) : (x % L);
}

// Wrap x into (-L/2,L/2].
// TODO improve
long periodicBCsym(const long x, const long L) {
	long y = periodicBC(x, L);
	return (y > L / 2) ? (y - L) : y;
}

// Add 1 in periodic boundary conditions (assuming 0 <= x < L)
long periodicAdd1(const long x, const long L) {
	return (x + 1) % L;
}

// Substract 1 in periodic boundary conditions (assuming 0 <= x < L)
long periodicSubs1(const long x, const long L) {
	return (x - 1 + L) % L;
}

// Increment in periodic boundary conditions (assuming 0 <= x < L)
void periodicInc(long &x, const long L) {
	if (x < L - 1) {
		x++;
	} else {
		x = 0;
	}
}

// Decrement in periodic boundary conditions (assuming 0 <= x < L)
void periodicDec(long &x, const long L) {
	if (x > 0) {
		x--;
	} else {
		x = L - 1;
	}
}
