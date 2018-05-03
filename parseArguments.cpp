#include <iostream>
#include <exception>
#include <boost/program_options.hpp>
#include "parseArguments.h"

namespace po = boost::program_options;

// Parse command-line arguments and store the values into Parameters.
// Return 1 if displaying help, 2 if a problem occurs, 0 otherwise.
int parseArguments(int argc, char **argv, Parameters &p) {
	po::options_description opts("Options");
	opts.add_options()
		("sites,n", po::value<long>(&p.nbSites)->required(), "Number of sites")
		("particles,m", po::value<long>(&p.nbParticles)->required(),
		 "Number of particles")
		("duration,T", po::value<double>(&p.duration)->required(),
		 "Duration of the simulation")
		("dt,t", po::value<double>(&p.dt)->required(),
		 "Timestep for export")
		("simuls,s", po::value<long>(&p.nbSimuls)->required(),
		 "Number of repetitions of the simulation")
		("initPos,d",
		 po::value< std::vector<long> >(&p.initPos)->multitoken()->required(),
		 "Initial positions of the tracers. It fixes the number of tracers.")
		("probas,p",
		 po::value< std::vector<double> >(&p.probas)->multitoken()->required(),
		 "Probabilities to jump to the right.")
		("moments,M", po::value<int>(&p.nbMoments)->default_value(
			DEFAULT_NB_MOMENTS), "Number of moments to compute")
		("threads,c", po::value<int>(&p.nbThreads)->default_value(
			DEFAULT_THREADS), "Number of threads")
		("output,o", po::value<std::string>(&p.output)->default_value(
			DEFAULT_OUTPUT_FILE), "Output file")
        ("alt,A", po::bool_switch(&p.alt),
		 "Use alternative algorithm for high density")
        ("single,S", po::bool_switch(&p.computeObs1TP),
		 "Compute and export occupations for a single particle")
        ("occ,O", po::bool_switch(&p.computeOcc),
		 "Compute and export occupations at final time")
        ("opn", po::bool_switch(&p.computeOccPosNeg),
		 "Compute and export occupations for sites +1 and -1")
        ("var", po::bool_switch(&p.computeVars),
		 "Compute and export variances of the displacements")
        ("W", po::bool_switch(&p.computeW),
		 "Compute and export W")
        ("visu,V", po::bool_switch(&p.visu), "Visualisation")
        ("sleep", po::value<int>(&p.sleep)->default_value(DEFAULT_VISU_SLEEP),
		 "Number of ms between visualizations")
        ("verbose,v", po::bool_switch(&p.verbose), "Verbose mode")
        ("help,h", "Print help message and exit")
		;

	try {
		po::variables_map vars;
		po::store(po::parse_command_line(argc, argv, opts), vars);

        // Display help and exit
        if (vars.count("help")) {
			std::cout << "Usage: " << argv[0] << " options\n";
			std::cout << opts << std::endl;
            return 1;
        }

        po::notify(vars);

		// The number of tracers is given by the number of positions
		p.nbTracers = (long) p.initPos.size();
		// Number of timesteps
		p.nbSteps = (long) (p.duration / p.dt);
	} catch (std::exception &e) {
		std::cerr << "Error: " << e.what() << std::endl;
		return 2;
	}

	return 0;
}

