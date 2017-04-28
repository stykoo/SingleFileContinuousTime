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
		("iters,i", po::value<long>(&p.nbIters)->required(),
		 "Number of time iterations")
		("simuls,s", po::value<long>(&p.nbSimuls)->required(),
		 "Number of repetitions of the simulation")
		("moments,M", po::value<int>(&p.nbMoments)->default_value(
			DEFAULT_NB_MOMENTS), "Number of moments to compute")
		("threads,t", po::value<int>(&p.nbThreads)->default_value(
			DEFAULT_THREADS), "Number of threads")
		("output,o", po::value<std::string>(&p.output)->default_value(
			DEFAULT_OUTPUT_FILE), "Output file")
        ("alt,A", po::bool_switch(&p.alt),
		 "Use alternative algorithm for high density")
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
	} catch (std::exception &e) {
		std::cerr << "Error: " << e.what() << std::endl;
		return 2;
	}

	return 0;
}

