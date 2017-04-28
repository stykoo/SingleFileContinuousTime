#include "parameters.h"

// Check if the parameters are valid. Return 0 if they are, 1 otherwise.
int checkParameters(const Parameters &p) {
	if (p.nbSites <= 0) {
		std::cerr << "Error: nbSites should be strictly positive."
			<< std::endl;
		return 1;
	}
	if (p.nbParticles <= 0) {
		std::cerr << "Error: nbParticles should be strictly positive."
			<< std::endl;
		return 1;
	}
	if (p.nbIters <= 0) {
		std::cerr << "Error: nbIters should be strictly positive."
			<< std::endl;
		return 1;
	}
	if (p.nbSimuls <= 0) {
		std::cerr << "Error: nbSimuls should be strictly positive."
			<< std::endl;
		return 1;
	}
	if (p.nbMoments <= 0) {
		std::cerr << "Error: nbMoments should be strictly positive."
			<< std::endl;
		return 1;
	}
	if (p.alt && 2 * p.nbParticles <= p.nbSites) {
		std::cerr << "Warning: You shouldn't use the alternative algorithm "
			<< " at such a low density." << std::endl;
	} else if (!p.alt && 10 * p.nbParticles >= 9 * p.nbSites) {
		std::cerr << "Warning: You may want to use the alternative algorithm "
			<< " for high density." << std::endl;
	}
	return 0;
}

// Print the parameters to stream.
void printParameters(const Parameters &p, std::ostream &stream) {
	stream << "sites=" << p.nbSites << ", particles=" << p.nbParticles
		<< ", iters=" << p.nbIters << ", simuls=" << p.nbSimuls
		<< ", moments=" << p.nbMoments;
}

