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
	if (p.nbTracers > p.nbParticles) {
		std::cerr << "Error: The number of tracers should be smaller than"
			<< " the number of particles." << std::endl;
		return 1;
	}
	if ((long) p.probas.size() != p.nbTracers) {
		std::cerr << "Critical error: the size of the vector of probabilies"
			<< " is wrong." << std::endl;
		return 1;
	}
	if ((long) p.dists.size() != p.nbTracers - 1) {
		std::cerr << "Error: the number of distances should be the number"
			<< " of tracers minus one." << std::endl;
		return 1;
	}
	for (auto pr : p.probas) {
		if (pr < 0. || pr > 1.) {
			std::cerr << "Error: the probabilities should be between 0"
				<< " and 1." << std::endl;
			return 1;
		}
	}
	long sumL = 0;
	for (auto d : p.dists) {
		if (d <= 0) {
			std::cerr << "Error: The distances should be strictly positive."
				<< std::endl;
			return 1;
		}
		sumL += d;
	}
	if (sumL >= p.nbSites) {
		std::cerr << "Error: The sum of the distances should be inferior to"
			<< " the size of the system." << std::endl;
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
		<< ", tracers=" << p.nbTracers << ", iters=" << p.nbIters
		<< ", simuls=" << p.nbSimuls << ", moments=" << p.nbMoments;
	stream << ", probas=";
	for (auto pr : p.probas) {
		stream << pr << ":";
	}
	if (p.nbTracers > 1) {
		stream << ", dists=";
		for (auto d : p.dists) {
			stream << d << ":";
		}
	}
}

