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
	if (p.duration <= 0.) {
		std::cerr << "Error: duration should be strictly positive."
			<< std::endl;
		return 1;
	}
	if (p.dt <= 0.) {
		std::cerr << "Error: dt should be strictly positive."
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
	if ((long) p.initPos.size() != p.nbTracers) {
		std::cerr << "Critical error: the size of the vector of positions"
			<< " is wrong." << std::endl;
		return 1;
	}
	for (auto po : p.initPos) {
		if (po < 0 || po >= p.nbSites) {
			std::cerr << "Error: the positions should in the right range."
				<< std::endl;
			return 1;
		}
	}
	for (long i = 0 ; i < p.nbTracers - 1 ; ++i) {
		if (p.initPos[i] >= p.initPos[i+1]) {
			std::cerr << "Error: please sort the positions."
				<< std::endl;
			return 1;
		}
	}
	if ((long) p.probas.size() != p.nbTracers) {
		std::cerr << "Error: the number of probabilities and the number"
			<< " of tracers should be equal." << std::endl;
		return 1;
	}
	for (auto pr : p.probas) {
		if (pr < 0. || pr > 1.) {
			std::cerr << "Error: the probabilities should be between 0"
				<< " and 1." << std::endl;
			return 1;
		}
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
		<< ", tracers=" << p.nbTracers << ", duration=" << p.duration
		<< ", simuls=" << p.nbSimuls << ", moments=" << p.nbMoments;
	stream << ", initPos=";
	for (auto po : p.initPos) {
		stream << po << ":";
	}
	stream << ", probas=";
	for (auto pr : p.probas) {
		stream << pr << ":";
	}
}

