#include <vector>
#include <complex>

// from the shared library
#include <test/tests.h>

#include "../fft-v4.h"

template<typename T>
bool closeEnough(std::vector<std::complex<T>> vectorA, std::vector<std::complex<T>> vectorB) {
	double totalEnergy = 0;
	double totalError = 0;
	if (vectorA.size() != vectorB.size()) return false;

	for (unsigned int i = 0; i < vectorA.size(); ++i) {
		T error = norm(vectorA[i] - vectorB[i]);
		T energy = norm(vectorA[i]*vectorA[i]) + norm(vectorB[i]*vectorB[i]);
		totalEnergy += energy;
		totalError += error;
	}

	T errorRatio = sqrt(totalError/totalEnergy);
	return errorRatio < 1e-6;
}