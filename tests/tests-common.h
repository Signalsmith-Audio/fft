#include <vector>
#include <complex>

// from the shared library
#include <test/tests.h>

#include "../signalsmith-fft.h"

template<typename T>
std::complex<T> randomComplex() {
	std::complex<T> r;
	r.real(rand()/(double)RAND_MAX - 0.5);
	r.imag(rand()/(double)RAND_MAX - 0.5);
	return r;
}

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

	if (!totalEnergy) return true;
	T errorRatio = sqrt(totalError/totalEnergy);
	return errorRatio < 1e-6;
}

template<typename T>
static void printArray(std::vector<T> array, bool newline=true) {
	for (unsigned int i = 0; i < array.size(); ++i) {
		if (i > 0) std::cout << "\t";
		std::cout << array[i];
	}
	if (newline) std::cout << std::endl;
}
