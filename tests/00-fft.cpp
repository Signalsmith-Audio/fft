#include <iostream>
#include <vector>
#include <cmath>
#include <complex>

#include "tests-common.h"

TEST("2^N bins", test_2N_bins) {
	using signalsmith::FFT;
	using std::vector;
	using std::complex;

	std::vector<int> sizes = {1, 2, 4, 8, 16, 32, 64, 128, 256};
	
	for (int size : sizes) {
		vector<complex<double>> input(size);
		vector<complex<double>> inputCopy(size);
		vector<complex<double>> output(size);
		vector<complex<double>> expected(size);

		FFT<double> fft(size);

		// Test each bin
		for (int bin = 0; bin < size; bin++) {
			for (int i = 0; i < size; i++) {
				double phase = 2*M_PI*i*bin/size;
				input[i] = inputCopy[i] = complex<double>{cos(phase), sin(phase)};
				expected[i] = (i == bin) ? size : 0;
			}

			fft.fft(input, output);

			if (!closeEnough(input, inputCopy)) {
				return test.fail("input was changed");
			}
			if (!closeEnough(output, expected)) {
				return test.fail("output != expected");
			}
		}
	}
}

template<typename T>
std::complex<T> randomComplex() {
	std::complex<T> r;
	r.real(rand()/(double)RAND_MAX - 0.5);
	r.imag(rand()/(double)RAND_MAX - 0.5);
	return r;
}

TEST("2^N linearity", test_2N_linearity) {
	using signalsmith::FFT;
	using std::vector;
	using std::complex;

	std::vector<int> sizes = {1, 2, 4, 8, 16, 32, 64, 128, 256};
	
	for (int size : sizes) {
		vector<complex<double>> inputA(size);
		vector<complex<double>> inputB(size);
		vector<complex<double>> inputAB(size);
		vector<complex<double>> outputA(size);
		vector<complex<double>> outputB(size);
		vector<complex<double>> outputAB(size);
		vector<complex<double>> outputSummed(size);

		FFT<double> fft(size);

		// Test linearity
		for (int i = 0; i < size; i++) {
			inputA[i] = randomComplex<double>();
			inputB[i] = randomComplex<double>();
			inputAB[i] = inputA[i] + inputB[i];
		}

		fft.fft(inputA, outputA);
		fft.fft(inputB, outputB);
		fft.fft(inputAB, outputAB);

		for (int i = 0; i < size; i++) {
			outputSummed[i] = outputA[i] + outputB[i];
		}

		if (!closeEnough(outputAB, outputSummed)) {
			return test.fail("result was not linear");
		}
	}
}
