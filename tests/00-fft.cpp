#include <iostream>
#include <vector>
#include <cmath>
#include <complex>

#include "tests-common.h"

std::vector<int> testSizes() {
	return {
		1, 2, 4, 8, 16, 32, 64, 128, 256,
		3, 6, 9, 12, 18, 24,
		5, 10, 15, 20, 25,
		7, 14, 21, 28, 49,
		11, 13, 17, 19, 22, 23
	};
}


TEST("Individual bins", test_2N_bins) {
	using signalsmith::FFT;
	using std::vector;
	using std::complex;

	std::vector<int> sizes = testSizes();
	
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
				std::cout << "N = " << size << "\n";
				std::cout << "   input:\t";
				printArray(input);
				std::cout << "  output:\t";
				printArray(output);
				std::cout << "expected:\t";
				printArray(expected);

				return test.fail("output != expected");
			}
		}
	}
}

TEST("Linearity", test_2N_linearity) {
	using signalsmith::FFT;
	using std::vector;
	using std::complex;

	std::vector<int> sizes = testSizes();
	
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

template<int fixedHarmonics=-1>
void inverseTest(Test test) {
	using signalsmith::FFT;

	std::vector<int> sizes = testSizes();
	
	for (int size : sizes) {
		std::vector<std::complex<double>> input(size);
		std::vector<std::complex<double>> spectrum(size);
		std::vector<std::complex<double>> output(size);
		std::vector<std::complex<double>> expected(size);

		for (int i = 0; i < size; i++) {
			if (fixedHarmonics >= 0) {
				double freq = fixedHarmonics;
				double phase = 2*M_PI*i*freq/size;
				input[i] = {cos(phase), sin(phase)};
			} else {
				input[i] = randomComplex<double>();
			}
			expected[i] = input[i]*(double)size;
		}

		FFT<double> fft(size);
		fft.fft(input, spectrum);
		fft.ifft(spectrum, output);

		if (!closeEnough(output, expected)) {
			printArray(input);
			printArray(spectrum);
			std::cout << size << "\n";
			printArray(output);
			printArray(expected);
			return test.fail("inverse did not match");
		}
	}
}

TEST("Inverse (first harmonic)", inverse_f1) {
	return inverseTest<1>(test);
}

TEST("Inverse (random)", inverse_random) {
	return inverseTest<-1>(test);
}

struct Powers {
	size_t two = 0, three = 0, five = 0, remainder = 1;
};
Powers factorise(size_t size) {
	Powers powers;
	while (size%2 == 0) {
		size /= 2;
		powers.two++;
	}
	while (size%3 == 0) {
		size /= 3;
		powers.three++;
	}
	while (size%5 == 0) {
		size /= 5;
		powers.five++;
	}
	powers.remainder = size;
	return powers;
}

TEST("Sizes", sizes) {
	using signalsmith::FFT;

	for (size_t i = 1; i < 1000; ++i) {
		size_t above = FFT<double>::fastSizeAbove(i);
		size_t below = FFT<double>::fastSizeBelow(i);

		if (above < i) return test.fail("above < i");
		if (below > i) return test.fail("below > i");

		auto factorsAbove = factorise(above);
		auto factorsBelow = factorise(below);

		if (factorsAbove.remainder != 1) return test.fail("non-fast above remainder");
		if (factorsBelow.remainder != 1) return test.fail("non-fast below remainder");

		if (factorsAbove.three + factorsAbove.five > 2) return test.fail("above is too complex");
		if (factorsBelow.three + factorsBelow.five > 2) return test.fail("below is too complex");
	}
}
