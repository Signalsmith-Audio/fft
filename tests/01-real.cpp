#include <iostream>
#include <vector>
#include <cmath>
#include <complex>

#include "tests-common.h"

#define LOG_VALUE(expr) \
	(std::cout << #expr << " = " << (expr) << "\n")

#define FAIL_VALUE_PAIR(expr1, expr2) \
	( \
		LOG_VALUE(expr1), LOG_VALUE(expr2), \
		test.fail(#expr1 " and " #expr2) \
	)

template<bool modified=false>
void test_real(Test &test) {
	using signalsmith::FFT;
	using signalsmith::RealFFT;
	using signalsmith::ModifiedRealFFT;
	using std::vector;
	using std::complex;

	for (int size = 2; size < 100; size += 2) {
		vector<complex<double>> complexInput(size);
		vector<complex<double>> complexMid(size);
		vector<complex<double>> complexOutput(size);
		vector<double> realInput(size);
		vector<complex<double>> realMid(size); // Only need half, but check it's undisturbed
		vector<double> realOutput(size);
		
		FFT<double> fft(size);
		typename std::conditional<modified, ModifiedRealFFT<double>, RealFFT<double>>::type realFft(size);

		// Random inputs
		for (int i = 0; i < size; ++i) {
			double v = rand()/(double)RAND_MAX - 0.5;
			complexInput[i] = v;
			realInput[i] = v;
		}
		if (modified) {
			for (int i = 0; i < size; ++i) {
				double rotPhase = -M_PI*i/size;
				complex<double> rot = {cos(rotPhase), sin(rotPhase)};
				complexInput[i] *= rot;
			}
		}
		for (int i = size/2; i < size; ++i) {
			// Should be undisturbed - fill with known value
			realMid[i] = complex<double>{52, 21};
		}
		
		fft.fft(complexInput, complexMid);
		realFft.fft(realInput, realMid);
		
		// Check complex spectrum matches
		if (!modified) {
			if (complexMid[0].imag() > 1e-6) return test.fail("complexMid[0].imag()");
			if (abs(complexMid[0].real() - realMid[0].real()) > 1e-6) return FAIL_VALUE_PAIR(complexMid[0].real(), realMid[0].real());
			if (abs(complexMid[size/2].real() - realMid[0].imag()) > 1e-6) return FAIL_VALUE_PAIR(complexMid[size/2].real(), realMid[0].imag());
		}
		for (int i = modified ? 0 : 1; i < size/2; ++i) {
			complex<double> diff = complexMid[i] - realMid[i];
			if (abs(diff) > size*1e-6) {
				LOG_VALUE(i);
				return FAIL_VALUE_PAIR(complexMid[i], realMid[i]);
			}
		}
		for (int i = size/2; i < size; ++i) {
			// It should have left the second half of realMid completely alone
			if (realMid[i] != complex<double>{52, 21}) return test.fail("realMid second half");
		}
		
		fft.ifft(complexMid, complexOutput);
		realFft.ifft(realMid, realOutput);

		if (modified) {
			for (int i = 0; i < size; ++i) {
				double rotPhase = M_PI*i/size;
				complex<double> rot = {cos(rotPhase), sin(rotPhase)};
				complexOutput[i] *= rot;
			}
		}

		for (int i = 0; i < size; ++i) {
			if (complexOutput[i].imag() > size*1e-6) return test.fail("complexOutput[i].imag");

			if (abs(complexOutput[i].real() - realOutput[i]) > size*1e-6) {
				LOG_VALUE(size);
				LOG_VALUE(i);
				return FAIL_VALUE_PAIR(complexOutput[i], realOutput[i]);
			}
		}
	}
}

TEST("Random real", random_real) {
	test_real<false>(test);
}

TEST("Modified real", random_modified_real) {
	test_real<true>(test);
}
