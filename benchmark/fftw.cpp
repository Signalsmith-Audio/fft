#include <string>

#include "benchmark.h"

#include <fftw3.h>

template<int searchLevel>
struct FFTW_Estimate {
	static std::string name() {
		if (searchLevel == FFTW_ESTIMATE) {
			return "FFTW (estimate)";
		} else if (searchLevel == FFTW_MEASURE) {
			return "FFTW (measure)";
		}
		return "FFTW";
	}
	static std::string resultTag() {
		if (searchLevel == FFTW_ESTIMATE) {
			return "fftw-estimate";
		} else if (searchLevel == FFTW_MEASURE) {
			return "fftw-measure";
		}
		return "fftw";
	}
	static std::string version() {
		return "3.3.8";
	}

	template <typename T>
	struct OutOfPlace : public OutOfPlaceRunner<T> {
		fftw_plan fftForward;
		fftw_complex *input, *output;
		
		OutOfPlace(size_t size) : OutOfPlaceRunner<T>(size) {
			input = (fftw_complex*)fftw_malloc(size*sizeof(fftw_complex));
			output = (fftw_complex*)fftw_malloc(size*sizeof(fftw_complex));
			fftForward = fftw_plan_dft_1d(size, input, output, FFTW_FORWARD, searchLevel);
		}
		~OutOfPlace() {
			fftw_destroy_plan(fftForward);
			fftw_free(input);
			fftw_free(output);
		}

		virtual void getPointers(std::complex<T> **inPointer, std::complex<T> **outPointer) {
			*inPointer = (std::complex<T>*)input;
			*outPointer = (std::complex<T>*)output;
		}

		void forward(std::complex<T> *input, std::complex<T> *output) {
			fftw_execute(fftForward);
		}
	};

	using DoubleOutOfPlace = OutOfPlace<double>;
};
static Benchmark<FFTW_Estimate<FFTW_ESTIMATE>> benchEstimate;
static Benchmark<FFTW_Estimate<FFTW_MEASURE>, 65536> benchMeasure;
