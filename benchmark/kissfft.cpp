#include <string>

#include "benchmark.h"

#include "others/kissfft/kissfft.hh"

struct Kiss {
	static std::string name() {
		return "KissFFT";
	}
	static std::string resultTag() {
		return "kissfft";
	}
	static std::string version() {
		return "(?)";
	}

	template <typename T>
	struct OutOfPlace : public OutOfPlaceRunner<T> {
		kissfft<T> fftForward;

		OutOfPlace(size_t size) : OutOfPlaceRunner<T>(size), fftForward(size, false) {}

		void forward(std::complex<T> *input, std::complex<T> *output) {
			fftForward.transform(input, output);
		}
	};

	using DoubleOutOfPlace = OutOfPlace<double>;
};
Benchmark<Kiss> benchKiss;
