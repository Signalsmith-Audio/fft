#include <string>

#include "benchmark.h"

#include "previous/fft-v4.h"

struct Signalsmith {
	static std::string name() {
		return "previous-v4";
	}
	static std::string resultTag() {
		return "previous-v4";
	}
	static std::string version() {
		return "(development)";
	}

	template <typename T>
	struct OutOfPlace : public OutOfPlaceRunner<T> {
		signalsmith_v4::FFT<T> fft;
		OutOfPlace(size_t size) : OutOfPlaceRunner<T>(size), fft(size) {}

		void forward(std::complex<T> *input, std::complex<T> *output) {
			fft.fft(input, output);
		}
	};

	using DoubleOutOfPlace = OutOfPlace<double>;
};
Benchmark<Signalsmith> benchSignalsmith;
