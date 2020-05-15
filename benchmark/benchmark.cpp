#include <string>

#include "benchmark.h"

struct Signalsmith {
	static std::string name() {
		return "Signalsmith";
	}
	static std::string resultTag() {
		return "signalsmith";
	}
	static std::string version() {
		return "(development)";
	}

	template <typename T>
	struct OutOfPlace {
		signalsmith::FFT<T> fft;
		OutOfPlace(size_t size) : fft(size) {}

		void forward(std::complex<T> *input, std::complex<T> *output) {
			fft.fft(input, output);
		}
	};

	using DoubleOutOfPlace = OutOfPlace<double>;
};
Benchmark<Signalsmith> benchSignalsmith;
