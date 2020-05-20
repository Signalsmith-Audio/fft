#include <string>

#include "benchmark.h"

#include "previous/signalsmith-fft.h"

using signalsmith::fft::FFT;
using signalsmith::fft::getFft;

struct SignalsmithPrevious {

	static std::string name() {
		return "Previous";
	}
	static std::string resultTag() {
		return "previous";
	}
	static std::string version() {
		return "1";
	}

	template <typename T>
	struct OutOfPlace {
		std::shared_ptr<FFT<double>> fft;

		OutOfPlace(size_t targetSize) : fft(getFft<double>(targetSize)) {
			size_t size = targetSize;
			while (fft->size() > targetSize) {
				--size;
				fft->setSize(size);
			}
		}

		void forward(std::complex<T> *input, std::complex<T> *output) {
			fft->fft((double*)input, (double*)output);
			fft->permute((double*)output);
		}
	};

	using DoubleOutOfPlace = OutOfPlace<double>;
};
Benchmark<SignalsmithPrevious> benchPrevious;