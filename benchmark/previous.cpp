#include <string>

#include "benchmark.h"

#include "previous/signalsmith-fft.h"

using signalsmith::fft::FFT;
using signalsmith::fft::getFft;

template <bool permute>
struct SignalsmithPrevious {

	static std::string name() {
		return std::string("Previous") + (permute ? " permute" : "");
	}
	static std::string resultTag() {
		return std::string("previous")+ (permute ? "-permute" : "");
	}
	static std::string version() {
		return "1";
	}

	template <typename T>
	struct OutOfPlace : OutOfPlaceRunner<T> {
		std::shared_ptr<FFT<double>> fft;

		OutOfPlace(size_t targetSize) : OutOfPlaceRunner<T>(targetSize), fft(getFft<double>(targetSize)) {
			size_t size = targetSize;
			while (fft->size() > targetSize) {
				--size;
				fft->setSize(size);
			}
		}

		void forward(std::complex<T> *input, std::complex<T> *output) {
			fft->fft((double*)input, (double*)output);
			if (permute) fft->permute((double*)output);
		}
	};

	using DoubleOutOfPlace = OutOfPlace<double>;
};
Benchmark<SignalsmithPrevious<true>> bench1;
Benchmark<SignalsmithPrevious<false>> bench2;
