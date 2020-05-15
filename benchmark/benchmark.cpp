#include <stdio.h>
#include <string>
#include <fstream>
#include <iostream>

#include "../tests/tests-common.h"

#include "dev-history/direct.h"
struct DevHistoryDirect {
	static std::string name() {
		return "Direct";
	}
	static std::string resultTag() {
		return "dev-history-direct";
	}
	static std::string version() {
		return "(initial)";
	}

	template <typename T>
	struct OutOfPlace {
		dev_direct::FFT<T> fft;
		OutOfPlace(size_t size) : fft(size) {}

		void forward(std::complex<T> *input, std::complex<T> *output) {
			fft.fft(input, output);
		}
	};

	using DoubleOutOfPlace = OutOfPlace<double>;
};

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

template<typename Implementation, int maxSize=INT_MAX>
struct Benchmark {
	static std::vector<int> getSizes() {
		std::vector<int> possibleSizes = {1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024};
		std::vector<int> result;
		for (int size : possibleSizes) {
			if (size <= maxSize) {
				result.push_back(size);
			}
		}
		return result;
	}

	template <typename Runner>
	static void runBenchmark(Test &test, std::string name, std::string resultFile) {
		auto sizes = getSizes();

		for (size_t i = 0; i < name.size(); i++) std::cout << "-";
		std::cout << "-\t";
		BenchmarkRate::print(sizes);
		std::cout << name << ":\n\t";
		std::vector<double> rates = BenchmarkRate::map<int>(sizes, [](int size, int repeats, Timer &timer) {
			timer.scaleRate(std::max(1.0, size*log(size))*1e-6);

			auto runner = Runner(size);

			std::vector<std::complex<double>> input(size);
			std::vector<std::complex<double>> output(size);
			for (int i = 0; i < size; i++) {
				input[i] = randomComplex<double>();
			}

			timer.start();
			for (int repeat = 0; repeat < repeats; ++repeat) {
				runner.forward(input.data(), output.data());
			}
			timer.stop();
		}, true);

		std::ofstream output;
		output.open(resultFile);

              output << "addResults(\"" << name << "\", [";
              for (size_t i = 0; i < sizes.size(); i++) {
              	if (i > 0) output << ",";
              	output << "\n\t{size: " << sizes[i] << ", rate: " << rates[i] << "}";
              }
              output << "\n]);";

		return test.pass();
	}

	TEST_METHOD("double out-of-place", double_out) {
		std::string resultFile = std::string("results/") + Implementation::resultTag() + ".js";
		runBenchmark<typename Implementation::DoubleOutOfPlace>(test, Implementation::name(), resultFile);
	}
};

Benchmark<Signalsmith> benchSignalsmith;
Benchmark<DevHistoryDirect, 64> benchDirect;
