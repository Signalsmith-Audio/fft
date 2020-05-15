#include <stdio.h>
#include <fstream>
#include <iostream>

#include "../tests/tests-common.h"

template<typename Implementation, int constMaxSize=INT_MAX>
struct Benchmark {
	static std::vector<int> getSizes() {
		std::vector<int> possibleSizes(0);
		int size = 1;
		int maxSize = std::min(constMaxSize, (1<<16));
		while (size <= maxSize) {
			possibleSizes.push_back(size);
			size *= 2;
		}

		std::vector<int> result;
		for (int size : possibleSizes) {
			if (size <= maxSize) {
				result.push_back(size);
			}

			std::vector<int> mults = {3, 5, 7, 11, 15, 23};
			for (int mult : mults) {
				if (size*mult < maxSize) {
					result.push_back(size*mult);
				}
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