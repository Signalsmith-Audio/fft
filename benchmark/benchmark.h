#include <stdio.h>
#include <fstream>
#include <iostream>
#include <algorithm>

#include "../tests/tests-common.h"

template<typename T>
class OutOfPlaceRunner {
	std::vector<std::complex<T>> inVector, outVector;
public:
	size_t size;

	OutOfPlaceRunner(size_t size) : size(size) {}
	virtual ~OutOfPlaceRunner() {}
	
	virtual void getPointers(std::complex<T> **inPointer, std::complex<T> **outPointer) {
		inVector.resize(size);
		outVector.resize(size);
		*inPointer = inVector.data();
		*outPointer = outVector.data();
	}
};

template<typename Implementation, int constMaxSize=INT_MAX>
struct Benchmark {
	static std::vector<int> getSizes(int customMaxSize=16777216) {
		std::vector<int> possibleSizes(0);
		int size = 1;
		int maxSize = std::min(constMaxSize, customMaxSize);
		while (size <= maxSize) {
			possibleSizes.push_back(size);
			size *= 2;
		}

		std::vector<int> result;
		for (int size : possibleSizes) {
			if (size <= maxSize) {
				result.push_back(size);
			}

			std::vector<int> mults = {3, 5, 7, 9, 11, 15, 23};
			mults = {3, 5, 9, 15, 25};
			mults = {3, 9};
			for (int mult : mults) {
				if (size*mult < maxSize) {
					result.push_back(size*mult);
				}
			}
		}
		return result;
	}

	template <typename Runner>
	static void runBenchmark(Test &test, std::string name, std::string resultPrefix) {
		auto sizes = getSizes();
		std::sort(sizes.begin(), sizes.end());

		std::cout << name << ":\n";
		std::vector<double> rates = BenchmarkRate::map<int>(sizes, [](int size, int repeats, Timer &timer) {
			auto runner = Runner(size);

			// The runner handles input/output allocation
			std::complex<double> *input, *output;
			runner.getPointers(&input, &output);

			for (int i = 0; i < size; i++) {
				input[i] = randomComplex<double>();
			}

			timer.start();
			for (int repeat = 0; repeat < repeats; ++repeat) {
				runner.forward(input, output);
			}
			timer.stop();
		}, true);
		std::vector<double> scaled = rates;
		for (size_t i = 0; i < sizes.size(); i++) {
			double size = sizes[i];
			double scaling = std::max(1.0, size*log(size))*1e-6;
			scaled[i] = rates[i]*scaling;
		}

		std::ofstream outputJs, outputCsv;
		outputJs.open(resultPrefix + ".js");
		outputCsv.open(resultPrefix + ".csv");

		outputJs << "addResults(\"" << name << "\", [";
		outputCsv << "size,ops/sec," << name << "\n";
		outputCsv.precision(15);

		for (size_t i = 0; i < sizes.size(); i++) {
			if (i > 0) outputJs << ",";
			outputJs << "\n\t{size: " << sizes[i] << ", rate: " << scaled[i] << "}";
			outputCsv << sizes[i] << "," << rates[i] << "," << scaled[i] << "\n";
		}
		outputJs << "\n]);";

		return test.pass();
	}

	TEST_METHOD("double out-of-place", double_out) {
		std::string resultPrefix = std::string("results/") + Implementation::resultTag();
		runBenchmark<typename Implementation::DoubleOutOfPlace>(test, Implementation::name(), resultPrefix);
	}
};
