#include <iostream>
#include <string>
#include <vector>

// from the shared library
#include <test/tests.h>

TEST("Example test", example_test) {
	if (false) test.fail("it failed");
}

struct TestObject {
	TEST_METHOD("Example method-test", example_test) {
		if (false) test.fail("it failed");
	}
};
TestObject myObject; // Can be customised/templated

TEST("Example benchmark", example_benchmarks) {
	std::vector<int> configs = {10, 100, 1000};

	std::cout << "size:\t";
	BenchmarkRate::print(configs);

	std::cout << "rate:\t";
	std::vector<double> rates = BenchmarkRate::map<int>(configs, [](int config, int repeats, Timer &timer) {
		timer.scaleRate(config/1e6); // scale to mega-ops/second
		timer.start();

		for (int repeat = 0; repeat < repeats; ++repeat) {
			// Actual test code
			int sum = 0;
			for (int i = 0; i < config; ++i) {
				sum++;
			}
		}

		timer.stop();
	}, true);

	return test.pass();
}