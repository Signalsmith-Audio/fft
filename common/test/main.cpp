#include "../simple-args.h"
#include "../console-colours.h"

#include "tests.h"
#include <cstdlib> // srand, rand
#include <ctime> // time

TestList _globalTestList;

void Test::run(int depth, bool silent) {
	if (running) return fail("Re-entered test function");
	if (!silent) {
		std::cerr << Console::Dim;
		for (int i = 0; i < depth - 1; i++) {
			std::cerr << "  >  ";
		}
		std::cerr << Console::Cyan << "Test: "
			<< Console::Reset << Console::Cyan << testName
			<< Console::White << " (" << codeLocation << ")" << Console::Reset << std::endl;
	}
	running = true;

	runFn(*this);

	running = false;
}

void TestList::add(Test& test) {
	if (currentlyRunning.size() > 0) {
		Test *latest = currentlyRunning[0];
		if (!latest->success) return;
		// This is a sub-test, run it immediately instead of adding it
		currentlyRunning.push_back(&test);
		test.run(currentlyRunning.size(), currentlySilent);
		currentlyRunning.pop_back();
		return;
	}
	tests.push_back(test);
}

void TestList::fail(std::string reason) {
	for (auto testPtr : currentlyRunning) {
		testPtr->fail(reason);
	}
}

int TestList::run(int repeats) {
	currentlySilent = false;
	for (int repeat = 0; repeat < repeats; repeat++) {
		for (unsigned int i = 0; i < tests.size(); i++) {
			Test& test = tests[i];
			currentlyRunning = {&test};
			test.run(0, currentlySilent);
			if (!test.success) {
				std::cerr << Console::Red << Console::Bright << "\nFailed: "
					<< Console::Reset << test.reason << "\n\n";
				return 1;
			}
		}
		currentlySilent = true;
	}
	currentlyRunning.resize(0);
	return 0;
}

double defaultBenchmarkTime = 1;
int defaultBenchmarkDivisions = 5;

int main(int argc, char* argv[]) {
	SimpleArgs args(argc, argv);
	args.helpFlag("help");

	int repeats = args.flag<int>("repeats", "loop the tests a certain number of times", 1);
	defaultBenchmarkTime = args.flag<double>("test-time", "target per-test duration for benchmarks (excluding setup)", 1);
	defaultBenchmarkDivisions = args.flag<double>("test-divisions", "target number of sub-divisions for benchmarks", 5);
	if (args.error()) return 1;

	int randomSeed = args.flag<int>("seed", "random seed", time(NULL));
	srand(randomSeed);
	std::cout << Console::Dim << "random seed: " << randomSeed << Console::Reset << "\n";
	return _globalTestList.run(repeats);
}
