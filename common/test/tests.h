#ifndef _TEST_FRAMEWORK_TESTS_H
#define _TEST_FRAMEWORK_TESTS_H

#include <iostream> // std::cout
#include <string>
#include <chrono>
#include <vector>
#include <functional>

class Test;

// A test-runner, defined statically
class TestList {
	std::vector<Test> tests;
	std::vector<Test*> currentlyRunning;
	bool currentlySilent = false;
public:
	void add(Test& test);
	int run(int repeats=1);
	void fail(std::string reason);
};

// A test object with an associated function, which adds itself to the above list
class Test {
	using TestFn = std::function<void(Test&)>;

	TestList& testList;
	std::string codeLocation;
	std::string testName;
	TestFn runFn;

	bool running = false;
public:
	Test(TestList& testList, std::string codeLocation, std::string testName, TestFn fn) : testList(testList), codeLocation(codeLocation), testName(testName), runFn(fn) {
		testList.add(*this);
	}
	void run(int depth, bool silent=false);

	bool success = true;
	std::string reason;
	void fail(std::string r="") {
		if (!success) return;
		success = false;
		reason = r;
		testList.fail(reason);
	}
	void pass() {}
};

#define TEST_VAR_NAME test
// A macro to define a new test
// Use like:	TEST(some_unique_name) {...}
#define TEST(description, uniqueName) \
	static void test_##uniqueName (Test &test); \
	static Test Test_##uniqueName {_globalTestList, std::string(__FILE__ ":") + std::to_string(__LINE__), description, test_##uniqueName}; \
	static void test_##uniqueName (Test &TEST_VAR_NAME)
// Use if defining test inside a struct (e.g. for templating)
#define TEST_METHOD(description, uniqueName) \
	static void test_##uniqueName (Test &test) { \
		testbody_##uniqueName(test); \
	} \
	Test Test_##uniqueName {_globalTestList, std::string(__FILE__ ":") + std::to_string(__LINE__), description, test_##uniqueName}; \
	static void testbody_##uniqueName (Test &TEST_VAR_NAME)

#define FAIL(reason) TEST_VAR_NAME.fail(reason)

extern TestList _globalTestList;

/***** Benchmarking stuff *****/

class Timer {
	std::chrono::high_resolution_clock::time_point startTime;
	double totalTime = 0;
	int segmentCount = 0;
	double scaleFactor = 1;
public:
	void start() {
		startTime = std::chrono::high_resolution_clock::now();
	}
	double stop() {
		std::chrono::duration<double> duration = std::chrono::high_resolution_clock::now() - startTime;
		segmentCount++;
		return (totalTime += duration.count());
	}
	void clear() {
		totalTime = 0;
		segmentCount = 0;
		scaleFactor = 1;
	}
	void scale(double scale) {
		scaleFactor *= scale;
	}
	void scaleRate(double scale) {
		scaleFactor /= scale;
	}
	double time() const {
		return totalTime;
	}
	double scaledTime() const {
		return totalTime*scaleFactor;
	}
	double segments() const {
		return segmentCount;
	}
};

/*
Executes a test function with an increasing number of repeats
until a certain amount of time is spent on the computation.

Performs a few repeated measurements with shorter periods,
and collects the fastest.

Example use:

	BenchmarkRate trial([](int repeats, Timer &timer) {
		timer.start();
		// test code
		timer.stop();
	});
	trial.run(1); // spend at least a second on it
	trial.fastest // also returned from .run()

Use for a range of configurations

	std::vector<int> configSize = {1, 2, 4, 8, 16};
	std::vector<double> rates = BenchmarkRate::map(configSize, [](int configSize, int repeats, Timer &timer) {
	
	});

*/

extern double defaultBenchmarkTime;
extern int defaultBenchmarkDivisions;

struct BenchmarkRate {
	using TestFunction = std::function<void(int, Timer&)>;

	TestFunction fn;
	std::vector<double> rates;
	double fastest = 0;
	double optimistic = 0;

	BenchmarkRate(TestFunction fn) : fn(fn) {}

	void clear() {
		rates.resize(0);
		fastest = 0;
	}

	double run(double targetTotalTime=0, int divisions=0) {
		if (targetTotalTime == 0) targetTotalTime = defaultBenchmarkTime;
		if (divisions == 0) divisions = defaultBenchmarkDivisions;

		Timer timer;
		double totalTime = 0;

		int repeats = 1;
		double targetBlockTime = std::min(targetTotalTime/(divisions + 1), 0.05); // 50ms blocks or less
		while (repeats < 1e10) {
			timer.clear();
			fn(repeats, timer);
			if (timer.segments() == 0) {
				std::cerr << "Benchmark function didn't call timer.start()/.stop()\n";
				// The test isn't calling the timer
				return 0;
			}

			double time = timer.time();
			totalTime += time;
			if (time >= targetBlockTime) {
				break;
			} else {
				int estimatedRepeats = repeats*targetBlockTime/(time + targetBlockTime*0.01);
				repeats = std::max(repeats*2, (int)estimatedRepeats);
			}
		}

		rates.push_back(repeats/timer.scaledTime());

		while (totalTime < targetTotalTime) {
			timer.clear();
			fn(repeats, timer);

			double time = timer.time();
			totalTime += time;
			rates.push_back(repeats/timer.scaledTime());
		}

		double sum = 0;
		for (double rate : rates) {
			fastest = std::max(fastest, rate);
			sum += rate;
		}
		double mean = sum/rates.size();

		double optimisticSum = 0;
		int optimisticCount = 0;
		for (double rate : rates) {
			if (rate >= mean) {
				optimisticSum += rate;
				optimisticCount++;
			}
		}
		optimistic = optimisticSum/optimisticCount;
		return optimistic;
	}

	template<typename Arg>
	static std::vector<double> map(std::vector<Arg> &args, std::function<void(Arg,int,Timer&)> fn, bool print=false) {
		std::vector<double> results;
		for (Arg& arg : args) {
			BenchmarkRate trial([&, arg, fn](int repeats, Timer &timer) {
				fn(arg, repeats, timer);
			});
			results.push_back(trial.run());
		}
		if (print) {
			BenchmarkRate::print(results, true);
		}
		return results;
	}

	template<typename T>
	static void print(std::vector<T> array, bool newline=true) {
		for (unsigned int i = 0; i < array.size(); ++i) {
			if (i > 0) std::cout << "\t";
			std::cout << array[i];
		}
		if (newline) std::cout << std::endl;
	}
};

#endif