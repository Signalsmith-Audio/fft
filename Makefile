.PHONY: test
ifndef VERBOSE
.SILENT:
endif

SHARED_PATH := "geraint"

############## Testing ##############

TEST_CPP_FILES := $(shell find tests -iname "*.cpp" | sort)

test: out/test
	./out/test

out/test: *.h $(shell find tests -iname "*.h") $(shell find tests -iname "*.cpp")
	echo "building tests: ${TEST_CPP_FILES}"
	mkdir -p out
	g++ -std=c++11 -Wall -Wextra -Wfatal-errors -g -O0 \
 		-Wpedantic -pedantic-errors \
		"${SHARED_PATH}/test/main.cpp" -I "${SHARED_PATH}" \
		-I tests/ ${TEST_CPP_FILES} \
		-o out/test

############## Benchmark ##############

BENCHMARK_CPP_FILES := benchmark/benchmark.cpp

benchmark: out/benchmark
	mkdir -p out/results
	cd out && ./benchmark --test-time=0.01

out/benchmark: *.h $(shell find benchmark -iname "*.h") $(shell find benchmark -iname "*.cpp")
	echo "building benchmark: ${BENCHMARK_CPP_FILES}"
	mkdir -p out
	g++ -std=c++11 -msse2 -mavx -g -O3 \
		"${SHARED_PATH}/test/main.cpp" -I "${SHARED_PATH}" \
		-I benchmark/ ${BENCHMARK_CPP_FILES} \
		-o out/benchmark

############## Benchmark against others ##############

BENCHMARK_ALL_CPP_FILES := benchmark/*.cpp

benchmark-all: out/benchmark-all
	mkdir -p out/results
	cd out && ./benchmark-all --test-time=0.02

out/benchmark-all: *.h $(shell find benchmark -iname "*.h") $(shell find benchmark -iname "*.cpp")
	echo "building benchmark: ${BENCHMARK_ALL_CPP_FILES}"
	mkdir -p out
	g++ -std=c++11 -msse2 -mavx -Wfatal-errors -g -O3 \
		"${SHARED_PATH}/test/main.cpp" -I "${SHARED_PATH}" \
		-I benchmark/ ${BENCHMARK_ALL_CPP_FILES} \
		-o out/benchmark-all

############## Clean ##############

clean:
	rm -rf out
