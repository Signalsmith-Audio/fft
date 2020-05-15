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

BENCHMARK_CPP_FILES := benchmark/*.cpp

benchmark: out/benchmark
	mkdir -p out/results
	cd out && ./benchmark --test-time=0.01

out/benchmark: *.h $(shell find benchmark -iname "*.h") $(shell find benchmark -iname "*.cpp")
	echo "building benchmark: ${BENCHMARK_CPP_FILES}"
	mkdir -p out
	g++ -std=c++11 -Wall -Wextra -Wfatal-errors -g -O3 \
 		-Wpedantic -pedantic-errors \
		"${SHARED_PATH}/test/main.cpp" -I "${SHARED_PATH}" \
		-I benchmark/ ${BENCHMARK_CPP_FILES} \
		-o out/benchmark

############## Clean ##############

clean:
	rm -rf out
