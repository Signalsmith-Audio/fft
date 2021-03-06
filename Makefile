.PHONY: test
ifndef VERBOSE
.SILENT:
endif

SHARED_PATH := "common"

clean:
	rm -rf out

############## Testing ##############

TEST_CPP_FILES := $(shell find tests -iname "*.cpp" | sort)

test: out/test
	./out/test

out/test: *.h $(shell find tests -iname "*.h") $(shell find tests -iname "*.cpp")
	echo "building tests: ${TEST_CPP_FILES}"
	mkdir -p out
	g++ -std=c++11 -Wall -Wextra -Wfatal-errors -O0 \
 		-Wpedantic -pedantic-errors \
		"${SHARED_PATH}/test/main.cpp" -I "${SHARED_PATH}" \
		-I tests/ ${TEST_CPP_FILES} \
		-o out/test

############## Benchmarking ##############

graphs:
	python benchmark/graphs.py

benchmarks: test benchmark-main benchmark-kissfft benchmark-fftw graphs

BENCHMARK_TEST_TIME := 0.05

# Generic versions

benchmark-%: out/benchmark-%
	mkdir -p out/results
	cd out && ./benchmark-$* --test-time=${BENCHMARK_TEST_TIME}

out/benchmark-%: *.h $(shell find benchmark -iname "*.h") $(shell find benchmark -iname "*.cpp")
	mkdir -p out
	g++ -std=c++11 -msse2 -mavx -Wfatal-errors -O3 \
		"${SHARED_PATH}/test/main.cpp" -I "${SHARED_PATH}" \
		-I benchmark/ benchmark/$*.cpp \
		-o out/benchmark-$*

# Custom versions which need more config

out/benchmark-fftw: $(shell find benchmark -iname "*.h") $(shell find benchmark -iname "*.cpp")
	mkdir -p out
	g++ -std=c++11 -msse2 -mavx -Wfatal-errors -g -O3 \
		"${SHARED_PATH}/test/main.cpp" -I "${SHARED_PATH}" \
		-I benchmark/ benchmark/fftw.cpp \
		-lfftw3 \
		-o out/benchmark-fftw

############## Development ##############

# After an initial "make benchmark", you can continue
dev: test benchmark-main graphs
