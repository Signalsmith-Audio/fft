.PHONY: test
ifndef VERBOSE
.SILENT:
endif

SHARED_PATH := "geraint"

clean:
	rm -rf out

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

############## Benchmarking ##############

graphs:
	python benchmark/graphs.py

benchmark: test benchmark-main benchmark-kissfft benchmark-previous benchmark-dev-history

BENCHMARK_TEST_TIME := 0.05

# Generic versions

benchmark-%: out/benchmark-%
	mkdir -p out/results
	cd out && ./benchmark-$* --test-time=${BENCHMARK_TEST_TIME}

out/benchmark-%: *.h $(shell find benchmark -iname "*.h") $(shell find benchmark -iname "*.cpp")
	mkdir -p out
	g++ -std=c++11 -msse2 -mavx -Wfatal-errors -g -O3 \
		"${SHARED_PATH}/test/main.cpp" -I "${SHARED_PATH}" \
		-I benchmark/ benchmark/$*.cpp \
		-o out/benchmark-$*

out/benchmark-fftw: $(shell find benchmark -iname "*.h") $(shell find benchmark -iname "*.cpp")
	mkdir -p out
	g++ -std=c++11 -msse2 -mavx -Wfatal-errors -g -O3 \
		"${SHARED_PATH}/test/main.cpp" -I "${SHARED_PATH}" \
		-I benchmark/ benchmark/fftw.cpp \
		-lfftw3 \
		-o out/benchmark-fftw
