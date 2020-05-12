.PHONY: main
ifndef VERBOSE
.SILENT:
endif

############## Main build ##############

main: out/hello
	./out/hello

out/hello: *.h hello-world.cpp
	mkdir -p out
	g++ -std=c++11 -Wall -Wextra -Wfatal-errors -g -O3 \
 		-Wpedantic -pedantic-errors \
		hello-world.cpp \
		-o out/hello

############## Testing ##############

TEST_CPP_FILES := $(shell find tests -iname "*.cpp" | sort)
SHARED_PATH := "geraint"

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

############## Clean ##############

clean:
	rm -rf out
