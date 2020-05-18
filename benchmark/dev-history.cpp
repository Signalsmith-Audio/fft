#include <string>

#include "benchmark.h"

template<template<typename> class Historical, char const *constName, char const *constTag>
struct DevHistory {

	static std::string name() {
		return constName;
	}
	static std::string resultTag() {
		return constTag;
	}
	static std::string version() {
		return "(dev)";
	}

	template <typename T>
	struct OutOfPlace {
		Historical<T> fft;
		OutOfPlace(size_t size) : fft(size) {}

		void forward(std::complex<T> *input, std::complex<T> *output) {
			fft.fft(input, output);
		}
	};

	using DoubleOutOfPlace = OutOfPlace<double>;
};

#include "dev-history/direct.h"
char directName[] = "Direct", directTag[] = "dev-history-direct";
Benchmark<DevHistory<dev_direct::FFT, directName, directTag>, 256> benchDirect;

#include "dev-history/unspecialised-factors.h"
char unspecialisedName[] = "Unspecialised", unspecialisedTag[] = "dev-history-unspecialised";
Benchmark<DevHistory<dev_unspecialised_factors::FFT, unspecialisedName, unspecialisedTag>, (1<<16)> benchUnspecialised;

#include "dev-history/specialised-23.h"
char specialised23Name[] = "Radix-2/3", specialised23Tag[] = "dev-history-specialised23";
Benchmark<DevHistory<dev_specialised23::FFT, specialised23Name, specialised23Tag>, (1<<16)> benchSpecialised23;

#include "dev-history/specialised-235.h"
char specialised235Name[] = "Radix-2/3/5", specialised235Tag[] = "dev-history-specialised235";
Benchmark<DevHistory<dev_specialised235::FFT, specialised235Name, specialised235Tag>, (1<<20)> benchSpecialised235;

#include "dev-history/memory-ping-pong.h"
char pingpongName[] = "Ping-pong", pingpongTag[] = "dev-history-ping-pong";
Benchmark<DevHistory<dev_memory_pingpong::FFT, pingpongName, pingpongTag>, (1<<20)> benchPingPong;

#include "dev-history/repeated-twiddles.h"
char repeatedTwiddlesName[] = "Repeated twiddles", repeatedTwiddlesTag[] = "dev-history-repeated-twiddles";
Benchmark<DevHistory<dev_repeated_twiddles::FFT, repeatedTwiddlesName, repeatedTwiddlesTag>, (1<<20)> benchRepeatedTwiddles;

