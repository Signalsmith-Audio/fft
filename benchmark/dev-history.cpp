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
	struct OutOfPlace : public OutOfPlaceRunner<T> {
		Historical<T> fft;
		OutOfPlace(size_t size) : OutOfPlaceRunner<T>(size), fft(size) {}

		void forward(std::complex<T> *input, std::complex<T> *output) {
			fft.fft(input, output);
		}
	};

	using DoubleOutOfPlace = OutOfPlace<double>;
};

#undef SIGNALSMITH_FFT_V5
#undef SIGNALSMITH_FFT_NAMESPACE
#define SIGNALSMITH_FFT_NAMESPACE signalsmith_direct
#include "dev-history/direct.h"
char directName[] = "Direct", directTag[] = "dev-history-direct";
Benchmark<DevHistory<signalsmith_direct::FFT, directName, directTag>, 256> benchDirect;
