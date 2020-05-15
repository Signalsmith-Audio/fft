#include <string>

#include "benchmark.h"

struct Signalsmith {
	static std::string name() {
		return "Signalsmith";
	}
	static std::string resultTag() {
		return "signalsmith";
	}
	static std::string version() {
		return "(development)";
	}

	template <typename T>
	struct OutOfPlace {
		signalsmith::FFT<T> fft;
		OutOfPlace(size_t size) : fft(size) {}

		void forward(std::complex<T> *input, std::complex<T> *output) {
			fft.fft(input, output);
		}
	};

	using DoubleOutOfPlace = OutOfPlace<double>;
};
Benchmark<Signalsmith> benchSignalsmith;

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
Benchmark<DevHistory<dev_unspecialised_factors::FFT, unspecialisedName, unspecialisedTag>, (1<<18)> benchUnspecialised;

