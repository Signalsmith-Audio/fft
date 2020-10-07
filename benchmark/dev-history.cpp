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

#undef SIGNALSMITH_FFT_V5
#undef SIGNALSMITH_FFT_NAMESPACE
#define SIGNALSMITH_FFT_NAMESPACE signalsmith_factorise
#include "dev-history/factorise.h"
char factoriseName[] = "Factorise", factoriseTag[] = "dev-history-factorise";
Benchmark<DevHistory<signalsmith_factorise::FFT, factoriseName, factoriseTag>, 8192> benchFactorise;

#undef SIGNALSMITH_FFT_V5
#undef SIGNALSMITH_FFT_NAMESPACE
#define SIGNALSMITH_FFT_NAMESPACE signalsmith_radix23
#include "dev-history/radix23.h"
char radix23Name[] = "Radix 2/3", radix23Tag[] = "dev-history-radix23";
Benchmark<DevHistory<signalsmith_radix23::FFT, radix23Name, radix23Tag>, 1048576> benchRadix23;

#undef SIGNALSMITH_FFT_V5
#undef SIGNALSMITH_FFT_NAMESPACE
#define SIGNALSMITH_FFT_NAMESPACE signalsmith_radix234
#include "dev-history/radix234.h"
char radix234Name[] = "Radix 2/3/4", radix234Tag[] = "dev-history-radix234";
Benchmark<DevHistory<signalsmith_radix234::FFT, radix234Name, radix234Tag>> benchRadix234;
