#ifndef SIGNALSMITH_FFT_V5
#define SIGNALSMITH_FFT_V5
// So we can easily include multiple versions by redefining this
#ifndef SIGNALSMITH_FFT_NAMESPACE
#define SIGNALSMITH_FFT_NAMESPACE signalsmith
#endif

#include <vector>
#include <complex>
#include <cmath>
#include <array>
#include <memory>

#ifndef SIGNALSMITH_INLINE
#ifdef __GNUC__
#define SIGNALSMITH_INLINE __attribute__((always_inline)) inline
#elif defined(__MSVC__)
#define SIGNALSMITH_INLINE __forceinline inline
#else
#define SIGNALSMITH_INLINE inline
#endif
#endif

namespace SIGNALSMITH_FFT_NAMESPACE {

	namespace perf {
		// Complex multiplication has edge-cases around Inf/NaN - handling those properly makes std::complex non-inlineable, so we use our own
		template <bool conjugateSecond, typename V>
		SIGNALSMITH_INLINE std::complex<V> complexMul(const std::complex<V> &a, const std::complex<V> &b) {
			return conjugateSecond ? std::complex<V>{
				b.real()*a.real() + b.imag()*a.imag(),
				b.real()*a.imag() - b.imag()*a.real()
			} : std::complex<V>{
				a.real()*b.real() - a.imag()*b.imag(),
				a.real()*b.imag() + a.imag()*b.real()
			};
		}

		template<bool flipped, typename V>
		SIGNALSMITH_INLINE std::complex<V> complexAddI(const std::complex<V> &a, const std::complex<V> &b) {
			return flipped ? std::complex<V>{
				a.real() + b.imag(),
				a.imag() - b.real()
			} : std::complex<V>{
				a.real() - b.imag(),
				a.imag() + b.real()
			};
		}
	}

	template<typename V>
	class FFT {
		using complex = std::complex<V>;
		size_t _size;
		std::vector<complex> workingVector;
		
		enum class StepType {
			permute, generic
		};
		struct Step {
			StepType type;
			size_t factor;
			size_t startIndex;
			size_t innerRepeats;
			size_t twiddleIndex;
		};
		std::vector<size_t> factors;
		std::vector<Step> plan;
		std::vector<complex> twiddleVector;
		
		struct PermutationPair {size_t from, to;};
		std::vector<PermutationPair> permutation;
		
		void addPlanSteps(size_t factorIndex, size_t start, size_t length) {
			if (factorIndex >= factors.size()) return;
			
			size_t factor = factors[factorIndex], subLength = length/factor;
			Step mainStep{StepType::generic, factor, start, subLength, twiddleVector.size()};
			for (size_t i = 0; i < subLength; ++i) {
				for (size_t f = 0; f < factor; ++f) {
					V phase = 2*M_PI*i*f/length;
					complex twiddle = {cos(phase), -sin(phase)};
					twiddleVector.push_back(twiddle);
				}
			}
			
			plan.push_back(mainStep);
			for (size_t i = 0; i < factor; ++i) {
				addPlanSteps(factorIndex + 1, start + i*subLength, subLength);
			}
		}
		void setPlan() {
			factors.resize(0);
			size_t size = _size, f = 2;
			while (size > 1) {
				if (size%f == 0) {
					factors.push_back(f);
					size /= f;
				} else if (f > sqrt(size)) {
					f = size;
				} else {
					++f;
				}
			}

			plan.resize(0);
			twiddleVector.resize(0);
			addPlanSteps(0, 0, _size);
			
			permutation.resize(0);
			permutation.push_back(PermutationPair{0, 0});
			size_t inputStep = _size, outputStep = 1;
			for (size_t f : factors) {
				inputStep /= f;
				size_t oldSize = permutation.size();
				for (size_t i = 1; i < f; ++i) {
					for (size_t j = 0; j < oldSize; ++j) {
						PermutationPair pair = permutation[j];
						pair.from += i*inputStep;
						pair.to += i*outputStep;
						permutation.push_back(pair);
					}
				}
				outputStep *= f;
			}
			
			Step step;
			step.type = StepType::permute;
			plan.push_back(step);
		}

		template<bool inverse>
		void fftStepGeneric(complex *data, const Step &step) {
			complex *working = workingVector.data();
			for (size_t i = 0; i < step.factor*step.innerRepeats; ++i) {
				working[i] = data[i];
			}

			const complex *twiddles = twiddleVector.data() + step.twiddleIndex;
			const size_t factor = step.factor;
			const size_t stride = step.innerRepeats;
			for (size_t repeat = 0; repeat < step.innerRepeats; ++repeat) {
				for (size_t f = 0; f < factor; ++f) {
					complex sum = working[0];
					for (size_t i = 1; i < factor; ++i) {
						V phase = 2*M_PI*f*i/factor;
						complex factor = {cos(phase), -sin(phase)};
						complex value = inverse ? perf::complexMul<true>(working[i*stride], twiddles[i]) : working[i*stride];
						sum += perf::complexMul<inverse>(value, factor);
					}
					data[f*stride] = inverse ? sum : perf::complexMul<false>(sum, twiddles[f]);
				}
				++data;
				++working;
				twiddles += factor;
			}
		}

		template<bool inverse>
		void permute(complex *data) {
			complex *working = workingVector.data();
			for (size_t i = 0; i < _size; ++i) {
				working[i] = data[i];
			}

			for (auto pair : permutation) {
				if (inverse) {
					data[pair.from] = working[pair.to];
				} else {
					data[pair.to] = working[pair.from];
				}
			}
		}

		template<bool inverse>
		void run(complex *data) {
			using std::swap;
			
//			std::cout << "start: ";
//			for (size_t i = 0; i < _size; ++i) {
//				std::cout << " " << data[i];
//			}
//			std::cout << std::endl;
			
			Step *step = inverse ? plan.data() + plan.size() - 1 : plan.data();
			Step *end = inverse ? plan.data() - 1 : plan.data() + plan.size();
			for (; step != end; inverse ? --step : ++step) {
				switch (step->type) {
					case StepType::permute:
						permute<inverse>(data);
						break;
					case StepType::generic:
						fftStepGeneric<inverse>(data + step->startIndex, *step);
						break;
				}

//				std::cout << "\t";
//				for (size_t i = 0; i < _size; ++i) {
//					std::cout << " " << data[i];
//				}
//				std::cout << std::endl;
			}
		}

		static bool validSize(size_t size) {
			constexpr static bool filter[32] = {
				1, 1, 1, 1, 1, 1, 1, 0, 1, 1, // 0-9
				1, 0, 1, 0, 0, 1, 1, 0, 1, 0, // 10-19
				1, 0, 0, 0, 1, 1, 0, 0/*27*/, 0, 0, // 20-29
				1, 0
			};
			return filter[size];
		}
	public:
		static size_t fastSizeAbove(size_t size) {
			size_t power2 = 1;
			while (size >= 32) {
				size = (size - 1)/2 + 1;
				power2 *= 2;
			}
			while (size < 32 && !validSize(size)) {
				++size;
			}
			return power2*size;
		}
		static size_t fastSizeBelow(size_t size) {
			size_t power2 = 1;
			while (size >= 32) {
				size /= 2;
				power2 *= 2;
			}
			while (size > 1 && !validSize(size)) {
				--size;
			}
			return power2*size;
		}

		FFT(size_t size, int fastDirection=0) : _size(0) {
			if (fastDirection > 0) size = fastSizeAbove(size);
			if (fastDirection < 0) size = fastSizeBelow(size);
			this->setSize(size);
		}

		size_t setSize(size_t size) {
			if (size != _size) {
				_size = size;
				workingVector.resize(size);
				setPlan();
			}
			return _size;
		}
		size_t setSizeMinimum(size_t size) {
			return setSize(fastSizeAbove(size));
		}
		size_t setSizeMaximum(size_t size) {
			return setSize(fastSizeBelow(size));
		}
		const size_t & size() const {
			return size;
		}

		void fft(std::vector<complex> const &input, std::vector<complex> &output) {
			return fft(input.data(), output.data());
		}
		void fft(complex const *input, complex *output) {
			if (output != input) {
				for (size_t i = 0; i < _size; ++i) {
					output[i] = input[i];
				}
			}
			return run<false>(output);
		}

		void ifft(std::vector<complex> const &input, std::vector<complex> &output) {
			return ifft(input.data(), output.data());
		}
		void ifft(complex const *input, complex *output) {
			if (output != input) {
				for (size_t i = 0; i < _size; ++i) {
					output[i] = input[i];
				}
			}
			return run<true>(output);
		}
	};
}

#undef SIGNALSMITH_FFT_NAMESPACE
#endif // SIGNALSMITH_FFT_V5
