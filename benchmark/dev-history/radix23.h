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
			generic, step2, step3
		};
		struct Step {
			StepType type;
			size_t factor;
			size_t startIndex;
			size_t innerRepeats;
			size_t outerRepeats;
			size_t twiddleIndex;
		};
		std::vector<size_t> factors;
		std::vector<Step> plan;
		std::vector<complex> twiddleVector;
		
		struct PermutationPair {size_t from, to;};
		std::vector<PermutationPair> permutation;
		
		void addPlanSteps(size_t factorIndex, size_t start, size_t length, size_t repeats) {
			if (factorIndex >= factors.size()) return;

			size_t factor = factors[factorIndex], subLength = length/factor;
			Step mainStep{StepType::generic, factor, start, subLength, repeats, twiddleVector.size()};

			if (factor == 2) mainStep.type = StepType::step2;
			if (factor == 3) mainStep.type = StepType::step3;

			// Twiddles
			for (size_t i = 0; i < subLength; ++i) {
				for (size_t f = 0; f < factor; ++f) {
					V phase = 2*M_PI*i*f/length;
					complex twiddle = {cos(phase), -sin(phase)};
					twiddleVector.push_back(twiddle);
				}
			}

			if (repeats == 1 && sizeof(complex)*subLength > 65536) {
				for (size_t i = 0; i < factor; ++i) {
					addPlanSteps(factorIndex + 1, start + i*subLength, subLength, 1);
				}
			} else {
				addPlanSteps(factorIndex + 1, start, subLength, repeats*factor);
			}
			plan.push_back(mainStep);
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
			addPlanSteps(0, 0, _size, 1);
			
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
		}

		template<bool inverse>
		void fftStepGeneric(complex *origData, const Step &step) {
			complex *working = workingVector.data();
			const size_t stride = step.innerRepeats;

			for (size_t outerRepeat = 0; outerRepeat < step.outerRepeats; ++outerRepeat) {
				complex *data = origData;
				
				const complex *twiddles = twiddleVector.data() + step.twiddleIndex;
				const size_t factor = step.factor;
				for (size_t repeat = 0; repeat < step.innerRepeats; ++repeat) {
					for (size_t i = 0; i < step.factor; ++i) {
						working[i] = perf::complexMul<inverse>(data[i*stride], twiddles[i]);
					}
					for (size_t f = 0; f < factor; ++f) {
						complex sum = working[0];
						for (size_t i = 1; i < factor; ++i) {
							V phase = 2*M_PI*f*i/factor;
							complex factor = {cos(phase), -sin(phase)};
							sum += perf::complexMul<inverse>(working[i], factor);
						}
						data[f*stride] = sum;
					}
					++data;
					twiddles += factor;
				}
				origData += step.factor*step.innerRepeats;
			}
		}

		template<bool inverse>
		void fftStep2(complex *origData, const Step &step) {
			const size_t stride = step.innerRepeats;
			const complex *origTwiddles = twiddleVector.data() + step.twiddleIndex;
			for (size_t outerRepeat = 0; outerRepeat < step.outerRepeats; ++outerRepeat) {
				const complex* twiddles = origTwiddles;
				for (complex *data = origData; data < origData + stride; ++data) {
					complex A = data[0];
					complex B = perf::complexMul<inverse>(data[stride], twiddles[1]);
					
					data[0] = A + B;
					data[stride] = A - B;
					twiddles += 2;
				}
				origData += 2*stride;
			}
		}

		template<bool inverse>
		void fftStep3(complex *origData, const Step &step) {
			constexpr complex factor3 = {-0.5, inverse ? 0.8660254037844386 : -0.8660254037844386};
			const size_t stride = step.innerRepeats;
			const complex *origTwiddles = twiddleVector.data() + step.twiddleIndex;
			
			for (size_t outerRepeat = 0; outerRepeat < step.outerRepeats; ++outerRepeat) {
				const complex* twiddles = origTwiddles;
				for (complex *data = origData; data < origData + stride; ++data) {
					complex A = data[0];
					complex B = perf::complexMul<inverse>(data[stride], twiddles[1]);
					complex C = perf::complexMul<inverse>(data[stride*2], twiddles[2]);
					
					complex realSum = A + (B + C)*factor3.real();
					complex imagSum = (B - C)*factor3.imag();

					data[0] = A + B + C;
					data[stride] = perf::complexAddI<false>(realSum, imagSum);
					data[stride*2] = perf::complexAddI<true>(realSum, imagSum);

					twiddles += 3;
				}
				origData += 3*stride;
			}
		}

		void permute(const complex *input, complex *data) {
			for (auto pair : permutation) {
				data[pair.from] = input[pair.to];
			}
		}

		template<bool inverse>
		void run(const complex *input, complex *data) {
			permute(input, data);
			
//			std::cout << "input: ";
//			for (size_t i = 0; i < _size; ++i) {
//				std::cout << " " << input[i];
//			}
//			std::cout << std::endl;
//			std::cout << "permuted: ";
//			for (size_t i = 0; i < _size; ++i) {
//				std::cout << " " << data[i];
//			}
//			std::cout << std::endl;

			for (const Step &step : plan) {
				switch (step.type) {
					case StepType::generic:
						fftStepGeneric<inverse>(data + step.startIndex, step);
						break;
					case StepType::step2:
						fftStep2<inverse>(data + step.startIndex, step);
						break;
					case StepType::step3:
						fftStep3<inverse>(data + step.startIndex, step);
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
			return run<false>(input, output);
		}

		void ifft(std::vector<complex> const &input, std::vector<complex> &output) {
			return ifft(input.data(), output.data());
		}
		void ifft(complex const *input, complex *output) {
			return run<true>(input, output);
		}
	};
}

#undef SIGNALSMITH_FFT_NAMESPACE
#endif // SIGNALSMITH_FFT_V5
