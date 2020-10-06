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

namespace signalsmith_v4 {

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

	template<typename V, size_t size>
	class FixedFFT {
		
	};

	template<typename V>
	class FFT {
		using complex = std::complex<V>;
		size_t _size;
		std::vector<complex> working;

		size_t finalFixedSize;
		size_t finalFixedPingPong = 0;

		enum class StepType {
			butterfly2, butterfly3, butterfly4, butterfly5, butterflyGeneric
		};
		
		struct Step {
			StepType type;
			size_t N;
			size_t twiddleOffset;
			size_t twiddleRepeats;
		};
		std::vector<Step> plan;
		std::vector<complex> twiddles;
		struct PermutionPair {size_t from, to;};
		std::vector<PermutionPair> permutationSequence;

		void setPlan() {
			plan.resize(0);
			twiddles.resize(0);
			size_t size = _size/finalFixedSize;

			while (size > 1) {
				size_t stepSize = size;
				if (size%4 == 0) {
					stepSize = 4;
				} else {
					for (size_t divisor = 2; divisor <= sqrt(size); ++divisor) {
						if (size%divisor == 0) {
							stepSize = divisor;
							break;
						}
					}
				}
				size_t twiddleRepeats = _size/size;
				// Calculate twiddles
				size_t twiddleOffset = twiddles.size();
				double phaseStep = 2*M_PI/size;
				for (size_t i = 0; i < size/stepSize; i++) {
					for (size_t bin = 0; bin < stepSize; bin++) {
						double twiddlePhase = phaseStep*bin*i;
						twiddles.push_back({cos(twiddlePhase), -sin(twiddlePhase)});
					}
				}	

				StepType type = StepType::butterflyGeneric;
				if (stepSize == 2) type = StepType::butterfly2;
				if (stepSize == 3) type = StepType::butterfly3;
				if (stepSize == 4) type = StepType::butterfly4;
				if (stepSize == 5) type = StepType::butterfly5;
				plan.push_back({type, stepSize, twiddleOffset, twiddleRepeats});
				size /= stepSize;
			}

			// Construct optimal cache-invariant permutation order
			permutationSequence.resize(1);
			permutationSequence[0] = {0, 0};
			size_t lowProduct = 1, highDivisor = 1;
			size_t lowIndex = 0, highIndex = plan.size();
			while (lowProduct*highDivisor != _size) {
				size_t factor = 0, fromStep = 0;
				if (lowProduct <= highDivisor) {
					factor = plan[lowIndex].N;
					++lowIndex;
					lowProduct *= factor;
					fromStep = lowProduct/factor;
				} else {
					--highIndex;
					factor = plan[highIndex].N;
					highDivisor *= factor;
					fromStep = _size/highDivisor;
				}
				const size_t toStep = _size/fromStep/factor;
				const size_t oldSize = permutationSequence.size();
				for (size_t i = 1; i < factor; ++i) {
					for (size_t j = 0; j < oldSize; ++j) {
						auto pair = permutationSequence[j];
						permutationSequence.push_back({
							pair.from + i*fromStep,
							pair.to + i*toStep
						});
					}
				}
			}
		}

		template<bool inverse>
		void fftStepGeneric(complex const *input, complex *output, const Step &step) {
			size_t stepSize = step.N;
			const complex *twiddles = &this->twiddles[step.twiddleOffset];

			size_t stride = _size/stepSize;
			const complex *end = input + stride;
			while (input != end) {
				for (size_t repeat = 0; repeat < step.twiddleRepeats; ++repeat) {
					for (size_t bin = 0; bin < stepSize; ++bin) {
						complex sum = input[0];
						for (size_t i = 1; i < stepSize; ++i) {
							V phase = 2*M_PI*bin*i/stepSize;
							complex factor = {cos(phase), -sin(phase)};
							sum += perf::complexMul<inverse>(input[i*stride], factor);
						}

						output[bin] = perf::complexMul<inverse>(sum, twiddles[bin]);
					}
					++input;
					output += stepSize;
				}
				twiddles += stepSize;
			}
		}

		template<bool inverse>
		void fftStep2(complex const *input, complex *output, const Step &step) {
			const complex *twiddles = &this->twiddles[step.twiddleOffset];
			size_t stride = _size/2;

			complex const *end = input + stride;
			while (input != end) {
				for (size_t repeat = 0; repeat < step.twiddleRepeats; ++repeat) {
					complex A = input[0], B = input[stride];

					output[0] = A + B;
					output[1] = perf::complexMul<inverse>(A - B, twiddles[1]);
					++input;
					output += 2;
				}
				twiddles += 2;
			}
		}

		template<bool inverse>
		void fftStep3(complex const *input, complex *output, const Step &step) {
			const complex factor3 = {-0.5, inverse ? 0.8660254037844386 : -0.8660254037844386};

			const complex *twiddles = &this->twiddles[step.twiddleOffset];
			size_t stride = _size/3;

			complex const *end = input + stride;
			while (input != end) {
				for (size_t repeat = 0; repeat < step.twiddleRepeats; ++repeat) {
					complex A = input[0], B = input[stride], C = input[stride*2];
					complex realSum = A + (B + C)*factor3.real();
					complex imagSum = (B - C)*factor3.imag();

					output[0] = A + B + C;
					output[1] = perf::complexMul<inverse>(perf::complexAddI<false>(realSum, imagSum), twiddles[1]);
					output[2] = perf::complexMul<inverse>(perf::complexAddI<true>(realSum, imagSum), twiddles[2]);
					++input;
					output += 3;
				}
				twiddles += 3;
			}
		}

		template<bool inverse>
		void fftStep4(complex const *input, complex *output, const Step &step) {
			const complex *twiddles = &this->twiddles[step.twiddleOffset];
			size_t stride = _size/4;

			complex const *end = input + stride;
			while (input != end) {
				for (size_t repeat = 0; repeat < step.twiddleRepeats; ++repeat) {
					complex A = input[0], B = input[stride], C = input[stride*2], D = input[stride*3];

					complex sumAC = A + C, sumBD = B + D;
					complex diffAC = A - C, diffBD = B - D;

					output[0] = sumAC + sumBD;
					output[1] = perf::complexMul<inverse>(perf::complexAddI<!inverse>(diffAC, diffBD), twiddles[1]);
					output[2] = perf::complexMul<inverse>(sumAC - sumBD, twiddles[2]);
					output[3] = perf::complexMul<inverse>(perf::complexAddI<inverse>(diffAC, diffBD), twiddles[3]);
					++input;
					output += 4;
				}
				twiddles += 4;
			}
		}

		template<bool inverse>
		void fftStep5(complex const *input, complex *output, const Step &step) {
			const complex factor5a = {0.30901699437494745, inverse ? 0.9510565162951535 : -0.9510565162951535};
			const complex factor5b = {-0.8090169943749473, inverse ? 0.5877852522924732 : -0.5877852522924732};

			const complex *twiddles = &this->twiddles[step.twiddleOffset];
			size_t stride = _size/5;

			complex const *end = input + stride;
			while (input != end) {
				for (size_t repeat = 0; repeat < step.twiddleRepeats; ++repeat) {
					complex A = input[0], B = input[stride], C = input[stride*2], D = input[stride*3], E = input[stride*4];
					complex realSum1 = A + (B + E)*factor5a.real() + (C + D)*factor5b.real();
					complex imagSum1 = (B - E)*factor5a.imag() + (C - D)*factor5b.imag();
					complex realSum2 = A + (B + E)*factor5b.real() + (C + D)*factor5a.real();
					complex imagSum2 = (B - E)*factor5b.imag() + (D - C)*factor5a.imag();

					output[0] = A + B + C + D + E;
					output[1] = perf::complexMul<inverse>(perf::complexAddI<false>(realSum1, imagSum1), twiddles[1]);
					output[2] = perf::complexMul<inverse>(perf::complexAddI<false>(realSum2, imagSum2), twiddles[2]);
					output[3] = perf::complexMul<inverse>(perf::complexAddI<true>(realSum2, imagSum2), twiddles[3]);
					output[4] = perf::complexMul<inverse>(perf::complexAddI<true>(realSum1, imagSum1), twiddles[4]);
					++input;
					output += 5;
				}
				twiddles += 5;
			}
		}

		template<bool inverse>
		void run(complex const *input, complex *output) {
			using std::swap;

			const complex *A = input;
			// Choose the starting state for the ping-pong pattern, so the result ends up in the output
			bool oddSteps = (plan.size()%2);
			complex *B = oddSteps ? working.data() : output;
			complex *other = oddSteps ? output : working.data();
	
			// Go through the steps
			for (const Step& step : plan) {
				switch (step.type) {
				case StepType::butterfly2:
					fftStep2<inverse>(A, B, step);
					break;
				case StepType::butterfly3:
					fftStep3<inverse>(A, B, step);
					break;
				case StepType::butterfly4:
					fftStep4<inverse>(A, B, step);
					break;
				case StepType::butterfly5:
					fftStep5<inverse>(A, B, step);
					break;
				case StepType::butterflyGeneric:
					fftStepGeneric<inverse>(A, B, step);
					break;
				}

				A = B;
				swap(B, other);
			}

			for (auto &pair : permutationSequence) {
				B[pair.from] = A[pair.to];
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
				finalFixedSize = 1;
				working.resize(size);

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
