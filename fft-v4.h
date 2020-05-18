#include <vector>
#include <complex>
#include <cmath>
#include <array>

#ifndef SIGNALSMITH_INLINE
#define SIGNALSMITH_INLINE __attribute__((always_inline)) inline
#endif

namespace signalsmith {

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
		std::vector<complex> working;

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
		std::vector<int> permutation;
		void setPlan() {
			plan.resize(0);
			twiddles.resize(0);
			size_t size = _size;
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

			// Construct permutation from factorised sizes
			permutation.resize(1);
			permutation[0] = 0;
			for (int i = plan.size() - 1; i >= 0; --i) {
				const Step &step = plan[i];
				size_t pSize = permutation.size();
				size_t stride = _size/pSize/step.N;
				for (size_t i = 1; i < step.N; i++) {
					for (size_t j = 0; j < pSize; j++) {
						permutation.push_back(permutation[j] + i*stride);
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

			// Permutation
			for (size_t i = 0; i < _size; i++) {
				B[permutation[i]] = A[i];
			}
		}

	public:
		FFT(size_t size) : _size(0) {
			this->setSize(size);
		}

		size_t setSize(size_t size) {
			if (size != _size) {
				_size = size;
				working.resize(size);

				setPlan();
			}
			return _size;
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
