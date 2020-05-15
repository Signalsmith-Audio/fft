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
		template <bool conjugateSecond, typename sampleType>
		SIGNALSMITH_INLINE std::complex<sampleType> complexMul(const std::complex<sampleType> &A, const std::complex<sampleType> &B) {
			if (conjugateSecond) {
				return std::complex<sampleType>{
					B.real()*A.real() + B.imag()*A.imag(),
					B.real()*A.imag() - B.imag()*A.real()
				};
			} else {
				return std::complex<sampleType>{
					A.real()*B.real() - A.imag()*B.imag(),
					A.real()*B.imag() + A.imag()*B.real()
				};
			}
		}

	}

	// template<typename T>
	// static void printArray(std::complex<T>* array, size_t size, bool newline=true) {
	// 	for (unsigned int i = 0; i < size; ++i) {
	// 		std::complex<T> value = array[i];
	// 		if (abs(value.real()) < 1e-10) value.real(0);
	// 		if (abs(value.imag()) < 1e-10) value.imag(0);
	// 		if (i > 0) std::cout << "\t";
	// 		std::cout << value;
	// 	}
	// 	if (newline) std::cout << std::endl;
	// }

	template<typename V>
	class FFT {
		using complex = std::complex<V>;
		size_t _size;
		std::vector<complex> working;

		struct Step {
			size_t N;
			size_t twiddleOffset;
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
				for (size_t divisor = 2; divisor <= sqrt(size); ++divisor) {
					if (size%divisor == 0) {
						stepSize = divisor;
						break;
					}
				}
				size_t twiddleRepeats = _size/size;
				// Calculate twiddles
				size_t twiddleOffset = twiddles.size();
				double phaseStep = 2*M_PI/size;
				for (size_t i = 0; i < size/stepSize; i++) {
					for (size_t r = 0; r < twiddleRepeats; r++) {
						for (size_t bin = 0; bin < stepSize; bin++) {
							double twiddlePhase = phaseStep*bin*i;
							twiddles.push_back({cos(twiddlePhase), -sin(twiddlePhase)});
						}
					}
				}	

				plan.push_back({stepSize, twiddleOffset});
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
			for (size_t offset = 0; offset < stride; offset++) {
				for (size_t bin = 0; bin < stepSize; ++bin) {
					complex sum = input[0];
					for (size_t i = 1; i < stepSize; ++i) {
						V phase = 2*M_PI*bin*i/stepSize;
						complex factor = {cos(phase), -sin(phase)};
						sum += perf::complexMul<inverse>(input[i*stride], factor);
					}

					output[bin] = perf::complexMul<inverse>(sum, twiddles[bin]);
				}
				twiddles += stepSize;

				input += 1;
				output += stepSize;
			}
		}

		template<bool inverse, size_t stepSize>
		static std::array<complex, stepSize*stepSize> factorsFixed() {
			std::array<complex, stepSize*stepSize> result;
			for (size_t i = 0; i < stepSize; i++) {
				for (size_t j = 0; j < stepSize; j++) {
					double phase = 2*M_PI*i*j/stepSize;
					result[i*stepSize + j] = {cos(phase), inverse ? sin(phase) : -sin(phase)};
				}
			}
			return result;
		}

		template<bool inverse, int stepSize>
		void fftStepFixed(complex const *input, complex *output, const Step &step) {
			static const std::array<complex, stepSize*stepSize> factors = factorsFixed<inverse, stepSize>();

			const complex *twiddles = &this->twiddles[step.twiddleOffset];

			size_t stride = _size/stepSize;
			for (size_t offset = 0; offset < stride; offset++) {
				complex sum = input[0];
				for (size_t i = 1; i < stepSize; ++i) {
					sum += input[i*stride];
				}
				output[0] = sum;

				for (size_t bin = 1; bin < stepSize; ++bin) {
					complex sum = input[0];
					for (size_t i = 1; i < stepSize; ++i) {
						sum += perf::complexMul<false>(input[i*stride], factors[bin*stepSize + i]);
					}
					output[bin] = perf::complexMul<inverse>(sum, twiddles[bin]);
				}

				input += 1;
				output += stepSize;
				twiddles += stepSize;
			}
		}

		template<bool inverse>
		void fftStep2(complex const *input0, complex *output, const Step &step) {
			const complex *twiddles = &this->twiddles[step.twiddleOffset];
			size_t stride = _size/2;

			complex const *input1 = input0 + stride;
			complex const *end = input1;
			while (input0 != end) {
				complex sum = *input0 + *input1;
				complex diff = *input0 - *input1;

				const complex &twiddle = twiddles[1];
				twiddles += 2;

				output[0] = sum;
				output[1] = perf::complexMul<inverse>(diff, twiddle);
				++input0;
				++input1;
				output += 2;
			}
		}

		template<bool inverse>
		void fftStep3(complex const *input0, complex *output, const Step &step) {
			const complex factor3 = {-0.5, inverse ? 0.8660254037844386 : -0.8660254037844386};

			const complex *twiddles = &this->twiddles[step.twiddleOffset];
			size_t stride = _size/3;

			complex const *input1 = input0 + stride;
			complex const *input2 = input0 + stride*2;
			complex const *end = input1;
			while (input0 != end) {
				complex A = *input0, B = *input1, C = *input2;
				complex realSum = A + (B + C)*factor3.real();
				complex imagSum = (B - C)*factor3.imag();

				const complex &twiddle1 = twiddles[1];
				const complex &twiddle2 = twiddles[2];
				twiddles += 3;

				output[0] = A + B + C;
				output[1] = perf::complexMul<inverse>(complex{realSum.real() - imagSum.imag(), realSum.imag() + imagSum.real()}, twiddle1);
				output[2] = perf::complexMul<inverse>(complex{realSum.real() + imagSum.imag(), realSum.imag() - imagSum.real()}, twiddle2);
				++input0;
				++input1;
				++input2;
				output += 3;
			}
		}
		template<bool inverse>
		void run(complex const *input, complex *output) {
			using std::swap;

			// Copy the input in
			complex *A = output, *B = working.data();
			for (size_t i = 0; i < _size; i++) {
				A[i] = input[i];
			}
	
			// Go through the steps
			for (const Step& step : plan) {
				if (step.N == 2) {
					fftStepFixed<inverse, 2>(A, B, step);
					// fftStep2<inverse>(A, B, step);
				} else if (step.N == 3) {
					// fftStepFixed<inverse, 3>(A, B, step);
					fftStep3<inverse>(A, B, step);
				} else {
					fftStepGeneric<inverse>(A, B, step);
				}
				swap(A, B);
			}

			// Un-permute the result
			for (size_t i = 0; i < _size; i++) {
				B[permutation[i]] = A[i];
			}
			swap(A, B);

			// Copy the output
			if (A != output) {
				for (size_t i = 0; i < _size; i++) {
					output[i] = A[i];
				}
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
