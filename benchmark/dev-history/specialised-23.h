#include <vector>
#include <complex>
#include <cmath>

namespace dev_specialised23 {
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
						sum += (inverse ? conj(factor) : factor)*input[i*stride];
					}

					const complex &twiddle = *twiddles;
					output[bin] = sum*(inverse ? conj(twiddle) : twiddle);
					twiddles++;
				}

				input += 1;
				output += stepSize;
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
				output[1] = diff*(inverse ? conj(twiddle) : twiddle);
				++input0;
				++input1;
				output += 2;
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
					fftStep2<inverse>(A, B, step);
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