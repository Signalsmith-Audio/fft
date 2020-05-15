#include <vector>
#include <complex>
#include <cmath>

#ifndef LOG
#define LOG 0
#endif

namespace signalsmith {
	template<typename T>
	static void printArray(std::complex<T>* array, size_t size, bool newline=true) {
		for (unsigned int i = 0; i < size; ++i) {
			std::complex<T> value = array[i];
			if (abs(value.real()) < 1e-10) value.real(0);
			if (abs(value.imag()) < 1e-10) value.imag(0);
			if (i > 0) std::cout << "\t";
			std::cout << value;
		}
		if (newline) std::cout << std::endl;
	}

	template<typename V>
	class FFT {
		using complex = std::complex<V>;
		size_t _size;
		std::vector<complex> working;

		template<bool inverse>
		void fftStepGeneric(complex const *input, complex *output, size_t size) {
			size_t stride = _size/size;
			for (size_t offset = 0; offset < stride; offset++) {
				for (size_t bin = 0; bin < size; ++bin) {
					complex sum = 0;
					for (size_t i = 0; i < size; ++i) {
						V phase = 2*M_PI*bin*i/size;
						complex factor = {cos(phase), inverse ? sin(phase) : -sin(phase)};
						sum += factor*input[i*stride];
					}
					output[bin] = sum;
				}
				input += 1;
				output += size;
			}
		}

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
				if (size%2 == 0) {
					stepSize = 2;
				}
				plan.push_back({stepSize, twiddles.size()});
				double phaseStep = 2*M_PI/size;
				for (size_t i = 0; i < size/stepSize; i++) {
					for (size_t r = 0; r < _size/size; r++) {
						for (size_t bin = 0; bin < stepSize; bin++) {
							double twiddlePhase = phaseStep*bin*i;
							twiddles.push_back({cos(twiddlePhase), -sin(twiddlePhase)});
						}
					}
				}	
				size /= stepSize;
			}

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
		void run(complex const *input, complex *output) {
			using std::swap;

			// Copy the input in
			complex *A = output, *B = working.data();
			for (size_t i = 0; i < _size; i++) {
				A[i] = input[i];
			}
	
			// Go through the steps
			for (const Step& step : plan) {
				fftStepGeneric<inverse>(A, B, step.N);
				for (size_t i = 0; i < _size; i++) {
					complex &twiddle = twiddles[step.twiddleOffset + i];
					B[i] *= (inverse ? conj(twiddle) : twiddle);
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