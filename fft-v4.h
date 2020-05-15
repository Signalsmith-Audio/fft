#include <vector>
#include <complex>
#include <cmath>

namespace signalsmith {
	template<typename V>
	class FFT {
		using complex = std::complex<V>;
		size_t _size;
		std::vector<complex> working;

		template<bool inverse>
		void innerFft(complex const *input, complex *output) {
			for (size_t bin = 0; bin < _size; ++bin) {
				complex sum = 0;
				for (size_t i = 0; i < _size; ++i) {
					V phase = 2*M_PI*bin*i/_size;
					complex factor = {cos(phase), inverse ? sin(phase) : -sin(phase)};
					sum += factor*input[i];
				}
				output[bin] = sum;
			}
		}

	public:
		FFT(size_t size) : _size(size) {}

		const size_t & size() const {
			return size;
		}

		void fft(std::vector<complex> const &input, std::vector<complex> &output) {
			return fft(input.data(), output.data());
		}
		void fft(complex const *input, complex *output) {
			return innerFft<false>(input, output);
		}

		void ifft(std::vector<complex> const &input, std::vector<complex> &output) {
			return ifft(input.data(), output.data());
		}
		void ifft(complex const *input, complex *output) {
			return innerFft<true>(input, output);
		}
	};
}