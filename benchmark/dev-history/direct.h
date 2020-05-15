#include <vector>
#include <complex>
#include <cmath>

namespace dev_direct {
	template<typename V>
	class FFT {
		using complex = std::complex<V>;
		size_t _size;
	public:
		FFT(size_t size) : _size(size) {}

		const size_t & size() const {
			return size;
		}

		void fft(std::vector<complex> const &input, std::vector<complex> &output) {
			return fft(input.data(), output.data());
		}
		void fft(complex const *input, complex *output) {
			for (size_t bin = 0; bin < _size; ++bin) {
				complex sum = 0;
				for (size_t i = 0; i < _size; ++i) {
					V phase = 2*M_PI*i*bin/_size;
					complex factor = {cos(phase), -sin(phase)};
					sum += factor*input[i];
				}
				output[bin] = sum;
			}
		}
	};
}