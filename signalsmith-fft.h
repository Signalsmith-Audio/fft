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
		
		struct Step {
			size_t factor;
		};

		template<bool inverse>
		void fftStepGeneric(complex *data, const Step &step) {
			std::vector<complex> working(_size);
			for (size_t i = 0; i < _size; ++i) {
				working[i] = data[i];
			}
			for (size_t f = 0; f < _size; ++f) {
				complex sum = working[0];
				for (size_t i = 1; i < _size; ++i) {
					V phase = 2*M_PI*f*i/_size;
					complex factor = {cos(phase), -sin(phase)};
					sum += perf::complexMul<inverse>(working[i], factor);
				}
				data[f] = sum;
			}
		}

		template<bool inverse>
		void run(complex *data) {
			using std::swap;
			
			Step singleStep;
			singleStep.factor = _size;
			fftStepGeneric<inverse>(data, singleStep);
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
//				working.resize(size);

//				setPlan();
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
