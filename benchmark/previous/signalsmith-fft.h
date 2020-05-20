#pragma once
#ifndef SIGNALSMITH_FFT_H_
#define SIGNALSMITH_FFT_H_

#include <vector> // for std::vector
#include <memory> // for std::shared_ptr
#include <cstring> // for std::memcpy
#include <math.h> // for sin()/cos()

/*
Signalsmith Audio's FFT library

Recommended use:

	auto fft = getFft<double>(32768); // returns a std::shared_ptr<FFT<double>>

	fft->fft(buffer); // In-place (recommended)
	fft->ifft(buffer);

	fft->fft(input, spectrum);
	fft->ifft(spectrum, input);

The .fft() method produces permuted output, and .ifft() expects permuted input.
There are methods to re-order the computed spectrum:

	fft->permute(permuted, unpermuted); // Out-of-place (recommended)
	fft->ipermute(unpermuted, permuted);

	fft->permute(result); // In-place
	fft->ipermute(result);

Not all sizes are supported. The actual size may be up to 33% larger than requested.
The actual size is returned from .size() or .setSize():

	int currentSize = fft.size();
	int actualSize = fft.setSize(16000); // returns 16384

You can also provide a second argument to setSize()/getFft() for parallel FFTs:

	auto fft = getFft<double>(256, 16); // 16 interleaved 256-point FFTs
	fft.setSize(512, 8); // 8 interleaved 512-point FFTs

If .setSize() doesn't actually change the size, it doesn't reallocate/recalculate.
If it does change, it costs as much as constructing a new object with getFft().

In a multi-threaded environment, each thread should have its own instance,
because of temporary memory used for out-of-place/in-place shuffling.

*/

// required for usable performance in debug releases
#define INLINE __attribute__((always_inline)) inline

namespace signalsmith {
namespace fft {

// complex pair
template <typename sampleType>
struct complex {
	sampleType real;
	sampleType imag;
};

namespace _internal {
	// Pre-declare COBRA method so we can use it in the FFT interface methods
	template <typename complexValue, bool inverse, int fixedParallel=0>
	void cobraPermute(complexValue* input, complexValue* output, int N, int parallel, const std::vector<int>& factors, const int* indexes);
}

// Generic FFT interface - transforming complex pairs of samples
template <typename sampleType>
class FFT {
	std::vector<complex<sampleType>> scratch;
	complex<sampleType>* scratchCopy(complex<sampleType>* input) {
		int scratchSize = size()*parallel();
		scratch.resize(scratchSize);
		memcpy(scratch.data(), input, sizeof(complexValue)*scratchSize);
		return scratch.data();
	}
public:
	using complexValue = complex<sampleType>;

	virtual ~FFT() {};

	virtual int size() = 0;
	virtual int parallel() = 0;
	virtual int setSize(int N, int parallel=1) = 0;
	virtual const std::vector<int>& getFactors() = 0;
	virtual const int* getPermutation() = 0;
	virtual void fft(complexValue* inputPairs) = 0;
	virtual void ifft(complexValue* inputPairs) = 0;

	virtual void permute(complexValue* input, complexValue* output) {
		if (input == output) return permute(input);
		_internal::cobraPermute<complexValue, false>(input, output, size(), parallel(), getFactors(), getPermutation());
	}
	virtual void ipermute(complexValue* input, complexValue* output) {
		if (input == output) return ipermute(input);
		_internal::cobraPermute<complexValue, true>(input, output, size(), parallel(), getFactors(), getPermutation());
	}

	// Convenience methods - type-friendliness and out-of-place/in-place adapters
	void fft(sampleType* inputPairs) {
		fft((complexValue*)inputPairs);
	};
	void fft(sampleType* inputPairs, sampleType* outputPairs) {
		if (inputPairs != outputPairs) std::memcpy(outputPairs, inputPairs, size()*parallel()*sizeof(complexValue));
		fft(outputPairs);
	}
	void ifft(sampleType* inputPairs) {
		ifft((complexValue*)inputPairs);
	};
	void ifft(sampleType* inputPairs, sampleType* outputPairs) {
		if (inputPairs != outputPairs) std::memcpy(outputPairs, inputPairs, size()*parallel()*sizeof(complexValue));
		ifft(outputPairs);
	}

	void permute(sampleType* inputPairs, sampleType* outputPairs) {
		permute((complexValue*)inputPairs, (complexValue*)outputPairs);
	}
	void permute(sampleType* inputPairs) {
		permute((complexValue*)inputPairs);
	}
	void permute(complexValue* input) {
		permute(scratchCopy(input), input);
	}
	void ipermute(sampleType* inputPairs, sampleType* outputPairs) {
		ipermute((complexValue*)inputPairs, (complexValue*)outputPairs);
	}
	void ipermute(sampleType* inputPairs) {
		ipermute((complexValue*)inputPairs);
	}
	void ipermute(complexValue* input) {
		ipermute(scratchCopy(input), input);
	}

	// Finds next-highest value of form: 2^N * (1, 3 or 9)
	static int findGoodSize(int N) {
		int power2 = 1;
		while (N > 1) {
			if (N == 3 || N == 9) break;

			N = (N/2) + (N%2); // Divide by 2, rounding up
			power2 *= 2;
		}
		return power2*N;
	}
};

template <typename sampleType>
std::shared_ptr<FFT<sampleType>> getFft(int N, int parallel=1, int cpuFeatures=-1);

namespace _internal {

template<typename complexValue, typename complexTwiddle>
INLINE complexValue applyTwiddle(const complexValue& value, const complexTwiddle& twiddle) {
	return {
		.real = value.real*twiddle.real - value.imag*twiddle.imag,
		.imag = value.imag*twiddle.real + value.real*twiddle.imag
	};
}
template<typename complexTwiddle>
INLINE complexTwiddle conjugateTwiddle(const complexTwiddle& twiddle) {
	return {.real=twiddle.real, .imag=-twiddle.imag};
}

template<typename complexValue, typename complexTwiddle>
INLINE void fft3(complexValue& offsetBuffer0, complexValue& offsetBuffer1, complexValue& offsetBuffer2, const complexValue& A, const complexValue& B, const complexValue& C) {
	const complexTwiddle factor3 = {-0.5, -0.8660254037844386};
	offsetBuffer0.real = A.real + B.real + C.real;
	offsetBuffer0.imag = A.imag + B.imag + C.imag;
	offsetBuffer1.real = A.real + (B.real + C.real)*factor3.real + (C.imag - B.imag)*factor3.imag;
	offsetBuffer1.imag = A.imag + (B.real - C.real)*factor3.imag + (B.imag + C.imag)*factor3.real;
	offsetBuffer2.real = A.real + (B.real + C.real)*factor3.real + (B.imag - C.imag)*factor3.imag;
	offsetBuffer2.imag = A.imag + (C.real - B.real)*factor3.imag + (B.imag + C.imag)*factor3.real;
}

template<typename complexValue>
INLINE void fft4(complexValue& offsetBuffer0, complexValue& offsetBuffer1, complexValue& offsetBuffer2, complexValue& offsetBuffer3, const complexValue& A, const complexValue& B, const complexValue& C, const complexValue& D) {
	{
		complexValue ACsum = {A.real + C.real, A.imag + C.imag};
		complexValue BDsum = {B.real + D.real, B.imag + D.imag};
		offsetBuffer0.real = ACsum.real + BDsum.real;
		offsetBuffer0.imag = ACsum.imag + BDsum.imag;
		offsetBuffer2.real = ACsum.real - BDsum.real;
		offsetBuffer2.imag = ACsum.imag - BDsum.imag;
	}
	{
		complexValue ACdiff = {A.real - C.real, A.imag - C.imag};
		complexValue swappedBDdiff = {B.imag - D.imag, D.real - B.real};
		offsetBuffer1.real = ACdiff.real + swappedBDdiff.real;
		offsetBuffer1.imag = ACdiff.imag + swappedBDdiff.imag;
		offsetBuffer3.real = ACdiff.real - swappedBDdiff.real;
		offsetBuffer3.imag = ACdiff.imag - swappedBDdiff.imag;
	}
}

/* COBRA permutation pattern
from Carter and Gatlin's "Towards an Optimal Bit-Reversal Permutation Program" <http://www.cs.technion.ac.il/~itai/Courses/Cache/bit.pdf>

Item-by-item permutation has a poorly-localised memory access pattern, and therefore poor performance for larger sizes.  COBRA is an approach to mitigate this.

We first separate the input indexes using three factors, where A*B*C = N:
	index = a + b*A + c*B*A
(where 0 <= a < A, etc.)

Because our permutation function is based on factor-swapping, for suitable A, B, C:
	permute(index) = permute(a + b*A + c*B*A)
		= permute(a) + permute(b*A) + permute(c*B*A)
and also:
	0 <= permute(c*B*A) < C

This means that for fixed "b", iterating through "a" gives us a nice access pattern on our input, and iterating through "c2 = permute(c*B*A)" gives us a nice access pattern on our output.

We use "transfer" as an (A x C) table, picking A and C so that they are roughly equal, and so the table is small enough that it's likely to fit close in the cache (L1 or L2).
Since it's entirely cached, random-access is less expensive, so we can copy in row-by-row for the input (reading A contiguous values at a time), and then copy out column-by-column for the output (writing C contiguous values at a time).
*/
template <typename complexValue, bool inverse, int fixedParallel>
void cobraPermute(complexValue* input, complexValue* output, int N, int parallel, const std::vector<int>& factors, const int* indexes) {
#define PARALLEL (fixedParallel ? fixedParallel : parallel)
	if (!fixedParallel && parallel == 1) return cobraPermute<complexValue,inverse,1>(input, output, N, parallel, factors, indexes);

	int transferSize = 16384/sizeof(complexValue)/PARALLEL; // don't take up more than 16kB
	if (N <= transferSize || transferSize <= 9) {
		// Small enough to do directly
		for (int i = 0; i < N; i++) {
			if (inverse) {
				for (int p = 0; p < PARALLEL; p++) {
					output[i*PARALLEL + p] = input[indexes[i]*PARALLEL + p];
				}
			} else {
				for (int p = 0; p < PARALLEL; p++) {
					output[indexes[i]*PARALLEL + p] = input[i*PARALLEL + p];
				}
			}
		}
		return;
	}

	int sqrtSize = (int)sqrt(transferSize);
	int A = N, B = 1, C = 1;
	unsigned int factorIndex = 0;
	while (factorIndex < factors.size() && 2*C < sqrtSize && 2*C < A) {
		int f = factors[factorIndex++];
		A /= f;
		C *= f;
	}
	while (factorIndex < factors.size() && A*C > transferSize) {
		int f = factors[factorIndex++];
		A /= f;
		B *= f;
	}
	transferSize = A*C; // It might be smaller than our original goal
	complexValue* transfer = new complexValue[transferSize*PARALLEL];

	if (inverse) {
		for (int b = 0; b < B; b++) {
			for (int a = 0; a < A; a++) {
				int offset = indexes[a + b*A];
				for (int c2 = 0; c2 < C; c2++) {
					int transferIndex = (a + c2*A)*PARALLEL, inputIndex = (offset + c2)*PARALLEL;
					for (int p = 0; p < PARALLEL; p++) {
						transfer[transferIndex + p] = input[inputIndex + p];
					}
				}
			}
			for (int c = 0; c < C; c++) {
				int c2 = indexes[c*B*A];
				int offset = b*A + c*B*A;
				int outputIndex = (offset)*PARALLEL, transferIndex = (c2*A)*PARALLEL;
				// Copy a contiguous block to the output
				std::memcpy(output + outputIndex, transfer + transferIndex, A*PARALLEL*sizeof(complexValue));
			}
		}
	} else {
		for (int b = 0; b < B; b++) {
			for (int c = 0; c < C; c++) {
				int c2 = indexes[c*B*A];
				int offset = b*A + c*B*A;
				int transferIndex = (c2*A)*PARALLEL, inputIndex = (offset)*PARALLEL;
				// Copy a contiguous block from the input
				std::memcpy(transfer + transferIndex, input + inputIndex, A*PARALLEL*sizeof(complexValue));
			}
			for (int a = 0; a < A; a++) {
				int offset = indexes[a + b*A];
				for (int c2 = 0; c2 < C; c2++) {
					int outputIndex = (offset + c2)*PARALLEL, transferIndex = (a + c2*A)*PARALLEL;
					for (int p = 0; p < PARALLEL; p++) {
						output[outputIndex + p] = transfer[transferIndex + p];
					}
				}
			}
		}
	}
	delete[] transfer;
#undef PARALLEL
}

// The main implementation
template <typename sampleType, typename twiddleType>
class BasicFFT {
	using complexValue = complex<sampleType>;
	using complexTwiddle = complex<twiddleType>;
private:

	class PlanStep {
	public:
		int N;
		int parallel;
		int repeats;
		bool isFactorised;
		std::vector<complexTwiddle> twiddles;
	};
	class Plan {
	public:
		std::vector<PlanStep> steps;
		std::vector<int> permutation;
		std::vector<int> factors;
	};

	int N = 0;
	int parallelN = 1;
	Plan plan;

	// Defines our strategy - we always factor out 4s first
	static int getFactor(int N) {
		if (N > 4) {
			if (!(N&3)) return 4;
			if (!(N%3)) return 3;
		}
		return 0;
	}

	static Plan getPlan(int N, int parallel) {
		int splitFactor = getFactor(N);
		if (splitFactor) {
			int splitN = N/splitFactor;
			Plan plan = getPlan(splitN, parallel*splitFactor);

			PlanStep step = PlanStep{
				.N=splitFactor,
				.parallel=parallel,
				.repeats=splitN,
				.isFactorised=true
			};
			// Twiddles and permutations
			std::vector<int> newPermutation(N);
			step.twiddles.resize(splitN*(splitFactor - 1));
			for (int i = 0; i < splitN; i++) {
				int permutedI = plan.permutation[i];
				for (int k = 0; k < splitFactor; k++) {
					newPermutation[i*splitFactor + k] = permutedI + k*splitN;

					if (k > 0) { // The first twiddle factor is always 1, so we skip it
						double twiddlePhase = (permutedI*k)*(2*M_PI)/N;
						step.twiddles[i*(splitFactor - 1) + (k - 1)] = {
							.real=(twiddleType)cos(twiddlePhase),
							.imag=(twiddleType)-sin(twiddlePhase)
						};
					}
				}
			}

			plan.permutation = newPermutation;
			plan.steps.push_back(step);
			plan.factors.push_back(splitFactor);
			return plan;
		} else {
			Plan plan;
			plan.factors.push_back(N);
			plan.permutation.resize(N);
			for (int i = 0; i < N; i++) {
				plan.permutation[i] = i;
			}
			if (N > 1) {
				PlanStep step = PlanStep{
					.N=N,
					.parallel=parallel,
					.repeats=1,
					.isFactorised=false
				};
				plan.steps.push_back(step);
			}
			return plan;
		}
	}

	// Apply parallel interleaved 2-point FFTs
	static void fft2Step(complexValue* buffer, int fftStride, int parallel) {
		complexValue* offsetBuffer0 = buffer;
		complexValue* offsetBuffer1 = buffer + fftStride;
		for (int offset = 0; offset < parallel; offset++) {
			complexValue A = *offsetBuffer0;
			complexValue B = *offsetBuffer1;

			offsetBuffer0->real += B.real;
			offsetBuffer0->imag += B.imag;
			offsetBuffer1->real = A.real - B.real;
			offsetBuffer1->imag = A.imag - B.imag;

			++offsetBuffer0;
			++offsetBuffer1;
		}
	}

	// Apply parallel interleaved 3-point FFTs, repeated across the buffer
	template<bool inverse, bool applyTwiddles>
	static void fft3Step(complexValue* buffer, const complexTwiddle* twiddles, int fftStride, int parallel, int repeatCount, int repeatStride) {
		for (int repeat = 0; repeat < repeatCount; repeat++) {
			complexValue* offsetBuffer0 = buffer;
			complexValue* offsetBuffer1 = offsetBuffer0 + fftStride;
			complexValue* offsetBuffer2 = offsetBuffer1 + fftStride;
			buffer += repeatStride;
			const complexTwiddle twiddle1 = inverse ? conjugateTwiddle<complexTwiddle>(twiddles[0]) : twiddles[0];
			const complexTwiddle twiddle2 = inverse ? conjugateTwiddle<complexTwiddle>(twiddles[1]) : twiddles[1];
			twiddles += 2;

			complexValue A, B, C, A2, B2, C2;
			for (int offset = 0; offset < parallel; offset++) {
				A = *offsetBuffer0;
				B = (applyTwiddles && !inverse) ? applyTwiddle<complexValue, complexTwiddle>(*offsetBuffer1, twiddle1) : *offsetBuffer1;
				C = (applyTwiddles && !inverse) ? applyTwiddle<complexValue, complexTwiddle>(*offsetBuffer2, twiddle2) : *offsetBuffer2;

				if (inverse) {
					fft3<complexValue, complexTwiddle>(A2, B2, C2, A, C, B); // swap positions 1+2
				} else {
					fft3<complexValue, complexTwiddle>(A2, B2, C2, A, B, C);
				}

				*offsetBuffer0 = A2;
				*offsetBuffer1 = (applyTwiddles && inverse) ? applyTwiddle<complexValue, complexTwiddle>(B2, twiddle1) : B2;
				*offsetBuffer2 = (applyTwiddles && inverse) ? applyTwiddle<complexValue, complexTwiddle>(C2, twiddle2) : C2;

				++offsetBuffer0;
				++offsetBuffer1;
				++offsetBuffer2;
			}
		}
	}

	// Apply parallel interleaved 4-point FFTs, repeated across the buffer
	template<bool inverse, bool applyTwiddles, int fixedParallel=0>
	INLINE static void fft4Step(complexValue* buffer, const complexTwiddle* twiddles, int fftStride, int parallel, int repeatCount, int repeatStride) {
		if (!fixedParallel) {
			switch (parallel) {
				case 1: return fft4Step<inverse, applyTwiddles, 1>(buffer, twiddles, fftStride, parallel, repeatCount, repeatStride);
				case 2: return fft4Step<inverse, applyTwiddles, 2>(buffer, twiddles, fftStride, parallel, repeatCount, repeatStride);
				case 4: return fft4Step<inverse, applyTwiddles, 4>(buffer, twiddles, fftStride, parallel, repeatCount, repeatStride);
			}
		}
		for (int repeat = 0; repeat < repeatCount; repeat++) {
			complexValue* offsetBuffer0 = buffer;
			complexValue* offsetBuffer1 = offsetBuffer0 + fftStride;
			complexValue* offsetBuffer2 = offsetBuffer1 + fftStride;
			complexValue* offsetBuffer3 = offsetBuffer2 + fftStride;
			buffer += repeatStride;
			const complexTwiddle twiddle1 = inverse ? conjugateTwiddle<complexTwiddle>(twiddles[0]) : twiddles[0];
			const complexTwiddle twiddle2 = inverse ? conjugateTwiddle<complexTwiddle>(twiddles[1]) : twiddles[1];
			const complexTwiddle twiddle3 = inverse ? conjugateTwiddle<complexTwiddle>(twiddles[2]) : twiddles[2];
			twiddles += 3;

			complexValue A, B, C, D, A2, B2, C2, D2;
			for (int offset = 0; fixedParallel ? offset < fixedParallel : offset < parallel; offset++) {
				A = *offsetBuffer0;
				B = (applyTwiddles && !inverse) ? applyTwiddle<complexValue, complexTwiddle>(*offsetBuffer1, twiddle1) : *offsetBuffer1;
				C = (applyTwiddles && !inverse) ? applyTwiddle<complexValue, complexTwiddle>(*offsetBuffer2, twiddle2) : *offsetBuffer2;
				D = (applyTwiddles && !inverse) ? applyTwiddle<complexValue, complexTwiddle>(*offsetBuffer3, twiddle3) : *offsetBuffer3;

				if (inverse) {
					fft4<complexValue>(A2, B2, C2, D2, A, D, C, B); // swap positions 1+3
				} else {
					fft4<complexValue>(A2, B2, C2, D2, A, B, C, D);
				}

				*offsetBuffer0 = A2;
				*offsetBuffer1 = (applyTwiddles && inverse) ? applyTwiddle<complexValue, complexTwiddle>(B2, twiddle1) : B2;
				*offsetBuffer2 = (applyTwiddles && inverse) ? applyTwiddle<complexValue, complexTwiddle>(C2, twiddle2) : C2;
				*offsetBuffer3 = (applyTwiddles && inverse) ? applyTwiddle<complexValue, complexTwiddle>(D2, twiddle3) : D2;

				++offsetBuffer0;
				++offsetBuffer1;
				++offsetBuffer2;
				++offsetBuffer3;
			}
		}
	}

	template<bool inverse>
	void runPlan(complexValue* buffer) {
		if (N <= 1) return;

		auto iterator = inverse ? plan.steps.end() - 1 : plan.steps.begin();
		auto iteratorEnd = inverse ? plan.steps.begin() - 1 : plan.steps.end();

		for (; iterator != iteratorEnd; inverse ? --iterator : ++iterator) {
			PlanStep& step = *iterator;
			if (step.isFactorised) {
				if (step.N == 4) {
					fft4Step<inverse, true>(buffer, step.twiddles.data(), step.parallel, step.parallel, step.repeats, step.N*step.parallel);
				} else {
					fft3Step<inverse, true>(buffer, step.twiddles.data(), step.parallel, step.parallel, step.repeats, step.N*step.parallel);
				}
			} else {
				// Dummy twiddles to pass in - otherwise it might dereference a null pointer (particularly if optimisation is off)
				static constexpr complexTwiddle dummyTwiddles[3] = {{0, 0}, {0, 0}, {0, 0}};
				if (step.N == 4) {
					fft4Step<inverse, false>(buffer, dummyTwiddles, step.parallel, step.parallel, 1, 0);
				} else if (step.N == 3) {
					fft3Step<inverse, false>(buffer, dummyTwiddles, step.parallel, step.parallel, 1, 0);
				} else {
					fft2Step(buffer, step.parallel, step.parallel);
				}
			}
		}
	}

public:

	BasicFFT(int N, int parallel) {
		setSize(N, parallel);
	};

	~BasicFFT() {}

	const int* getPermutation() {
		return plan.permutation.data();
	}
	const std::vector<int>& getFactors() {
		return plan.factors;
	}

	int size() {
		return N;
	}

	int parallel() {
		return this->parallelN;
	}

	int setSize(int N, int parallel) {
		if (N != this->N) N = FFT<sampleType>::findGoodSize(N);
		if (N != this->N || parallel != this->parallelN) {
			this->N = N;
			this->parallelN = parallel;

			plan = getPlan(N, parallel);
		}
		return N;
	}

	void fft(complexValue* input) {
		runPlan<false>(input);
	};
	void ifft(complexValue* input) {
		runPlan<true>(input);
	};
};

// Wrapper from BasicFFT to FFT interface
template <typename sampleType>
class BasicWrapper : public FFT<sampleType> {
	BasicFFT<sampleType, sampleType> innerFft;
public:
	BasicWrapper(int N, int parallel) : innerFft(N, parallel) {}

	int size() override {return innerFft.size();};
	int parallel() override {return innerFft.parallel();};
	int setSize(int N, int parallel) override {return innerFft.setSize(N, parallel);};
	const std::vector<int>& getFactors() override {return innerFft.getFactors();};
	const int* getPermutation() override {return innerFft.getPermutation();};

	void fft(complex<sampleType>* inputPairs) override {innerFft.fft(inputPairs);};
	void ifft(complex<sampleType>* inputPairs) override {innerFft.ifft(inputPairs);};
};

} // namespace _internal

// Generic form of template - use BasicFFT with sampleType
template <typename sampleType>
std::shared_ptr<FFT<sampleType>> getFft(int N, int parallel, int cpuFeatures) {
	using namespace _internal;
	return std::make_shared<BasicWrapper<sampleType>>(N, parallel);
}

} // fft
} // signalsmith

#endif // SIGNALSMITH_FFT_H_
