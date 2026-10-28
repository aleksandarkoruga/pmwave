// PluginPMWave.hpp
// Aleksandar Koruga (aleksandar.koruga@gmail.com)

#pragma once
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include "SC_PlugIn.hpp"

// four lanes (a junction's rails, or four rails in the length update) map onto one SSE register; other targets
// (or -DPMWAVE_NO_SIMD) use the plain loops, which clang/gcc vectorise for NEON on their own
#if !defined(PMWAVE_NO_SIMD) && (defined(__SSE2__) || defined(_M_X64) || (defined(_M_IX86_FP) && _M_IX86_FP >= 2))
#include <emmintrin.h>
#define PMWAVE_SSE 1
#endif


namespace PMWave {

// build-time variants: grid size (power of two >= 2) and longest rail in ms, e.g. -DPMWAVE_GRID=2 -DPMWAVE_MAX_RAIL_MS=21
#ifndef PMWAVE_GRID
#define PMWAVE_GRID 4
#endif
#ifndef PMWAVE_MAX_RAIL_MS
#define PMWAVE_MAX_RAIL_MS 40
#endif

// junction grid (powers of two, the torus wraps with masks); junction J = i * NY + j
constexpr int NX = PMWAVE_GRID;
constexpr int NY = PMWAVE_GRID;
constexpr int NJ = NX * NY;
// every junction stores four rails (its west and north waveguides, both directions); rail R = J * 4 + lane
constexpr int NR = NJ * 4;
// rail length range in seconds; ring buffers are sized for the longest rail at the server's sample rate
constexpr double kMaxRailTime = PMWAVE_MAX_RAIL_MS * 0.001;
constexpr double kMinRailTime = 0.0001;
// the newest cubic tap sits one sample before the rail length, and must not reach the slot being written
constexpr float kMinRailSamples = 4.f;
// fastest change of a rail length in samples per sample (Doppler factor 0..2 on what is ringing)
constexpr float kMaxSlew = 1.f;
// loop filter coefficients are tabulated over log2(rail length)
constexpr int kLossTabSize = 256;
// slots of padding (16 bytes each) after every junction's ring buffer: one cache line
constexpr int kRailPad = 4;
// wander: per-rail smooth random modulation, new random target every 1/rate seconds
constexpr float kWanderRateMin = 0.1f;
constexpr float kWanderRateMax = 0.4f;

// 2^x for |x| < 60: nearest integer plus a degree-5 polynomial on [-0.5, 0.5] (relative error < 3e-6)
inline float exp2Fast(const float x)
{
	const int i = (int)(x + 64.5f) - 64;
	const float g = (x - i) * 0.69314718f;
	const float p = 1 + g * (1 + g * (0.5f + g * (1.f / 6 + g * (1.f / 24 + g * (1.f / 120)))));
	const uint32_t bits = (uint32_t)(i + 127) << 23;
	float scale;
	std::memcpy(&scale, &bits, sizeof(scale));
	return p * scale;
}

#ifdef PMWAVE_SSE
// same arithmetic, four values at a time
inline __m128 exp2Fast(const __m128 x)
{
	const __m128i i = _mm_sub_epi32(_mm_cvttps_epi32(_mm_add_ps(x, _mm_set1_ps(64.5f))), _mm_set1_epi32(64));
	const __m128 g = _mm_mul_ps(_mm_sub_ps(x, _mm_cvtepi32_ps(i)), _mm_set1_ps(0.69314718f));
	__m128 p = _mm_add_ps(_mm_set1_ps(1.f / 24), _mm_mul_ps(g, _mm_set1_ps(1.f / 120)));
	p = _mm_add_ps(_mm_set1_ps(1.f / 6), _mm_mul_ps(g, p));
	p = _mm_add_ps(_mm_set1_ps(0.5f), _mm_mul_ps(g, p));
	p = _mm_add_ps(_mm_set1_ps(1.f), _mm_mul_ps(g, p));
	p = _mm_add_ps(_mm_set1_ps(1.f), _mm_mul_ps(g, p));
	return _mm_mul_ps(p, _mm_castsi128_ps(_mm_slli_epi32(_mm_add_epi32(i, _mm_set1_epi32(127)), 23)));
}
#endif

class PMWave : public SCUnit {
public:
    PMWave();
    ~PMWave();

private:
    // Calc function
    void next(int nSamples);
	std::array<float, 2> _render(const float _in);

	// Junction J stores the rails of its west and north waveguides, both directions, interleaved. Lanes: 0 = W->J,
	// 1 = N->J, 2 = J->W, 3 = J->N. Every rail has its own length, glide time and wander, so all 64 delay lines
	// differ and the two directions of a waveguide are not the same length.
	// Junction-major ([NJ][_bufLength + kRailPad][4], RTAlloc'd): a fetched cache line serves four samples of a
	// junction (these buffers are far larger than the caches once several plates run). The pad staggers the
	// junctions across cache sets: unpadded they sit a power of two apart and every write slot lands in the same set.
	float* _rails;
	int _bufLength;
	int _mask;
	int _railStride; // floats from one junction's buffer to the next
	float _lp[NJ][4];   // loop filter state of each rail
	float _wave[NJ][4]; // filtered waves arriving this sample, same lanes
	float _p[NJ];       // junction pressures

	// Rail lengths live in log2(samples): `time` maps linearly onto that range, and the glide is a one-pole lag there,
	// so equal steps are equal musical intervals and every glide takes about the same time. A rail reads
	// 2^(target - remaining lag + wander) samples, computed as targetLength * 2^offset: the offset is small near the
	// target, so float keeps the read position smooth to the last bit of the length.
	float _logMin;
	float _logMax;
	float _minLength;
	float _maxLength;
	float _logTarget[NR];    // time + rand1 offset, clipped to the range
	float _targetLength[NR]; // 2^_logTarget
	// target - lagged length (without wander). Kept as the remaining distance so the lag decays in full precision:
	// stepping the length itself stalls once (target - length) * coefficient drops below its float resolution
	float _logRemain[NR];
	float _lagCoef[NR];   // lag: fraction of the remaining distance per sample; constant rate: octaves per sample
	bool _constantRate;   // glideMode 1: every rail moves at a steady rate in log length (glide = seconds per octave)
	bool _retarget;       // targets moved: recompute the read coefficients even if no lag is running
	float _logRead[NR];   // log2 of the length being read, for the loss tables
	float _length[NR];    // length being read, in samples

	// read coefficients of each rail
	int _tapDelay[NR];  // integer part of the length; the four taps sit at delays d-1 .. d+2
	float _tap[4][NR];  // 4-point Lagrange weights
	float _gain[NR];
	float _pole[NR];

	// loop filter coefficients over log2(rail length): y = gain * x + pole * y[n-1]
	float _gainTab[kLossTabSize + 1];
	float _poleTab[kLossTabSize + 1];
	float _lossScale; // table entries per octave

	float _rand1[NR];
	float _rand2[NR];
	// wander: smoothstep between random values in [-1, 1]
	float _wanderPhase[NR];
	float _wanderInc[NR];
	float _wanderFrom[NR];
	float _wanderTo[NR];

	int w;
	float _time;
	float _glide;    // 60 dB lag time in seconds (like Lag)
	float _decay;    // T60 in seconds at DC (raw input value, clipped in _updateLoss)
	float _hfDecay;  // T60 in seconds at Nyquist
	float _randAmt1; // rail length spread, ratio: lengths within x/÷ (1 + rand1)
	float _spread1;  // 2 * log2(1 + rand1), octaves per unit of _rand1
	float _randAmt2; // glide time spread, ratio: times within x/÷ (1 + rand2)
	float _wander;   // wander depth, ratio
	float _wanderOctaves;


	// NaN-safe (fmaxf(NaN, min) = min), for inputs; MSVC calls into the CRT for these, so the per-sample paths below
	// use _clamp on values that are already finite
	float _fclipf(const float _in, const float _min, const float _max)
	{
		return fminf(fmaxf(_in, _min), _max);
	}
	static float _clamp(const float x, const float lo, const float hi) { return x < lo ? lo : (x > hi ? hi : x); }

	void _updateTargets()
	{
		const float logTime = _logMax + _fclipf(_time, 0, 1) * (_logMin - _logMax);
		for (int R = 0; R < NR; ++R)
		{
			const float target = _fclipf(logTime + _spread1 * _rand1[R], _logMin, _logMax);
			_logRemain[R] += target - _logTarget[R]; // the lagged length stays where it is
			_logTarget[R] = target;
			_targetLength[R] = exp2Fast(target);
		}
		_retarget = true;
	}

	void _updateLag()
	{
		const float spread = 2 * std::log2(1 + _randAmt2);
		for (int R = 0; R < NR; ++R)
		{
			const double seconds = _glide * std::exp2(spread * _rand2[R]);
			if (_constantRate)
				_lagCoef[R] = seconds > 0 ? (float)(1 / (seconds * mRate->mSampleRate)) : 1e30f;
			else // 1 - 0.001^(1 / samples): the remaining interval shrinks by 60 dB in `seconds`
				_lagCoef[R] = seconds > 0 ? (float)-std::expm1(-6.907755278982137 / (seconds * mRate->mSampleRate)) : 1.f;
		}
	}

	float _wanderValue(const int R) const
	{
		const float ph = _wanderPhase[R];
		return _wanderFrom[R] + (_wanderTo[R] - _wanderFrom[R]) * (ph * ph * (3 - 2 * ph));
	}

	// Per-rail loss (Jot): a rail D samples long gets a one-pole loop filter g(1-b)/(1-b z^-1) with
	// |H(0)| = 10^(-3D/(fs*decay)) and |H(pi)| = 10^(-3D/(fs*hfDecay)), so T60 does not depend on plate size.
	void _updateLoss(const float decay, const float hfDecay)
	{
		_decay = decay;
		_hfDecay = hfDecay;
		const double kLo = -6.907755278982137 / (mRate->mSampleRate * _fclipf(decay, 0.01f, 1000.f));
		const double kHi = -6.907755278982137 / (mRate->mSampleRate * _fclipf(hfDecay, 0.001f, 1000.f));
		for (int i = 0; i <= kLossTabSize; ++i)
		{
			const double D = std::exp2(_logMin + i / (double)_lossScale); // rail delay at entry i
			const double r = std::exp((kHi - kLo) * D); // Nyquist / DC gain ratio
			const double b = (1 - r) / (1 + r);
			_gainTab[i] = (float)(std::exp(kLo * D) * (1 - b));
			_poleTab[i] = (float)b;
		}
		_updateGains();
	}

	void _updateGains()
	{
		for (int R = 0; R < NR; ++R)
		{
			const float x = (_logRead[R] - _logMin) * _lossScale;
			const int i = x < 0 ? 0 : (x >= kLossTabSize ? kLossTabSize - 1 : (int)x);
			const float fx = _clamp(x - i, 0, 1);
			_gain[R] = _gainTab[i] + (_gainTab[i + 1] - _gainTab[i]) * fx;
			_pole[R] = _poleTab[i] + (_poleTab[i + 1] - _poleTab[i]) * fx;
		}
	}

	// read lengths of all rails from their log2 offsets to the target, then the read coefficients
	void _updateLengths(const float* offset)
	{
		float D[NR];
		for (int R = 0; R < NR; ++R) _logRead[R] = _clamp(_logTarget[R] + offset[R], _logMin, _logMax);
#ifdef PMWAVE_SSE
		for (int R = 0; R < NR; R += 4)
		{
			const __m128 d = _mm_mul_ps(_mm_loadu_ps(_targetLength + R), exp2Fast(_mm_loadu_ps(offset + R)));
			_mm_storeu_ps(D + R, _mm_min_ps(_mm_max_ps(d, _mm_set1_ps(_minLength)), _mm_set1_ps(_maxLength)));
		}
#else
		for (int R = 0; R < NR; ++R) D[R] = _clamp(_targetLength[R] * exp2Fast(offset[R]), _minLength, _maxLength);
#endif
		// a read point moving faster than one sample per sample would overtake time itself (only very short glides)
		for (int R = 0; R < NR; ++R)
		{
			const float prev = _length[R];
			if (std::abs(D[R] - prev) > kMaxSlew)
			{
				D[R] = D[R] > prev ? prev + kMaxSlew : prev - kMaxSlew;
				_logRead[R] = std::log2(D[R]);
				// move the lag to where the rail actually is (offset = wander - remaining lag)
				_logRemain[R] += offset[R] - (_logRead[R] - _logTarget[R]);
			}
		}

		// 4-point Lagrange through delays d-1 .. d+2: flat to much higher frequencies than linear interpolation,
		// so the high-frequency decay barely depends on the fractional part of the length
#ifdef PMWAVE_SSE
		for (int R = 0; R < NR; R += 4)
		{
			const __m128 d = _mm_loadu_ps(D + R);
			const __m128i di = _mm_cvttps_epi32(d);
			const __m128 f = _mm_sub_ps(d, _mm_cvtepi32_ps(di));
			const __m128 fp1 = _mm_add_ps(f, _mm_set1_ps(1.f));
			const __m128 fm1 = _mm_sub_ps(f, _mm_set1_ps(1.f));
			const __m128 fm2 = _mm_sub_ps(f, _mm_set1_ps(2.f));
			const __m128 a = _mm_mul_ps(f, fm1);   // f (f-1)
			const __m128 c = _mm_mul_ps(fp1, fm2); // (f+1)(f-2)
			_mm_storeu_ps(_tap[0] + R, _mm_mul_ps(_mm_mul_ps(a, fm2), _mm_set1_ps(-1.f / 6)));
			_mm_storeu_ps(_tap[1] + R, _mm_mul_ps(_mm_mul_ps(c, fm1), _mm_set1_ps(0.5f)));
			_mm_storeu_ps(_tap[2] + R, _mm_mul_ps(_mm_mul_ps(c, f), _mm_set1_ps(-0.5f)));
			_mm_storeu_ps(_tap[3] + R, _mm_mul_ps(_mm_mul_ps(a, fp1), _mm_set1_ps(1.f / 6)));
			_mm_storeu_si128((__m128i*)(_tapDelay + R), di);
			_mm_storeu_ps(_length + R, d);
		}
#else
		for (int R = 0; R < NR; ++R)
		{
			const int d = (int)D[R];
			const float f = D[R] - d;
			const float a = f * (f - 1), c = (f + 1) * (f - 2);
			_tap[0][R] = a * (f - 2) * (-1.f / 6);
			_tap[1][R] = c * (f - 1) * 0.5f;
			_tap[2][R] = c * f * -0.5f;
			_tap[3][R] = a * (f + 1) * (1.f / 6);
			_tapDelay[R] = d;
			_length[R] = D[R];
		}
#endif
		_updateGains();
	}

};

} // namespace PMWave
