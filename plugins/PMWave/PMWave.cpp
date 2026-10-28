// PluginPMWave.cpp
// Aleksandar Koruga (aleksandar.koruga@gmail.com)

#include "SC_PlugIn.hpp"
#include "PMWave.hpp"

static InterfaceTable* ft;

namespace PMWave {

// excitation and pickups; the input sits at (0,0)
constexpr int kInput = 0;
// (1, NY/2) and (0, NY/2) lie at odd/even hop distance from the input, so their arrivals interleave and the channels
// are uncorrelated. On the 4x4 grid that is (1,2) and (0,2), at levels within ~2 dB: mirror pairs like (x,y)/(y,x)
// are near-identical, and (2,2) (antipodal focus) and (0,0) (input) are ~8 dB louder than every other junction.
constexpr int kPickupL = 1 * NY + NY / 2;
constexpr int kPickupR = 0 * NY + NY / 2;

PMWave::PMWave() {
	auto unit = this; // for ClearUnitIfMemFailed

	// 0	1				2					3			4				5		6		7		8
	//arg in, decay (T60 s), time (size, 0..1), glide (s), hfDecay (T60 s), rand1, rand2, wander, glideMode
	w = 0;
	const double fs = mRate->mSampleRate;
	const double longest = kMaxRailTime * fs;
	_bufLength = 1;
	while (_bufLength < longest + 4) _bufLength <<= 1;
	_mask = _bufLength - 1;
	_logMax = (float)std::log2(longest);
	_logMin = (float)std::log2(std::fmax(kMinRailTime * fs, (double)kMinRailSamples));
	_maxLength = exp2Fast(_logMax);
	_minLength = exp2Fast(_logMin);
	_lossScale = kLossTabSize / (_logMax - _logMin);

	_railStride = (_bufLength + kRailPad) * 4;
	_rails = (float*)RTAlloc(mWorld, sizeof(float) * NJ * _railStride);
	ClearUnitIfMemFailed(_rails);
	std::memset(_rails, 0, sizeof(float) * NJ * _railStride);
	std::memset(_lp, 0, sizeof(_lp));
	std::memset(_wave, 0, sizeof(_wave));
	std::memset(_p, 0, sizeof(_p));

	_time = in0(2);
	_glide = fmaxf(in0(3), 0.0f);
	_randAmt1 = fmaxf(in0(5), 0.0f);
	_spread1 = 2 * std::log2(1 + _randAmt1);
	_randAmt2 = fmaxf(in0(6), 0.0f);
	_wander = fmaxf(in0(7), 0.0f);
	_wanderOctaves = std::log2(1 + _wander);
	_constantRate = in0(8) > 0.5f;

	// per-synth generator: reproducible with RandSeed and safe on the audio thread
	RGen& rgen = *mParent->mRGen;
	for (int R = 0; R < NR; ++R)
	{
		_rand1[R] = rgen.frand() - 0.5f;
		_rand2[R] = rgen.frand() - 0.5f;
		_wanderPhase[R] = rgen.frand();
		_wanderInc[R] = (float)((kWanderRateMin + (kWanderRateMax - kWanderRateMin) * rgen.frand()) / fs);
		_wanderFrom[R] = rgen.frand2();
		_wanderTo[R] = rgen.frand2();
		_logTarget[R] = _logRemain[R] = 0.f;
	}
	_updateTargets();
	_updateLag();

	// the plate starts at its target size, wander included
	float offset[NR];
	for (int R = 0; R < NR; ++R)
	{
		_logRemain[R] = 0.f;
		offset[R] = _wanderOctaves * _wanderValue(R);
		_logRead[R] = _logTarget[R];
		_length[R] = _fclipf(_targetLength[R] * exp2Fast(offset[R]), _minLength, _maxLength);
	}
	_updateLoss(in0(1), in0(4)); // loss tables
	_updateLengths(offset);      // read coefficients
	_retarget = false;

	mCalcFunc = make_calc_function<PMWave, &PMWave::next>();
	// the plate starts silent; running next(1) here would advance the state and inject in(0)[0] twice
	out0(0) = 0.f;
	out0(1) = 0.f;
}

PMWave::~PMWave() {
	if (_rails) RTFree(mWorld, _rails);
}

void PMWave::next(int nSamples) {

	// scalar and control-rate inputs point at a single float, so only audio-rate inputs may be indexed per sample
	const float* ins[8];
	int step[8];
	for (int k = 0; k < 8; ++k)
	{
		ins[k] = in(k);
		step[k] = isAudioRateIn(k) ? 1 : 0;
	}

	// inputs whose coefficients cost exp() are read once per block
	if (in0(1) != _decay || in0(4) != _hfDecay) _updateLoss(in0(1), in0(4));
	const float glide = fmaxf(in0(3), 0.0f);
	const float randAmt2 = fmaxf(in0(6), 0.0f);
	const bool constantRate = in0(8) > 0.5f;
	if (glide != _glide || randAmt2 != _randAmt2 || constantRate != _constantRate)
	{
		_glide = glide;
		_randAmt2 = randAmt2;
		_constantRate = constantRate;
		_updateLag();
	}

	float* outL = out(0);
	float* outR = out(1);
	for (int cnt = 0; cnt < nSamples; ++cnt)
	{
		const float time = ins[2][cnt * step[2]];
		const float randAmt1 = fmaxf(ins[5][cnt * step[5]], 0.0f);
		if (randAmt1 != _randAmt1)
		{
			_randAmt1 = randAmt1;
			_spread1 = 2 * std::log2(1 + randAmt1);
			_time = time;
			_updateTargets();
		}
		else if (time != _time)
		{
			_time = time;
			_updateTargets();
		}

		const float wander = fmaxf(ins[7][cnt * step[7]], 0.0f);
		if (wander != _wander)
		{
			_wander = wander;
			_wanderOctaves = std::log2(1 + wander);
		}

		const auto res = _render(ins[0][cnt * step[0]]);
		outL[cnt] = res[0];
		outR[cnt] = res[1];
	}
}


std::array<float, 2> PMWave::_render(const float _in)
{
	//GLISSANDOS BETWEEN JUNCTIONS: lag in log2(length) plus wander, every junction first, so all rails below read
	//lengths of the same instant
	bool moving = _retarget || _wanderOctaves > 0.f;
	if (_constantRate)
	{
		for (int R = 0; R < NR; ++R)
		{
			const float remain = _logRemain[R];
			if (remain != 0.f)
			{
				const float step = _lagCoef[R];
				_logRemain[R] = remain > step ? remain - step : (remain < -step ? remain + step : 0.f);
				moving = true;
			}
		}
	}
	else
	{
		for (int R = 0; R < NR; ++R)
		{
			const float remain = _logRemain[R];
			if (remain != 0.f)
			{
				const float left = remain - remain * _lagCoef[R];
				_logRemain[R] = std::abs(left) < 1e-6f ? 0.f : left;
				moving = true;
			}
		}
	}
	if (moving)
	{
		float offset[NR];
		if (_wanderOctaves > 0.f)
		{
			for (int R = 0; R < NR; ++R)
			{
				float ph = _wanderPhase[R] + _wanderInc[R];
				if (ph >= 1.f)
				{
					ph -= 1.f;
					_wanderFrom[R] = _wanderTo[R];
					_wanderTo[R] = mParent->mRGen->frand2();
				}
				_wanderPhase[R] = ph;
				offset[R] = _wanderOctaves * _wanderValue(R) - _logRemain[R];
			}
		}
		else
		{
			for (int R = 0; R < NR; ++R) offset[R] = -_logRemain[R];
		}
		_updateLengths(offset);
		_retarget = false;
	}

	// read the four rails stored at each junction, each at its own length: cubic interpolation, then the per-rail
	// loop filter. A sample written at time t sits in slot t - 1, so delay k is slot w - 1 - k.
	for (int J = 0; J < NJ; ++J)
	{
		const float* rails = _rails + J * _railStride;
		const int R = J * 4;
		float a0[4], a1[4], a2[4], a3[4];
		for (int k = 0; k < 4; ++k)
		{
			const int s0 = w - _tapDelay[R + k]; // delay d - 1
			a0[k] = rails[(s0 & _mask) * 4 + k];
			a1[k] = rails[((s0 - 1) & _mask) * 4 + k];
			a2[k] = rails[((s0 - 2) & _mask) * 4 + k];
			a3[k] = rails[((s0 - 3) & _mask) * 4 + k];
		}
#ifdef PMWAVE_SSE
		__m128 x = _mm_mul_ps(_mm_loadu_ps(_tap[0] + R), _mm_loadu_ps(a0));
		x = _mm_add_ps(x, _mm_mul_ps(_mm_loadu_ps(_tap[1] + R), _mm_loadu_ps(a1)));
		x = _mm_add_ps(x, _mm_mul_ps(_mm_loadu_ps(_tap[2] + R), _mm_loadu_ps(a2)));
		x = _mm_add_ps(x, _mm_mul_ps(_mm_loadu_ps(_tap[3] + R), _mm_loadu_ps(a3)));
		const __m128 y = _mm_add_ps(_mm_mul_ps(_mm_loadu_ps(_gain + R), x), _mm_mul_ps(_mm_loadu_ps(_pole + R), _mm_loadu_ps(_lp[J])));
		_mm_storeu_ps(_lp[J], y);
		_mm_storeu_ps(_wave[J], y);
#else
		for (int k = 0; k < 4; ++k)
		{
			const float x = _tap[0][R + k] * a0[k] + _tap[1][R + k] * a1[k] + _tap[2][R + k] * a2[k] + _tap[3][R + k] * a3[k];
			_wave[J][k] = _lp[J][k] = _gain[R + k] * x + _pole[R + k] * _lp[J][k];
		}
#endif
	}

	// junction pressures: half the sum of the arriving waves (from W and N on own rails, from E and S on theirs)
	for (int i = 0; i < NX; ++i)
		for (int j = 0; j < NY; ++j)
		{
			const int J = i * NY + j, E = i * NY + ((j + 1) & (NY - 1)), S = ((i + 1) & (NX - 1)) * NY + j;
			_p[J] = 0.5f * (_wave[J][0] + _wave[J][1] + _wave[E][2] + _wave[S][3]);
		}
	// the input arrives through a matched fifth port (2/5). The junction's own waves keep their 1/2 scaling
	// (the former junction gain fixed at 1.25), so the excitation point does not absorb energy.
	_p[kInput] += 0.4f * _in;

	// lossless scattering into ring memory: outgoing wave = junction pressure - wave that arrived from that side
	const int slot = ((w - 1) & _mask) * 4;
	for (int i = 0; i < NX; ++i)
		for (int j = 0; j < NY; ++j)
		{
			const int J = i * NY + j, W = i * NY + ((j - 1) & (NY - 1)), N = ((i - 1) & (NX - 1)) * NY + j;
			float* o = _rails + J * _railStride + slot;
#ifdef PMWAVE_SSE
			const __m128 v = _mm_loadu_ps(_wave[J]);
			_mm_storeu_ps(o, _mm_sub_ps(_mm_setr_ps(_p[W], _p[N], _p[J], _p[J]), _mm_shuffle_ps(v, v, _MM_SHUFFLE(1, 0, 3, 2))));
#else
			o[0] = _p[W] - _wave[J][2];
			o[1] = _p[N] - _wave[J][3];
			o[2] = _p[J] - _wave[J][0];
			o[3] = _p[J] - _wave[J][1];
#endif
		}

	w = (w + 1) & _mask;
	return {{ _p[kPickupL], _p[kPickupR] }};
}


} // namespace PMWave


PluginLoad(PMWaveUGens) {
	// Plugin magic
	ft = inTable;
	registerUnit<PMWave::PMWave>(ft, "PMWave", false);
}
