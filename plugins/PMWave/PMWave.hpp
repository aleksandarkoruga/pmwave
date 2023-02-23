// PluginPMWave.hpp
// Aleksandar Koruga (aleksandar.koruga@gmail.com)

#pragma once
#include <vector>
#include "SC_PlugIn.hpp"
#define BUFLENGTH 256
#define NXY 4


namespace PMWave {

class PMWave : public SCUnit {
public:
    PMWave();
	std::vector<float> _render(const float _in);

    // Destructor
    // ~PMWave();

private:
    // Calc function
    void next(int nSamples);

	float m_junction;	
    // Member variables
	int _nJunctionsX;
	int _nJunctionsY;
	int n_samps;

	//initialize global vectors
	float _currentDistance[NXY][NXY];
	float _meanLength[NXY][NXY];
	float sens[NXY][NXY];
	float rsens[NXY][NXY];


	float _p0[NXY][NXY][BUFLENGTH];
	float _p1[NXY][NXY][BUFLENGTH];
	float _p2[NXY][NXY][BUFLENGTH];
	float _p3[NXY][NXY][BUFLENGTH];


	float _del[NXY];

	int _dwn[NXY][NXY];
	int _ddup[NXY][NXY];
	float _fractional[NXY][NXY];
	float _rand1[NXY][NXY];
	float _rand2[NXY][NXY];
	float _p[NXY][NXY];


	float _time;//=50;// TIME to n
	float _damping;//=0.99999;
	float _cutoff;//=0.2;					
	float _speed;//=0.8f;


	float  ddiv;//=0.5;
	int w;//=0;						                       
	float envch;//=0;
	int   in_x;//=0; // where to write x 
	int   in_y;//=0; //  where to write y 
	float envd;//=0.03; 
	int _x;//=2; //where to read x
	int _y;//=2; //where to read y
	int _xOld;//=0;
	int _yOld;//=1;

	float _randAmt1;//=0.;
	float _randAmt2;//=0.;



//void _setDamping(const float _damp);

//void _updateTimeMeanLength(const float _inTime);
//void _updateSpeed(const float _inSpeed);
//float _render(const float _in);
//void _updateRandAmt1(const float _randAmt);
//void _updateRandAmt2(const float _randAmt);
//void _resetRandom();
//void _setCutoff(const float _cut);
	float getInvSampleRate() { return (float)(1 / mRate->mSampleRate); }



	void _updateTimeMeanLength(const float _inTime)
	{


		_time = n * _fclipf(_inTime, 0, 1);

		for (int i = 0; i < _nJunctionsX; ++i)
		{
			for (int j = 0; j < _nJunctionsY; ++j)
			{
				_meanLength[i][j] = _fclipf((_time + (2 * _rand1[i][j] * _randAmt1)), 1, n);
			}

		}

	}

	void _updateSpeed(const float _inSpeed)
	{
		_speed = _fclipf(_inSpeed, 0, 1);
		for (int i = 0; i < _nJunctionsX; ++i)
		{
			for (int j = 0; j < _nJunctionsY; ++j)
			{
				sens[i][j] = _fclipf(_speed + (2.0f * _rand2[i][j] * _randAmt2), getInvSampleRate(), 1.0f);
				rsens[i][j] = sens[i][j];

			}

		}

	}



	void _updateRandAmt1(const float _randAmt)
	{
		_randAmt1 = _randAmt;

		for (int i = 0; i < _nJunctionsX; ++i)
		{
			for (int j = 0; j < _nJunctionsY; ++j)
			{
				_meanLength[i][j] = _fclipf(_time + (2 * _rand1[i][j] * _randAmt1), 1, (float)n);
			}
		}


	}

	void _updateRandAmt2(const float _randAmt)
	{
		_randAmt2 = _randAmt;
		for (int i = 0; i < _nJunctionsX; ++i)
		{
			for (int j = 0; j < _nJunctionsY; ++j)
			{
				sens[i][j] = _fclipf(_speed + (2 * _rand2[i][j] * _randAmt2), getInvSampleRate(), 1.0f);
				rsens[i][j] = sens[i][j];
			}

		}


	}
	void _setDamping(const float _damp) { _damping = _damp; }
	void _setCutoff(const float _cut) { _cutoff = _cut; }

	void _resetRandom() {
		//((rand() / (float)RAND_MAX)-0.5f)
		for (int i = 0; i < _nJunctionsX; ++i)
		{
			for (int j = 0; j < _nJunctionsY; ++j)
			{
				_rand1[i][j] = ((rand() / (float)RAND_MAX) - 0.5f);
				_rand2[i][j] = ((rand() / (float)RAND_MAX) - 0.5f);
				sens[i][j] = _fclipf(_speed + (2.0f * _rand2[i][j] * _randAmt2), getInvSampleRate(), 1.0f);
				rsens[i][j] = sens[i][j];
				_meanLength[i][j] = _fclipf(_time + (2 * _rand1[i][j] * _randAmt1), 1, (float)n);

			}

		}

	}



	float _fclipf(const float _in, const float _min, const float _max)
	{
		return fminf(fmaxf(_in, _min), _max);
	}







};

} // namespace PMWave
