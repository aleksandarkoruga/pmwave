// PluginPMWave.cpp
// Aleksandar Koruga (aleksandar.koruga@gmail.com)

#include "SC_PlugIn.hpp"
#include "PMWave.hpp"

static InterfaceTable* ft;

namespace PMWave {

PMWave::PMWave() {

	{

		auto unit = this;

		//set defaults 


		unit->ddiv = 0.5;
		unit->w = 0;
		unit->envch = 0;
		unit->in_x = 0; // where to write x 
		unit->in_y = 0; //  where to write y 
		unit->envd = 0.03;
		unit->_x = 2; //where to read x
		unit->_y = 2; //where to read y
		unit->_xOld = 0;
		unit->_yOld = 1;

		unit->_randAmt1 = 0.01;
		unit->_randAmt2 = 0.03;



		//1 damp     2 time       3 speed      4 cutoff      5 6 rand    7    trig


		// 0			1				2			3			4			 5	6		7				8				9					10
		//arg in = 0.0, damping=0.5, distance=0.5,speed=0.2 , ,cutoff=0.5,pAmt1,pAmt2, probTrig=0.0, nJunctionsX = 3.0, nJunctionsY = 3.0, delaySamples = 100
		//get constructor arguments
		unit->m_junction = 1.0;
		unit->n_samps = BUFLENGTH;// size of the delay in samples, set it once in constructor.
		//delay time
		unit->_time = static_cast<float>(unit->n_samps) * unit->_fclipf(in(2)[0], 0.0f, 1.0f);
		//speed from one delay time to another
		unit->_speed = unit->_fclipf(in(3)[0], unit->getInvSampleRate(), 1.0f);
		unit->_damping = unit->_fclipf(ZIN0(1), 0.0f, 1.0f);
		unit->_randAmt1 = unit->_fclipf(ZIN0(5), 0.0f, 1.0f);
		unit->_randAmt2 = unit->_fclipf(ZIN0(6), 0.0f, 1.0f);
		unit->_cutoff = unit->_fclipf(ZIN0(4), 0.0f, 1.0f);

		//minimum size 2x2
		unit->_nJunctionsX = NXY;//std::abs((int)std::max((int)ZIN0(8),2));
		unit->_nJunctionsY = NXY;//std::abs((int)std::max((int)ZIN0(9),2));


		for (int i = 0; i < NXY; ++i) { unit->_del[i] = 0; }


		for (int i = 0; i < unit->_nJunctionsX; ++i)
		{
			for (int j = 0; j < unit->_nJunctionsY; ++j)
			{
				for (int _z = 0; _z < unit->n_samps; _z++)
				{
					unit->_p0[i][j][_z] = 0;
					unit->_p1[i][j][_z] = 0;
					unit->_p2[i][j][_z] = 0;
					unit->_p3[i][j][_z] = 0;
				}
				unit->_currentDistance[i][j] = static_cast<float>( unit->n_samps) * 0.5f;
				unit->_p[i][j] = 0;
				unit->_fractional[i][j] = 0;
				unit->_ddup[i][j] = 0;
				unit->_dwn[i][j] = 0;





				unit->_rand2[i][j] = ((rand() / (float)RAND_MAX) - 0.5);
				//std::srand(std::time(nullptr));
				unit->_rand1[i][j] = ((rand() / (float)RAND_MAX) - 0.5);

				unit->_meanLength[i][j] = unit->_fclipf((unit->_time + (2 * unit->_rand1[i][j] * unit->_randAmt1)), 1, static_cast<float>(unit->n_samps));

				unit->sens[i][j] = unit->_fclipf(unit->_speed + (2 * unit->_rand2[i][j] * unit->_randAmt2), unit->getInvSampleRate(), 1.0f);
				unit->rsens[i][j] = unit->sens[i][j];
				// std::cout< <_meanLength[i][j]< <"\n";
			}

		}


		unit->_resetRandom();
		unit->_updateRandAmt1(unit->_randAmt1);
		unit->_updateRandAmt2(unit->_randAmt2);
		// Specify which function is used for audio calculation

		unit->_updateTimeMeanLength(unit->_time);
		unit->_updateSpeed(unit->_speed);




	}




    mCalcFunc = make_calc_function<PMWave, &PMWave::next>();
    next(1);
}

void PMWave::next(int nSamples) {
	auto unit = this;
	//float* output1 = ZOUT(0);
	//float* output2 = ZOUT(1);

	// Pointer to the input buffer
	//float* input = ZIN(0);

	// 0			1				2			3			4			 5	6		7			
	//arg in = 0.0, damping=0.5, distance=0.5,speed=0.2  ,cutoff=0.5,pAmt1,pAmt2, probTrig=0.0

	// Obtain the first argument
	//float thres = ZIN0(1);

	// Use LOOP macro to iterate and ZXP to advance pointers
	for(int cnt=0;cnt<nSamples;++cnt)
	{
	
	unit->_setDamping(in(1)[cnt]);
	if (std::abs(in(2)[cnt] * unit->n_samps - unit->_time) > 0.00001)unit->_updateTimeMeanLength(in(2)[cnt]);

	if (unit->_fclipf(in(3)[cnt], 0, 1) != unit->_speed)unit->_updateSpeed(in(3)[cnt]);

	unit->_setCutoff(in(4)[cnt]);
	unit->m_junction = in(7)[cnt];

	if (unit->_fclipf(in(5)[cnt], 0.0f, 1.0f) != unit->_randAmt1)unit->_updateRandAmt1(in(5)[cnt]);
	if (unit->_fclipf(in(6)[cnt], 0.0f, 1.0f) != unit->_randAmt2)unit->_updateRandAmt2(in(6)[cnt]);

	auto res = unit->_render(in(0)[cnt]);
	out(0)[cnt] =  res[0];
	out(1)[cnt] =  res[1];
		// ZXP(output1) = res[0];
		 // ZXP(output2) = res[1];

}

		
	
	
}


std::vector<float> PMWave::_render(const float _in) 
{


	//MAIN LOOP on plate

	for (int ii = 0; ii < _nJunctionsX; ++ii)
	{
		for (int jj = 0; jj < _nJunctionsY; ++jj)
		{



			//s=0 n=1 e= 2 w=3



			//GLISSANDOS BETWEEN JUNCTIONS

			//distance between value to reach and current value
			float dsens = std::abs(_meanLength[ii][jj] - _currentDistance[ii][jj]);

			sens[ii][jj] = dsens < sens[ii][jj] ? dsens : rsens[ii][jj];
			_currentDistance[ii][jj] += _currentDistance[ii][jj] < _meanLength[ii][jj] ? (sens[ii][jj]) : (_currentDistance[ii][jj] > _meanLength[ii][jj] ? (-sens[ii][jj]) : 0);



			//get int + fract distance from target pointer to current position & +1 (for interpolation)

			_dwn[ii][jj] = (int)_currentDistance[ii][jj];
			_ddup[ii][jj] = _dwn[ii][jj] + 1;
			_fractional[ii][jj] = sc_frac(_currentDistance[ii][jj]);
			
			// ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

									//junction update

			_del[0] = ((_p0[ii][((jj - 1) + _nJunctionsY) &(_nJunctionsY-1)][(w + _dwn[ii][jj]) &(n_samps-1)]) * (1.0 - _fractional[ii][jj]) + (_p0[ii][((jj - 1) + _nJunctionsY) &(_nJunctionsY-1)][(w + _ddup[ii][jj]) &(n_samps-1)]) * _fractional[ii][jj]);
			_del[1] = ((_p1[ii][(jj + 1) &(_nJunctionsY-1)][(w + _dwn[(ii + 1) &(_nJunctionsX-1)][jj]) &(n_samps-1)]) * (1 - _fractional[(ii + 1) &(_nJunctionsX-1)][jj]) + (_p1[ii][(jj + 1) &(_nJunctionsY-1)][(w + _ddup[(ii + 1) &(_nJunctionsX-1)][jj]) &(n_samps-1)]) * _fractional[(ii + 1) &(_nJunctionsX-1)][jj]);
			_del[2] = ((_p2[(ii + 1) &(_nJunctionsX-1)][jj][(w + _dwn[ii][(jj + 1) &(_nJunctionsY-1)]) &(n_samps-1)]) * (1.0 - _fractional[ii][(jj + 1) &(_nJunctionsY-1)]) + (_p2[(ii + 1) &(_nJunctionsX-1)][jj][(w + _ddup[ii][(jj + 1) &(_nJunctionsY-1)]) &(n_samps-1)]) * _fractional[ii][(jj + 1) &(_nJunctionsY-1)]);
			_del[3] = ((_p3[((ii + _nJunctionsX) - 1) &(_nJunctionsX-1)][jj][(w + _dwn[ii][jj]) &(n_samps-1)]) * (1.0 - _fractional[ii][jj]) + (_p3[((ii + _nJunctionsX) - 1) &(_nJunctionsX-1)][jj][(w + _ddup[ii][jj]) &(n_samps-1)]) * _fractional[ii][jj]);

			_p[ii][jj] = _del[0] + _del[1] + _del[2] + _del[3];





			ddiv = 0.5;
			//if we are at the junction selected for the input, add input and set ddiv 0.4
			if ((ii == (int)in_x && jj == (int)in_y))
			{
				_p[ii][jj] *= m_junction;
				_p[ii][jj] += _in;

				ddiv = 0.4;
			}

			//junction result											
			_p[ii][jj] *= ddiv;   //normalize  -  attenuate junction


// ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

						// filtering and ring memory 
						//s=0 n=1 e= 2 w=3
			int _wn = w + n_samps;
			_p0[ii][jj][(_wn - 1) &(n_samps-1)] = (((_damping * _p[ii][jj] - _del[1])) * (1.0f - _cutoff) + _cutoff * (_p0[ii][jj][(_wn - 2) &(n_samps-1)]));
			_p1[ii][jj][(_wn - 1) &(n_samps-1)] = (((_damping * _p[ii][jj] - _del[0])) * (1.0f - _cutoff) + _cutoff * (_p1[ii][jj][(_wn - 2) &(n_samps-1)]));
			_p2[ii][jj][(_wn - 1) &(n_samps-1)] = (((_damping * _p[ii][jj] - _del[3])) * (1.0f - _cutoff) + _cutoff * (_p2[ii][jj][(_wn - 2) &(n_samps-1)]));
			_p3[ii][jj][(_wn - 1) &(n_samps-1)] = (((_damping * _p[ii][jj] - _del[2])) * (1.0f - _cutoff) + _cutoff * (_p3[ii][jj][(_wn - 2) &(n_samps-1)]));


		}
	}

	//OUTPUT SAMPLES and control signals


	// xy are the coordinates given for reading from junctions, xold and yold are for passing with an envelope from one to another.
	//TODO: implement stereo output for supercollider, till then return the selected junction for output
	//float _out1 = (_p[_x][_y] * envch + (1.0f - envch) * _p[_xOld][_yOld]);
	//float _out2 = (_p[_y][_x] * envch + (1.0f - envch) * _p[_yOld][_xOld]);





	w++;
	w = w &( n_samps-1);

	envch = envch + envd;
	envch = envch >= 1 ? 1 : envch;
	// result
	//return({ _p[_x][_y] });
	return { _p[_x][_y],_p[_y][_x] };
}


} // namespace PMWave


PluginLoad(PMWaveUGens) {
	// Plugin magic
	ft = inTable;
	registerUnit<PMWave::PMWave>(ft, "PMWave", false);
}