PMWave : MultiOutUGen {
	*ar { |in = 0.0, decay = 4.0, time = 0.3, glide = 0.2, hfDecay = 1.0, rand1 = 0.1, rand2 = 0.2, wander = 0.0, glideMode = 0|
		^this.multiNew('audio', in, decay, time, glide, hfDecay, rand1, rand2, wander, glideMode);
	}
	init { |... theInputs|
		inputs = theInputs;
		^this.initOutputs(2, rate);
	}
}
