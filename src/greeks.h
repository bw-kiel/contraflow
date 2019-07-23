#ifndef GREEKS_H
#define GREEKS_H

#include "piping.h"
#include "utilities.h"

class Greeks
{
public:
	Greeks() {}
	Greeks(double _beta_1, double _beta_12, double _beta, double _gamma, double _delta) :
		beta_1(_beta_1), beta_12(_beta_12), beta(_beta), gamma(_gamma), delta(_delta)
	{
		//LOG("beta_1:  " << beta_1);
		//LOG("beta_12: " << beta_12);
		//LOG("gamma:   " << gamma);
		//LOG("delta:   " << delta);
	}

	double get_beta_1() { return beta_1; }
	double get_beta_12() { return beta_12; }
	double get_beta() { return beta; }
	double get_gamma() { return gamma; }
	double get_delta() { return delta; }

private:
	double beta_1;
	double beta_12;
	double beta;
	double gamma;
	double delta;
};

#endif
