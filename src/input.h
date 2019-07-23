#ifndef INPUT_H_
#define INPUT_H_

struct SegmentData
{
	int N;
	double L;
	double D;
	double lambda_g;
};


struct FluidData
{
	double lambda;
	double mu;
	double c;
	double rho;
};


struct PipingData
{
	double d_0_i;
	double d_0_o;
	double d_1_i;
	double d_1_o;
	double w;
	double lambda_0;
	double lambda_1;
};




#endif
