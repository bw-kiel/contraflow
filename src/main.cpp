#include "construct.h"

#define N_SEG 4

int main()
{
	std::vector<SegmentData> segmentDataVec;
	for(int i=0; i< N_SEG; ++i)
		segmentDataVec.push_back(
				{
					5,		// N
					100., 	// L
					0.13, 	// D
					2.3 	// lambda_g
				});

	Construct construct(0, segmentDataVec,
			{ // piping
					.0262,  // d_0_i
					.032,	// d_0_o
					-1.,	// d_1_i
					-1.,	// d_1_o
					.06,	// w
					.38,	// lambda_0
					.38		// lambda_1
			},
			{ // fluid
					.6405,		// lambda
					.54741e-3, 	// mu
					4180.95, 	// c
					988.1 		// rho
			});  // type, number of segments


	stru3::DVec T_s = stru3::DVec(21);  // number of nodes
	T_s =10;
	//for(int i=0; i<21;++i)
	//	T_s[i] = i;

	construct.set_variables(
			2.53e-4, // flow rate Q
			80, 	 // feed flow temperature T_in_0
			T_s		 // soil temperature
	);
	construct.set_functions();
	construct.calculate_temperatures();

	return 0;
}
