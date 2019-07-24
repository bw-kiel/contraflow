#include "contraflow.h"

#define N_SEG 4
#define N_ELE 5

int main()
{
	std::vector<SegmentData> segmentDataVec;
	for(int i=0; i< N_SEG; ++i)
		segmentDataVec.push_back(
				{
					N_ELE,		// N
					100., 	// L
					0.13, 	// D
					2.3 	// lambda_g
				});

	Contraflow contraflow(0, segmentDataVec,
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

	contraflow.calculate(
			2.53e-4, // flow rate Q
			80, 	 // feed flow temperature T_in_0
			T_s		 // soil temperature
	);

	LOG("T_in:");
	for(int i=0; i<N_SEG * N_ELE + 1; ++i)
		LOG(contraflow.get_result().T_in[i]);

	LOG("T_out:");
	for(int i=0; i<N_SEG * N_ELE + 1; ++i)
		LOG(contraflow.get_result().T_out[i]);

	
	LOG("Resistances:");
	for(int i=0; i<N_SEG; ++i)
	{
		LOG(contraflow.get_result().resistances_vec[i].R_1_Delta << " " << 
				contraflow.get_result().resistances_vec[i].R_1_Delta);
	}

	return 0;
}
