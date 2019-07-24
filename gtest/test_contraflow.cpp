#include "contraflow.h"



TEST(U, SEGMENTS_4)
{
	std::vector<SegmentData> segmentDataVec;
	for(int i=0; i< 4; ++i)
		segmentDataVec.push_back(
				{
					5,		// N
					100., 	// L
					0.13, 	// D
					2.3 	// lambda_g
				});

testing::internal::CaptureStdout();
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

EXPECT_EQ("U\n", testing::internal::GetCapturedStdout());	


	stru3::DVec T_s = stru3::DVec(21);  // number of nodes
	T_s =10;
	//for(int i=0; i<21;++i)
	//	T_s[i] = i;

	contraflow.calculate(
			2.53e-4, // flow rate Q
			80, 	 // feed flow temperature T_in_0
			T_s		 // soil temperature
	);

testing::internal::CaptureStdout();
	for(int i=0; i<21; ++i)
		LOG(contraflow.get_result().T_in[i]);
EXPECT_EQ("80\n72.7439\n66.2294\n60.3805\n55.1293\n50.8447\n46.5698\n42.7319\n39.2863\n36.193\n33.8462\n31.3112\n29.0357\n26.9933\n25.1601\n23.945\n22.4266\n21.0644\n19.8425\n18.7467\n18.1944\n", testing::internal::GetCapturedStdout());	


testing::internal::CaptureStdout();
	for(int i=0; i<21; ++i)
		LOG(contraflow.get_result().T_out[i]);
EXPECT_EQ("11.7221\n11.6195\n11.5228\n11.4308\n11.3426\n12.0657\n12.0728\n12.091\n12.1206\n12.1619\n13.0241\n13.1815\n13.3628\n13.5704\n13.8064\n14.8826\n15.2755\n15.7168\n16.2116\n16.7657\n18.1944\n", testing::internal::GetCapturedStdout());	



testing::internal::CaptureStdout();
	for(int i=0; i<4; ++i)
	{
		LOG(contraflow.get_result().resistances_vec[i].R_1_Delta << " " << 
				contraflow.get_result().resistances_vec[i].R_1_Delta);
	}
EXPECT_EQ("0.181525 0.181525\n0.181525 0.181525\n0.181525 0.181525\n0.181525 0.181525\n", testing::internal::GetCapturedStdout());	

}

