#include "contraflow.h"

TEST(U, SEGMENTS_4)
{
	int N_seg = 4;
	int N = 5;

	std::vector<contra::SegmentData> segmentDataVec;
	for(int i=0; i< N_seg; ++i)
		segmentDataVec.push_back(
				{
					N,		// N
					100., 	// L
					0.13, 	// D
					2.3 	// lambda_g
				});

testing::internal::CaptureStdout();
	contra::Contraflow contraflow(0, segmentDataVec,
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


	stru3::DVec T_s = stru3::DVec(N*N_seg+1);  // number of nodes
	T_s =10;
	//for(int i=0; i<21;++i)
	//	T_s[i] = i;

	contraflow.calculate(
			2.53e-4, // flow rate Q
			80, 	 // feed flow temperature T_in_0
			T_s		 // soil temperature
	);

	LOG("T_in");
	for(int i=0; i<N*N_seg+1; ++i)
		std::cout << contraflow.get_result().T_in[i] << " ";
	LOG("");
	LOG("T_out");
	for(int i=0; i<N*N_seg+1; ++i)
		std::cout << contraflow.get_result().T_out[i] << " ";
	LOG("");

testing::internal::CaptureStdout();
	for(int i=0; i<N*N_seg+1; ++i)
		std::cout << contraflow.get_result().T_in[i] << " ";
EXPECT_EQ("80 72.8487 66.4284 60.6646 55.4901 50.8447 46.6745 42.9309 39.5704 36.5538 33.8462 31.4159 29.2347 27.2773 25.5209 23.945 22.5314 21.2635 20.1266 19.1075 18.1944 ", testing::internal::GetCapturedStdout());	


testing::internal::CaptureStdout();
	for(int i=0; i<N*N_seg+1; ++i)
		std::cout << contraflow.get_result().T_in[i] << " " ;
EXPECT_EQ("80 72.8487 66.4284 60.6646 55.4901 50.8447 46.6745 42.9309 39.5704 36.5538 33.8462 31.4159 29.2347 27.2773 25.5209 23.945 22.5314 21.2635 20.1266 19.1075 18.1944 ", testing::internal::GetCapturedStdout());	



testing::internal::CaptureStdout();
	for(int i=0; i<4; ++i)
	{
		LOG(contraflow.get_result().resistances_vec[i].R_1_Delta << " " << 
				contraflow.get_result().resistances_vec[i].R_1_Delta);
	}
EXPECT_EQ("0.181525 0.181525\n0.181525 0.181525\n0.181525 0.181525\n0.181525 0.181525\n", testing::internal::GetCapturedStdout());	
}

///////////////////////////////////////////////////////////////////////

TEST(_2U, SEGMENTS_4)
{
	int N_seg = 4;
	int N = 5;

	std::vector<contra::SegmentData> segmentDataVec;
	for(int i=0; i< N_seg; ++i)
		segmentDataVec.push_back(
				{
					N,		// N
					100., 	// L
					0.13, 	// D
					2.3 	// lambda_g
				});

testing::internal::CaptureStdout();
	contra::Contraflow contraflow(1, segmentDataVec,
			{ // piping
					.0262,  // d_0_i
					.032,	// d_0_o
					.0262,	// d_1_i
					.032,	// d_1_o
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

EXPECT_EQ("2U\n", testing::internal::GetCapturedStdout());	


	stru3::DVec T_s = stru3::DVec(N*N_seg+1);  // number of nodes
	T_s =10;
	//for(int i=0; i<21;++i)
	//	T_s[i] = i;

	contraflow.calculate(
			2.53e-4, // flow rate Q
			80, 	 // feed flow temperature T_in_0
			T_s		 // soil temperature
	);

	LOG("T_in");
	for(int i=0; i<N*N_seg+1; ++i)
		std::cout << contraflow.get_result().T_in[i] << " ";
	LOG("");
	LOG("T_out");
	for(int i=0; i<N*N_seg+1; ++i)
		std::cout << contraflow.get_result().T_out[i] << " ";
	LOG("");

/*
testing::internal::CaptureStdout();
	for(int i=0; i<N*N_seg+1; ++i)
		std::cout << contraflow.get_result().T_in[i] << " ";
EXPECT_EQ("80 72.8487 66.4284 60.6646 55.4901 50.8447 46.6745 42.9309 39.5704 36.5538 33.8462 31.4159 29.2347 27.2773 25.5209 23.945 22.5314 21.2635 20.1266 19.1075 18.1944 ", testing::internal::GetCapturedStdout());	


testing::internal::CaptureStdout();
	for(int i=0; i<N*N_seg+1; ++i)
		std::cout << contraflow.get_result().T_in[i] << " " ;
EXPECT_EQ("80 72.8487 66.4284 60.6646 55.4901 50.8447 46.6745 42.9309 39.5704 36.5538 33.8462 31.4159 29.2347 27.2773 25.5209 23.945 22.5314 21.2635 20.1266 19.1075 18.1944 ", testing::internal::GetCapturedStdout());	



testing::internal::CaptureStdout();
	for(int i=0; i<4; ++i)
	{
		LOG(contraflow.get_result().resistances_vec[i].R_1_Delta << " " << 
				contraflow.get_result().resistances_vec[i].R_1_Delta);
	}
EXPECT_EQ("0.181525 0.181525\n0.181525 0.181525\n0.181525 0.181525\n0.181525 0.181525\n", testing::internal::GetCapturedStdout());	
*/
}


///////////////////////////////////////////////////////////////////////

TEST(CX, SEGMENTS_4)
{
	int N_seg = 4;
	int N = 5;

	std::vector<contra::SegmentData> segmentDataVec;
	for(int i=0; i< N_seg; ++i)
		segmentDataVec.push_back(
				{
					N,		// N
					100., 	// L
					0.1, 	// D
					2.3 	// lambda_g
				});

testing::internal::CaptureStdout();
	contra::Contraflow contraflow(2, segmentDataVec,
			{ // piping
					.018,  // d_0_i
					.024,	// d_0_o
					.042,	// d_1_i
					.050,	// d_1_o
					.0,	// w
					.38,	// lambda_0
					.38		// lambda_1
			},
			{ // fluid
					.6405,		// lambda
					.54741e-3, 	// mu
					4180.95, 	// c
					988.1 		// rho
			});  // type, number of segments

EXPECT_EQ("CX\n", testing::internal::GetCapturedStdout());	


	stru3::DVec T_s = stru3::DVec(N*N_seg+1);  // number of nodes
	T_s =10;
	//for(int i=0; i<21;++i)
	//	T_s[i] = i;

	contraflow.calculate(
			2.53e-4, // flow rate Q
			80, 	 // feed flow temperature T_in_0
			T_s		 // soil temperature
	);

	LOG("T_in");
	for(int i=0; i<N*N_seg+1; ++i)
		std::cout << contraflow.get_result().T_in[i] << " ";
	LOG("");
	LOG("T_out");
	for(int i=0; i<N*N_seg+1; ++i)
		std::cout << contraflow.get_result().T_out[i] << " ";
	LOG("");
}
