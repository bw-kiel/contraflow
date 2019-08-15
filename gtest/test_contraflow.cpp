#include "contraflow.h"



///////////////////////////////////////////////////////////////////////

TEST(U, SEGMENTS_1)
{
	int N_seg = 1;
	int N = 100;

	std::vector<contra::SegmentData> segmentDataVec;
	for(int i=0; i< N_seg; ++i)
		segmentDataVec.push_back(
				{
					N,		// N
					100., 	// L
					0.13, 	// D
					2.3 	// lambda_g
				});

//testing::internal::CaptureStdout();
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

//EXPECT_EQ("U\n", testing::internal::GetCapturedStdout());	


	stru3::DVec T_s = stru3::DVec(N*N_seg+1);  // number of nodes
	T_s = 10;
	//for(int i=0; i<21;++i)
	//	T_s[i] = i;

	contraflow.calculate(
			2.53e-4, // flow rate Q
			0, 	 //mode, T_in given
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
EXPECT_EQ("80 79.6266 79.2551 78.8857 78.5183 78.1529 77.7895 77.428 77.0685 76.7109 76.3553 76.0016 75.6498 75.3 74.952 74.6059 74.2617 73.9194 73.5789 73.2402 72.9034 72.5685 72.2353 71.904 71.5744 71.2467 70.9207 70.5965 70.274 69.9533 69.6344 69.3172 69.0017 68.6879 68.3758 68.0654 67.7568 67.4497 67.1444 66.8407 66.5387 66.2383 65.9395 65.6424 65.3469 65.053 64.7607 64.47 64.1808 63.8933 63.6073 63.3228 63.0399 62.7586 62.4788 62.2005 61.9237 61.6485 61.3747 61.1024 60.8317 60.5624 60.2945 60.0281 59.7632 59.4998 59.2377 58.9771 58.7179 58.4602 58.2038 57.9488 57.6953 57.4431 57.1923 56.9429 56.6948 56.4481 56.2027 55.9587 55.716 55.4747 55.2346 54.9959 54.7585 54.5224 54.2876 54.054 53.8218 53.5908 53.3611 53.1327 52.9055 52.6795 52.4548 52.2313 52.0091 51.788 51.5682 51.3496 51.1322 ", testing::internal::GetCapturedStdout());	


testing::internal::CaptureStdout();
	for(int i=0; i<N*N_seg+1; ++i)
		std::cout << contraflow.get_result().T_out[i] << " " ;
EXPECT_EQ("34.5081 34.6321 34.7569 34.8824 35.0086 35.1355 35.2632 35.3916 35.5207 35.6506 35.7812 35.9125 36.0447 36.1775 36.3112 36.4456 36.5808 36.7167 36.8534 36.9909 37.1292 37.2683 37.4081 37.5488 37.6903 37.8325 37.9756 38.1195 38.2642 38.4097 38.5561 38.7032 38.8513 39.0001 39.1498 39.3004 39.4517 39.604 39.7571 39.9111 40.0659 40.2216 40.3782 40.5357 40.6941 40.8533 41.0135 41.1745 41.3365 41.4994 41.6632 41.8279 41.9935 42.1601 42.3276 42.496 42.6654 42.8357 43.007 43.1792 43.3524 43.5266 43.7018 43.8779 44.055 44.2331 44.4122 44.5923 44.7734 44.9555 45.1386 45.3228 45.5079 45.6942 45.8814 46.0697 46.259 46.4494 46.6409 46.8334 47.0269 47.2216 47.4173 47.6142 47.8121 48.0111 48.2112 48.4125 48.6148 48.8183 49.0229 49.2286 49.4355 49.6435 49.8527 50.063 50.2745 50.4872 50.701 50.916 51.1322 ", testing::internal::GetCapturedStdout());	



testing::internal::CaptureStdout();
	for(int i=0; i<N_seg; ++i)
	{
		LOG(contraflow.get_result().resistances_vec[i].R_1_Delta << " " << 
				contraflow.get_result().resistances_vec[i].R_1_Delta);
	}
EXPECT_EQ("0.181515 0.181515\n", testing::internal::GetCapturedStdout());	
}

///////////////////////////////////////////////////////////////////////

TEST(U, SEGMENTS_4)
{
	int N_seg = 4;
	int N = 100 / N_seg;

	std::vector<contra::SegmentData> segmentDataVec;
	for(int i=0; i< N_seg; ++i)
		segmentDataVec.push_back(
				{
					N,		// N
					100./N_seg,	// L
					0.13, 	// D
					2.3 	// lambda_g
				});

//testing::internal::CaptureStdout();
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

//EXPECT_EQ("U\n", testing::internal::GetCapturedStdout());	


	stru3::DVec T_s = stru3::DVec(N*N_seg+1);  // number of nodes
	T_s = 10;
	//for(int i=0; i<21;++i)
	//	T_s[i] = i;

	contraflow.calculate(
			2.53e-4, // flow rate Q
			0, 	 //mode, T_in given
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
EXPECT_EQ("80 79.6266 79.2551 78.8857 78.5183 78.1529 77.7895 77.428 77.0685 76.7109 76.3553 76.0016 75.6498 75.3 74.952 74.6059 74.2617 73.9194 73.5789 73.2402 72.9034 72.5685 72.2353 71.904 71.5744 71.2467 70.9207 70.5965 70.274 69.9533 69.6344 69.3172 69.0017 68.6879 68.3758 68.0654 67.7568 67.4497 67.1444 66.8407 66.5387 66.2383 65.9395 65.6424 65.3469 65.053 64.7607 64.47 64.1808 63.8933 63.6073 63.3228 63.0399 62.7586 62.4788 62.2005 61.9237 61.6485 61.3747 61.1024 60.8317 60.5624 60.2945 60.0281 59.7632 59.4998 59.2377 58.9771 58.7179 58.4602 58.2038 57.9488 57.6953 57.4431 57.1923 56.9429 56.6948 56.4481 56.2027 55.9587 55.716 55.4747 55.2346 54.9959 54.7585 54.5224 54.2876 54.054 53.8218 53.5908 53.3611 53.1327 52.9055 52.6795 52.4548 52.2313 52.0091 51.788 51.5682 51.3496 51.1322 ", testing::internal::GetCapturedStdout());	


testing::internal::CaptureStdout();
	for(int i=0; i<N*N_seg+1; ++i)
		std::cout << contraflow.get_result().T_out[i] << " " ;
EXPECT_EQ("34.5081 34.6321 34.7569 34.8824 35.0086 35.1355 35.2632 35.3916 35.5207 35.6506 35.7812 35.9125 36.0447 36.1775 36.3112 36.4456 36.5808 36.7167 36.8534 36.9909 37.1292 37.2683 37.4081 37.5488 37.6903 37.8325 37.9756 38.1195 38.2642 38.4097 38.5561 38.7032 38.8513 39.0001 39.1498 39.3004 39.4517 39.604 39.7571 39.9111 40.0659 40.2216 40.3782 40.5357 40.6941 40.8533 41.0135 41.1745 41.3365 41.4994 41.6632 41.8279 41.9935 42.1601 42.3276 42.496 42.6654 42.8357 43.007 43.1792 43.3524 43.5266 43.7018 43.8779 44.055 44.2331 44.4122 44.5923 44.7734 44.9555 45.1386 45.3228 45.5079 45.6942 45.8814 46.0697 46.259 46.4494 46.6409 46.8334 47.0269 47.2216 47.4173 47.6142 47.8121 48.0111 48.2112 48.4125 48.6148 48.8183 49.0229 49.2286 49.4355 49.6435 49.8527 50.063 50.2745 50.4872 50.701 50.916 51.1322 ", testing::internal::GetCapturedStdout());	




testing::internal::CaptureStdout();
	for(int i=0; i<N_seg; ++i)
	{
		LOG(contraflow.get_result().resistances_vec[i].R_1_Delta << " " << 
				contraflow.get_result().resistances_vec[i].R_1_Delta);
	}
EXPECT_EQ("0.181515 0.181515\n0.181515 0.181515\n0.181515 0.181515\n0.181515 0.181515\n", testing::internal::GetCapturedStdout());	


}

///////////////////////////////////////////////////////////////////////

TEST(_2U, SEGMENTS_1)
{
	int N_seg = 1;
	int N = 100;

	std::vector<contra::SegmentData> segmentDataVec;
	for(int i=0; i< N_seg; ++i)
		segmentDataVec.push_back(
				{
					N,		// N
					100., 	// L
					0.13, 	// D
					2.3 	// lambda_g
				});

//testing::internal::CaptureStdout();
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

//EXPECT_EQ("2U\n", testing::internal::GetCapturedStdout());	


	stru3::DVec T_s = stru3::DVec(N*N_seg+1);  // number of nodes
	T_s =10;
	//for(int i=0; i<21;++i)
	//	T_s[i] = i;

	contraflow.calculate(
			2.53e-4, // flow rate Q
			0, 	 //mode, T_in given
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


//
//	contraflow.calculate(
//			2.53e-4, // flow rate Q
//			1, 	 //mode, DT given
//			10., 	 // temperature difference
//			T_s		 // soil temperature
//	);
//
//	LOG("T_in");
//	for(int i=0; i<N*N_seg+1; ++i)
//		std::cout << contraflow.get_result().T_in[i] << " ";
//	LOG("");
//	LOG("T_out");
//	for(int i=0; i<N*N_seg+1; ++i)
//		std::cout << contraflow.get_result().T_out[i] << " ";
//	LOG("");
//


testing::internal::CaptureStdout();
	for(int i=0; i<N*N_seg+1; ++i)
		std::cout << contraflow.get_result().T_in[i] << " ";
EXPECT_EQ("80 79.2587 78.5256 77.8006 77.0837 76.3747 75.6736 74.9802 74.2945 73.6165 72.946 72.2829 71.6273 70.9789 70.3377 69.7037 69.0768 68.4569 67.8439 67.2377 66.6383 66.0457 65.4596 64.8802 64.3072 63.7407 63.1805 62.6267 62.079 61.5376 61.0022 60.4729 59.9496 59.4321 58.9206 58.4148 57.9148 57.4204 56.9317 56.4485 55.9708 55.4986 55.0317 54.5702 54.114 53.663 53.2171 52.7764 52.3408 51.9101 51.4845 51.0637 50.6478 50.2368 49.8304 49.4289 49.0319 48.6396 48.2519 47.8687 47.49 47.1158 46.7459 46.3804 46.0192 45.6622 45.3095 44.961 44.6166 44.2763 43.9401 43.6079 43.2796 42.9553 42.635 42.3184 42.0058 41.6969 41.3917 41.0903 40.7925 40.4985 40.208 39.9211 39.6377 39.3579 39.0815 38.8086 38.539 38.2729 38.0101 37.7507 37.4945 37.2416 36.9919 36.7454 36.502 36.2618 36.0247 35.7908 35.5598 ", testing::internal::GetCapturedStdout());	


testing::internal::CaptureStdout();
	for(int i=0; i<N*N_seg+1; ++i)
		std::cout << contraflow.get_result().T_out[i] << " " ;
EXPECT_EQ("24.5045 24.5184 24.5341 24.5515 24.5706 24.5914 24.614 24.6383 24.6643 24.6921 24.7216 24.7528 24.7858 24.8206 24.8571 24.8954 24.9354 24.9772 25.0208 25.0661 25.1133 25.1622 25.2129 25.2654 25.3198 25.3759 25.4339 25.4937 25.5553 25.6188 25.6841 25.7513 25.8204 25.8913 25.9641 26.0388 26.1154 26.1939 26.2743 26.3567 26.441 26.5272 26.6154 26.7056 26.7977 26.8918 26.988 27.0861 27.1863 27.2885 27.3927 27.499 27.6074 27.7179 27.8305 27.9451 28.0619 28.1809 28.302 28.4253 28.5507 28.6784 28.8082 28.9403 29.0746 29.2112 29.3501 29.4913 29.6347 29.7805 29.9287 30.0792 30.232 30.3873 30.545 30.7051 30.8677 31.0327 31.2003 31.3703 31.5429 31.718 31.8957 32.076 32.2589 32.4445 32.6327 32.8236 33.0172 33.2135 33.4125 33.6144 33.819 34.0265 34.2368 34.45 34.666 34.885 35.107 35.3319 35.5598 ", testing::internal::GetCapturedStdout());	


testing::internal::CaptureStdout();
	for(int i=0; i<N_seg; ++i)
	{
		LOG(contraflow.get_result().resistances_vec[i].R_1_Delta << " " << 
				contraflow.get_result().resistances_vec[i].R_1_Delta);
	}
EXPECT_EQ("0.106591 0.106591\n", testing::internal::GetCapturedStdout());	
}


///////////////////////////////////////////////////////////////////////
/*

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
			0, 	 //mode, T_in given
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

*/
