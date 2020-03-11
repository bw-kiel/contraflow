#include "contraflow.h"


TEST(U, Multipole)
{
	int N_seg = 1;
	int N = 100;

	std::vector<contra::SegmentData> segmentDataVec;
	for(int i=0; i< N_seg; ++i)
		segmentDataVec.push_back(
				contra::SegmentData({
					N,		// N
					100., 	// L
					0.11, 	// D
					3. 	// lambda_g
				}));

//testing::internal::CaptureStdout();
	contra::Contraflow contraflow(0, segmentDataVec,
			{ // piping
					.036,  // d_0_i
					.04,	// d_0_o
					-1.,	// d_1_i
					-1.,	// d_1_o
					.04,	// w
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

//testing::internal::CaptureStdout();
//	for(int i=0; i<N*N_seg+1; ++i)
//		std::cout << contraflow.get_result().T_in[i] << " ";
//EXPECT_EQ("", testing::internal::GetCapturedStdout());	


//testing::internal::CaptureStdout();
//	for(int i=0; i<N*N_seg+1; ++i)
//		std::cout << contraflow.get_result().T_out[i] << " " ;
//EXPECT_EQ("", testing::internal::GetCapturedStdout());	


/*
testing::internal::CaptureStdout();
	for(int i=0; i<N_seg; ++i)
	{
		LOG(contraflow.get_result().resistances_vec[i].R_1_Delta << " " << 
				contraflow.get_result().resistances_vec[i].R_1_Delta);
	}
EXPECT_EQ("0.181515 0.181515\n", testing::internal::GetCapturedStdout());	
*/

}

///////////////////////////////////////////////////////////////////////






TEST(U, SEGMENTS_1)
{
	int N_seg = 1;
	int N = 100;

	std::vector<contra::SegmentData> segmentDataVec;
	for(int i=0; i< N_seg; ++i)
		segmentDataVec.push_back(
				contra::SegmentData({
					N,		// N
					100., 	// L
					0.13, 	// D
					2.3 	// lambda_g
				}));

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

TEST(U, SEGMENTS_1_HighResolution)
{
	int N_seg = 1;
	int N = 200;

	std::vector<contra::SegmentData> segmentDataVec;
	for(int i=0; i< N_seg; ++i)
		segmentDataVec.push_back(
				contra::SegmentData({
					N,		// N
					100., 	// L
					0.13, 	// D
					2.3 	// lambda_g
				}));

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
EXPECT_EQ("80 79.813 79.6266 79.4406 79.2551 79.0702 78.8857 78.7018 78.5183 78.3354 78.1529 77.9709 77.7895 77.6085 77.428 77.248 77.0685 76.8895 76.7109 76.5329 76.3553 76.1782 76.0016 75.8255 75.6498 75.4747 75.3 75.1257 74.952 74.7787 74.6059 74.4336 74.2617 74.0903 73.9194 73.7489 73.5789 73.4093 73.2402 73.0716 72.9034 72.7357 72.5685 72.4017 72.2353 72.0694 71.904 71.739 71.5744 71.4103 71.2467 71.0835 70.9207 70.7584 70.5965 70.435 70.274 70.1135 69.9533 69.7937 69.6344 69.4756 69.3172 69.1592 69.0017 68.8446 68.6879 68.5316 68.3758 68.2204 68.0654 67.9109 67.7568 67.603 67.4497 67.2969 67.1444 66.9923 66.8407 66.6895 66.5387 66.3883 66.2383 66.0887 65.9395 65.7908 65.6424 65.4944 65.3469 65.1997 65.053 64.9066 64.7607 64.6151 64.47 64.3252 64.1808 64.0368 63.8933 63.7501 63.6073 63.4649 63.3228 63.1812 63.0399 62.8991 62.7586 62.6185 62.4788 62.3395 62.2005 62.0619 61.9237 61.7859 61.6485 61.5114 61.3747 61.2384 61.1024 60.9669 60.8317 60.6968 60.5624 60.4283 60.2945 60.1612 60.0281 59.8955 59.7632 59.6313 59.4998 59.3686 59.2377 59.1072 58.9771 58.8473 58.7179 58.5889 58.4602 58.3318 58.2038 58.0761 57.9488 57.8219 57.6953 57.569 57.4431 57.3175 57.1923 57.0674 56.9429 56.8187 56.6948 56.5713 56.4481 56.3252 56.2027 56.0805 55.9587 55.8372 55.716 55.5952 55.4747 55.3545 55.2346 55.1151 54.9959 54.877 54.7585 54.6403 54.5224 54.4048 54.2876 54.1707 54.054 53.9378 53.8218 53.7061 53.5908 53.4758 53.3611 53.2467 53.1327 53.0189 52.9055 52.7923 52.6795 52.567 52.4548 52.3429 52.2313 52.12 52.0091 51.8984 51.788 51.678 51.5682 51.4588 51.3496 51.2408 51.1322 ", testing::internal::GetCapturedStdout());	


testing::internal::CaptureStdout();
	for(int i=0; i<N*N_seg+1; ++i)
		std::cout << contraflow.get_result().T_out[i] << " " ;
EXPECT_EQ("34.5081 34.57 34.6321 34.6944 34.7569 34.8195 34.8824 34.9454 35.0086 35.0719 35.1355 35.1992 35.2632 35.3273 35.3916 35.456 35.5207 35.5855 35.6506 35.7158 35.7812 35.8468 35.9125 35.9785 36.0447 36.111 36.1775 36.2443 36.3112 36.3783 36.4456 36.5131 36.5808 36.6486 36.7167 36.785 36.8534 36.9221 36.9909 37.06 37.1292 37.1986 37.2683 37.3381 37.4081 37.4784 37.5488 37.6194 37.6903 37.7613 37.8325 37.904 37.9756 38.0474 38.1195 38.1917 38.2642 38.3368 38.4097 38.4828 38.5561 38.6296 38.7032 38.7771 38.8513 38.9256 39.0001 39.0749 39.1498 39.225 39.3004 39.3759 39.4517 39.5278 39.604 39.6804 39.7571 39.834 39.9111 39.9884 40.0659 40.1437 40.2216 40.2998 40.3782 40.4569 40.5357 40.6148 40.6941 40.7736 40.8533 40.9333 41.0135 41.0939 41.1745 41.2554 41.3365 41.4178 41.4994 41.5812 41.6632 41.7454 41.8279 41.9106 41.9935 42.0767 42.1601 42.2437 42.3276 42.4117 42.496 42.5806 42.6654 42.7504 42.8357 42.9212 43.007 43.093 43.1792 43.2657 43.3524 43.4394 43.5266 43.6141 43.7018 43.7897 43.8779 43.9663 44.055 44.1439 44.2331 44.3225 44.4122 44.5021 44.5923 44.6827 44.7734 44.8643 44.9555 45.0469 45.1386 45.2306 45.3228 45.4152 45.5079 45.6009 45.6942 45.7876 45.8814 45.9754 46.0697 46.1642 46.259 46.3541 46.4494 46.545 46.6409 46.737 46.8334 46.93 47.0269 47.1241 47.2216 47.3193 47.4173 47.5156 47.6142 47.713 47.8121 47.9115 48.0111 48.111 48.2112 48.3117 48.4125 48.5135 48.6148 48.7164 48.8183 48.9204 49.0229 49.1256 49.2286 49.3319 49.4355 49.5393 49.6435 49.7479 49.8527 49.9577 50.063 50.1686 50.2745 50.3807 50.4872 50.5939 50.701 50.8084 50.916 51.024 51.1322 ", testing::internal::GetCapturedStdout());	







testing::internal::CaptureStdout();
	for(int i=0; i<N_seg; ++i)
	{
		LOG(contraflow.get_result().resistances_vec[i].R_1_Delta << " " << 
				contraflow.get_result().resistances_vec[i].R_1_Delta);
	}
EXPECT_EQ("0.181515 0.181515\n", testing::internal::GetCapturedStdout());	


}

///////////////////////////////////////////////////////////////////////

TEST(_2U, SEGMENTS_1)
{
	int N_seg = 1;
	int N = 100;

	std::vector<contra::SegmentData> segmentDataVec;
	for(int i=0; i< N_seg; ++i)
		segmentDataVec.push_back(
				contra::SegmentData({
					N,		// N
					100., 	// L
					0.13, 	// D
					2.3 	// lambda_g
				}));

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
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
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
