#include "construct.h"
#include "stru3_matrix.h"
#include "stru3_gauss.cpp"



Construct::Construct(int type, std::vector<SegmentData> segmentData_vec,
		PipingData pipingData, FluidData fluidData) :
					L_tot(0.), N_seg(segmentData_vec.size()),
					piping(pipingData, fluidData)
{
	N_tot = 1;

	for(int i=0; i<N_seg;++i)
	{
		segment_vec.push_back(Segment(segmentData_vec[i]));
		L_tot += segment_vec[i].get_casing().get_L();
		N_tot += segment_vec[i].get_casing().get_N();
	}

	T_in = stru3::DVec(N_tot);
	T_out = stru3::DVec(N_tot);

	T_fin = stru3::DVec(N_tot);
	T_fout = stru3::DVec(N_tot);

	T_s = stru3::DVec(N_tot);

	piping.configure(type);

	int j=N_tot-1;
	for(int i=0; i<N_seg;++i)
	{
		j -= segment_vec[i].get_casing().get_N();
		segment_vec[i].configure(&T_in[j+1], &T_out[j+1], &T_s[j]);
	}
}

void Construct::set_variables(double _Q, double _T_in_0, stru3::DVec _T_s)
{
	piping.set_Q(_Q);
	T_in_0 = _T_in_0;
	T_s = _T_s;
}

void Construct::set_functions()
{
	piping.set_flow(L_tot);
	for(int i=0; i<segment_vec.size(); ++i)
	{	
		segment_vec[i].set_resistances(piping.get_configuration());
		segment_vec[i].set_functions(&piping);
	}
}

void Construct::calculate_temperatures()
{
	int dim = 2 * N_seg;
	stru3::DMat m = assemble_matrix();
	stru3::DVec b = assemble_RHS();

	// solve matrix
	stru3::classical_elimination(m, b);
	stru3::DVec T = stru3::back_substitution(m, b);  // temperatures on sceleton

	for(int i=0; i < N_seg-1; ++i)
		segment_vec[i].set_temperatures(T[2+2*i], T[1+2*i]);

	segment_vec[N_seg-1].set_temperatures(T_in_0, T[1+2*(N_seg-1)]);

	T_in[N_tot-1] = T[0];
	T_out[N_tot-1] = T[0];
	for(int i=1; i<N_seg; ++i)
	{
		T_in[N_tot-1-i*5] = T[2*i];
		T_out[N_tot-1-i*5] = T[2*i-1];
	}

	T_in[0] = T_in_0;
	T_out[0] = T[2*N_seg-1];

	//for(int i=0; i < N_tot; ++i)
	//	LOG(T_out[i]);
	//LOG("....");
	//for(int i=0; i < N_seg; ++i)
	//segment_vec[i].log();


	LOG("T_in:");
	for(int i=0; i < N_seg; ++i)
		segment_vec[i].calculate_temperatures();
	
	for(int i=0; i < N_tot; ++i)
		LOG(T_in[i]);

	LOG("T_out:");
	for(int i=0; i < N_tot; ++i)
		LOG(T_out[i]);

	//for(int i=0; i < N_tot; ++i)
	//	LOG(T_s[i]);
;
}

stru3::DMat Construct::assemble_matrix()
{
	int dim = 2 * N_seg;
	int N;
	stru3::DMat m(dim, dim);  // initializes with 0.

	m(0, 0) = -2.;
	for(int i=1; i < N_seg; ++i)
	{
		m(2*i, 2*i-1) = -1.;
		m(2*i, 2*i) = -1.;

		m(2*i+1, 2*i-1) = -1.;
		m(2*i+1, 2*i) = 1.;
	}

	for(int i=0; i < N_seg-1; ++i)
	{
		N = segment_vec[i].get_casing().get_N();

		m(2*i, 2*i+1) = segment_vec[i].get_f2(N-1) + segment_vec[i].get_f3(N-1);
		m(2*i, 2*i+2) = segment_vec[i].get_f1(N-1) - segment_vec[i].get_f2(N-1);

		m(2*i+1, 2*i+1) = segment_vec[i].get_f3(N-1) - segment_vec[i].get_f2(N-1);
		m(2*i+1, 2*i+2) = -(segment_vec[i].get_f1(N-1) + segment_vec[i].get_f2(N-1));
	}

	N = segment_vec[N_seg-1].get_casing().get_N();
	m(dim-2, dim-1) = segment_vec[N_seg-1].get_f2(N-1) + segment_vec[N_seg-1].get_f3(N-1);
	m(dim-1, dim-1) = segment_vec[N_seg-1].get_f3(N-1) - segment_vec[N_seg-1].get_f2(N-1);

	return m;
}

stru3::DVec Construct::assemble_RHS()
{
	int dim = 2 * N_seg;
	int N;
	stru3::DVec b(dim);

	int ii = 0;
	for(int i=0; i < N_seg; ++i)
	{
		N = segment_vec[i].get_casing().get_N();
		double L = segment_vec[i].get_casing().get_L();
		double dz = L / N;	
		double z = 0;
		for(int j=0; j < N; ++j)
		{
			b[2*i] -= (segment_vec[i].F4(L, z, z+dz) - segment_vec[i].F5(L, z, z+dz)) * (T_s[ii] + T_s[ii+1]) / 2;
			b[2*i+1] += (segment_vec[i].F4(L, z, z+dz) + segment_vec[i].F5(L, z, z+dz)) * (T_s[ii] + T_s[ii+1]) / 2;

			z += dz;
			ii++;
		}
	}

	N = segment_vec[N_seg-1].get_casing().get_N();
	b[dim-2] -= (segment_vec[N_seg-1].get_f1(N-1) - segment_vec[N_seg-1].get_f2(N-1)) * T_in_0;
	b[dim-1] += (segment_vec[N_seg-1].get_f1(N-1) + segment_vec[N_seg-1].get_f2(N-1)) * T_in_0;

	return b;
}

