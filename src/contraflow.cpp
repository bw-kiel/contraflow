#include "contraflow.h"
#include "stru3_matrix.h"
#include "stru3_gauss.cpp"

namespace contra
{


Contraflow::Contraflow(int type, std::vector<SegmentData> segmentData_vec,
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

	result.T_in.resize(N_tot);
	result.T_out.resize(N_tot);

	result.T_fin.resize(N_tot);
	result.T_fout.resize(N_tot);

	T_s.resize(N_tot);

	piping.configure(type);

	int j=N_tot-1;
	for(int i=0; i<N_seg;++i)
	{
		j -= segment_vec[i].get_casing().get_N();
		segment_vec[i].configure(&(result.T_in[j]), &(result.T_out[j]), &T_s[j]);
		result.resistances_vec.push_back(Resistances({0., 0.}));  // initialize
	}
}


void Contraflow::set_functions()
{
	piping.set_flow(L_tot);
	for(int i=0; i<segment_vec.size(); ++i)
	{	
		result.resistances_vec[i] = segment_vec[i].set_resistances(piping.get_configuration());
		segment_vec[i].set_functions(&piping);
	}
}

void Contraflow::calculate(double _Q, int mode, double var, stru3::DVec _T_s)
{
	if(mode == 0)
		T_in_0 = var;

	piping.set_Q(_Q);
	T_s = _T_s;

	set_functions();

	DEBUG("assemble matrix");
	stru3::DMat m = assemble_matrix(mode);
	DEBUG("assemble RHS");
	stru3::DVec b = assemble_RHS(mode, var);

	//LOG(m);
	//for(int i=0; i< 2*N_seg; ++i)
	//	LOG(b[i]);
	//LOG("");
	// solve matrix
	DEBUG("classical elimination");
	stru3::classical_elimination(m, b);
	DEBUG("back substitution");
	stru3::DVec T_scel = stru3::back_substitution(m, b);  // temperatures on sceleton

	result.T_in[N_tot-1] = T_scel[0];
	result.T_out[N_tot-1] = T_scel[0];

	//LOG("result:");
	//for(int i=0; i< 2*N_seg; ++i)
	//	LOG(T_scel[i]);

	int j=N_tot-1;

	if(mode == 0)
	{
		for(int i=0; i<N_seg-1;++i)
		{
			j -= segment_vec[i].get_casing().get_N();
	
			result.T_in[j] = T_scel[2*(i+1)];
			result.T_out[j] = T_scel[2*(i+1)-1];
		}
	
		result.T_in[0] = T_in_0;
		result.T_out[0] = T_scel[2*N_seg-1];
	}
	else
	{
		for(int i=0; i<N_seg;++i)
		{
			j -= segment_vec[i].get_casing().get_N();
	
			result.T_in[j] = T_scel[2*(i+1)];
			result.T_out[j] = T_scel[2*(i+1)-1];
		}

	}

	for(int i=0; i < N_seg; ++i)
	{
		segment_vec[i].calculate_temperatures(piping.get_configuration());
	}

	
}

stru3::DMat Contraflow::assemble_matrix(int mode)
{
	int dim;

	if(mode == 0)
		dim = 2 * N_seg;
	else
		dim = 2 * N_seg + 1;

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

	if(mode == 0)
	{
		for(int i=0; i < N_seg-1; ++i)
		{
			N = segment_vec[i].get_casing().get_N();
	
			m(2*i, 2*i+1) = segment_vec[i].get_f2(N-1) + segment_vec[i].get_f3(N-1);  // a
			m(2*i, 2*i+2) = segment_vec[i].get_f1(N-1) - segment_vec[i].get_f2(N-1);  // b
	
			m(2*i+1, 2*i+1) = segment_vec[i].get_f3(N-1) - segment_vec[i].get_f2(N-1);  // ta
			m(2*i+1, 2*i+2) = -(segment_vec[i].get_f1(N-1) + segment_vec[i].get_f2(N-1));  // tb
		}
	
		N = segment_vec[N_seg-1].get_casing().get_N();
		m(dim-2, dim-1) = segment_vec[N_seg-1].get_f2(N-1) + segment_vec[N_seg-1].get_f3(N-1);  // a
		m(dim-1, dim-1) = segment_vec[N_seg-1].get_f3(N-1) - segment_vec[N_seg-1].get_f2(N-1);  // ta
	}
	else
	{
		for(int i=0; i < N_seg; ++i)
		{
			N = segment_vec[i].get_casing().get_N();
	
			m(2*i, 2*i+1) = segment_vec[i].get_f2(N-1) + segment_vec[i].get_f3(N-1);
			m(2*i, 2*i+2) = segment_vec[i].get_f1(N-1) - segment_vec[i].get_f2(N-1);
	
			m(2*i+1, 2*i+1) = segment_vec[i].get_f3(N-1) - segment_vec[i].get_f2(N-1);
			m(2*i+1, 2*i+2) = -(segment_vec[i].get_f1(N-1) + segment_vec[i].get_f2(N-1));
		}
	
		m(dim-1, dim-2) = 1.;
		m(dim-1, dim-1) = -1.;
	}
	return m;
}

stru3::DVec Contraflow::assemble_RHS(int mode, double var)
{
	int dim;
	if(mode == 0)
		dim = 2 * N_seg;
	else
		dim = 2 * N_seg + 1;
	int N;
	stru3::DVec b(dim);
	Configuration* configuration = piping.get_configuration();
	if(configuration == NULL)
		std::runtime_error("no configuration");

	int ii = 0;
	for(int i=0; i < N_seg; ++i)
	{
		Greeks greeks = segment_vec[i].get_greeks();
		N = segment_vec[i].get_casing().get_N();
		double L = segment_vec[i].get_casing().get_L();
		double dz = L / N;	
		double z = 0;

		for(int j=0; j < N; ++j)
		{
			b[2*i] -= (configuration->F4(L, z, z+dz, greeks) - configuration->F5(L, z, z+dz, greeks)) * (T_s[ii] + T_s[ii+1]) / 2;
			b[2*i+1] += (configuration->F4(L, z, z+dz, greeks) + configuration->F5(L, z, z+dz, greeks)) * (T_s[ii] + T_s[ii+1]) / 2;

			z += dz;
			ii++;
		}
	}

	if(mode == 0)
	{
		N = segment_vec[N_seg-1].get_casing().get_N();
		b[dim-2] -= (segment_vec[N_seg-1].get_f1(N-1) - segment_vec[N_seg-1].get_f2(N-1)) * T_in_0;
		b[dim-1] += (segment_vec[N_seg-1].get_f1(N-1) + segment_vec[N_seg-1].get_f2(N-1)) * T_in_0;
	}
	else
	{
		b[dim-1] = var;  // temperature difference	
	}

	return b;
}

}
