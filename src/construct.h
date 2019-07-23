#ifndef CONSTRUCT_H
#define CONSTRUCT_H


#include "segment.h"
#include "piping.h"
#include "input.h"
#include <vector>

#include "stru3_matrix.h"



class Construct
{
public:
	Construct(int type, std::vector<SegmentData> segmentData_vec,
			PipingData pipingData, FluidData fluidData);
	void set_variables(double _Q, double _T_in_0, stru3::DVec _T_s);
	void set_functions();
	void calculate_temperatures();
	stru3::DMat assemble_matrix();
	stru3::DVec assemble_RHS();
	
private:
	double L_tot;	// total length
	int N_seg;
	int N_tot; // number of total nodes
	Piping piping;

	std::vector<Segment> segment_vec;

	double T_in_0;

	stru3::DVec T_in;
	stru3::DVec T_out;
	stru3::DVec T_fin;
	stru3::DVec T_fout;
	stru3::DVec T_s;
};

#endif
