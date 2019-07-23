#ifndef SEGMENT_H
#define SEGMENT_H

#include <stdexcept>
#include <vector>
#include "casing.h"
#include "greeks.h"
#include "piping.h"
#include "stru3_matrix.h"
#include "input.h"




class Segment{
public:
	Segment() {}
	Segment(SegmentData segmentData) : casing(segmentData)
	{
		f1 = stru3::DVec(casing.get_N());
		f2 = stru3::DVec(casing.get_N());
		f3 = stru3::DVec(casing.get_N());
	}
	Casing get_casing() { return casing; }
	void set_resistances(Configuration* configuration);
	void set_functions(Piping* piping);

	void set_temperatures(double _T_in_0, double _T_out_0) // on sceleton
	{ T_in_0 = _T_in_0; T_out_0 = _T_out_0; }

	void calculate_temperatures();
	void configure(double* _T_in, double* _T_out, double* _T_s)
	{
		T_in = _T_in; // to second element
		T_out = _T_out; // to second element
		T_s = _T_s; // to first element !!!!!
	}
	void log()
	{
		LOG(T_out[0]);
	}

	double F4(const double &z, const double &a, const double &b);
	double F5(const double &z, const double &a, const double &b);

	double get_f1(int n) { return f1[n]; }
	double get_f2(int n) { return f2[n]; }
	double get_f3(int n) { return f3[n]; }
private:
	Casing casing;
	Greeks greeks;
	double T_in_0, T_out_0;
	double* T_in;
	double* T_out;
	double* T_s;
	
	stru3::DVec f1, f2, f3;	// array

};

#endif
