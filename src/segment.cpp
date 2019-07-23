#include "segment.h"
#include "configuration.h"
#include "utilities.h"
#include <math.h>


void Segment::set_resistances(Configuration* configuration)
{
	configuration->set_resistances(casing.get_D(), casing.get_lambda_g());
}


void Segment::set_functions(Piping* piping)
{
	greeks = piping->get_configuration()->set_greeks(piping);

	const double dgz = greeks.get_gamma() * casing.get_L() / casing.get_N();
	double gz = dgz;

	for(int i=0; i<casing.get_N(); ++i)
	{
		const double ch = cosh(gz);
		const double sh = sinh(gz);
		const double ebz = exp(greeks.get_beta() * gz / greeks.get_gamma()); // 1. if not CX

		f1[i] = ebz * (ch - sh * greeks.get_delta());
		f2[i] = ebz * (sh * greeks.get_beta_12() / greeks.get_gamma()); 
		f3[i] = ebz * (ch + sh * greeks.get_delta());

		//LOG("z:  " << gz/greeks.get_gamma());
		//LOG("	f1: " << f1[i]);
		//LOG("	f2: " << f2[i]);
		//LOG("	f3: " << f3[i]);

		gz += dgz;
	} 
}

double Segment::F4(const double &z, const double &a, const double &b)
{
// not valid for CX
	double beta_1 = greeks.get_beta_1();
	double beta_12 = greeks.get_beta_12();
	double gamma = greeks.get_gamma();
	double delta = greeks.get_delta();

	double res = (sinh(gamma*(z-a)) - sinh(gamma*(z-b))) * beta_1 / gamma;
	res += (cosh(gamma*(z-b)) - cosh(gamma*(z-a))) * 
				(beta_1 / gamma) * (delta + beta_12 / gamma);
	return res;
}

double Segment::F5(const double &z, const double &a, const double &b)
{
// not valid for CX
	double beta_1 = greeks.get_beta_1();
	double beta_12 = greeks.get_beta_12();
	double gamma = greeks.get_gamma();
	double delta = greeks.get_delta();

	double res = (sinh(gamma*(z-a)) - sinh(gamma*(z-b))) * beta_1 / gamma;
	res -= (cosh(gamma*(z-b)) - cosh(gamma*(z-a))) * 
				(beta_1 / gamma) * (delta + beta_12 / gamma);
	return res;
}

void Segment::calculate_temperatures()
{
	const int N = casing.get_N();

	for(int i=0; i<N-1; ++i)  // last node from sceleton
	{
		const double dz = casing.get_L() / N;
		
		T_in[i] = T_in_0 * f1[i] + T_out_0 * f2[i];
		T_out[i] = -T_in_0 * f2[i] + T_out_0 * f3[i];

		for(int j=0; j<i+1; ++j)
		{
			T_in[i] += F4((i+1)*dz, j*dz, (j+1)*dz)*(T_s[i]+T_s[i+1])/2;
			T_out[i] -= F5((i+1)*dz, j*dz, (j+1)*dz)*(T_s[i]+T_s[i+1])/2;
		}
	}
}
