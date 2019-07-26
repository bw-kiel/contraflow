#include "segment.h"
#include "configuration.h"
#include "utilities.h"
#include <math.h>

namespace contra
{

Resistances Segment::set_resistances(Configuration* configuration)
{
	return configuration->set_resistances(casing.get_D(), casing.get_lambda_g());
}


void Segment::set_functions(Piping* piping)
{
	Configuration* configuration = piping->get_configuration();
	greeks = configuration->set_greeks(piping);

	const double dgz = greeks.get_gamma() * casing.get_L() / casing.get_N();
	double gz = dgz;

	for(int i=0; i<casing.get_N(); ++i)
	{
		configuration->set_functions(f1[i], f2[i], f3[i], gz, greeks);		
		//LOG("z:  " << gz/greeks.get_gamma());
		//LOG("	f1: " << f1[i]);
		//LOG("	f2: " << f2[i]);
		//LOG("	f3: " << f3[i]);

		gz += dgz;
	} 
}

void Segment::calculate_temperatures(Configuration* configuration)
{
	const int N = casing.get_N();

	for(int i=1; i<N; ++i)  // last node from sceleton
	{
		const double dz = casing.get_L() / N;
		
		T_in[i] = T_in[0] * f1[i-1] + T_out[0] * f2[i-1];
		T_out[i] = -T_in[0] * f2[i-1] + T_out[0] * f3[i-1];

		for(int j=0; j<i; ++j)
		{
			T_in[i] += configuration->F4((i+1)*dz, j*dz, (j+1)*dz, greeks)*(T_s[i]+T_s[i+1])/2;
			T_out[i] -= configuration->F5((i+1)*dz, j*dz, (j+1)*dz, greeks)*(T_s[i]+T_s[i+1])/2;
		}
	}
}

}

