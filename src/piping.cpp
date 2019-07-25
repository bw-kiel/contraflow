#include "piping.h"
#include "configuration.h"
#include "interface.h"

namespace contra
{

Piping::Piping(PipingData pipingData, FluidData fluidData) : fluid(fluidData),
		d_0_i(pipingData.d_0_i), d_0_o(pipingData.d_0_o),
		d_1_i(pipingData.d_1_i), d_1_o(pipingData.d_1_o),
		w(pipingData.w),
		lambda_0(pipingData.lambda_0), lambda_1(pipingData.lambda_1),
		d(0.)
{}

void Piping::configure(int type)
{
	if(type == 0)
		configuration = new Configuration_U(this);
	else if(type == 1)
		configuration = new Configuration_2U(this);
	else if(type == 2)
		configuration = new Configuration_CX(this);
	else
		throw std::runtime_error("Segment factory: type unknown"); 
}

void Piping::set_flow(double L)
{	
	configuration->set_flow(L);
}

Configuration* Piping::get_configuration() { return configuration; }

}
