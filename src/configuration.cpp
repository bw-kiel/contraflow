#include "configuration.h"
#include "utilities.h"
#include "fluid.h"
#include "greeks.h"
#include <math.h>

namespace contra
{

Configuration_U::Configuration_U(Piping* _piping) : Configuration(_piping)
{
	LOG("U");
	piping->A_0 = M_PI * piping->d_0_i * piping->d_0_i / 4; 
	piping->A_1 = piping->A_0; 

}

void Configuration_U::set_flow(double L)
{
	piping->u_0 = piping->Q / piping->A_0; 
	piping->u_1 = piping->u_0;

	piping->Re_0 = piping->fluid.Reynolds(piping->u_0, piping->d_0_i);
	piping->Re_1 = piping->Re_0;

	piping->Nu_0 = piping->fluid.Nusselt(piping->Re_0, L, piping->d_0_i);
	piping->Nu_1 = piping->Nu_0;
	
	//LOG("u:  " << piping->u_0);
	//LOG("Re: " << piping->Re_0);
	//LOG("Nu: " << piping->Nu_0);
}

Resistances Configuration_U::set_resistances(double D, double lambda_g)
{
	double D2 = D * D;
	double d2 = piping->d_0_o * piping->d_0_o;
	double w2 = piping->w * piping->w;
	double l = 2 * M_PI * lambda_g;

	double x = log(sqrt(D*D + 2 * d2)/(2*piping->d_0_o)) / 
		log(D/(1.41421356237*piping->d_0_o));
	double R_g = (1.601 - 0.888 * piping->w / D) * acosh((D2 + d2 - w2)/(2 * D * piping->d_0_o)) / l;
	double R_ar = acosh((2*w2 - d2)/d2) / l;

	R_adv = 1 / (piping->Nu_0 * piping->fluid.get_lambda() * M_PI);
	R_con_a = log(piping->d_0_o / piping->d_0_i) /(2*piping->lambda_0*M_PI);
	R_con_b = x * R_g;

	R_gs = (1 - x) * R_g;
	R_fg = R_adv + R_con_a + R_con_b;

	double A = 2 * R_gs;
	double B = R_ar - 2 * x * R_g;
	R_gg = A * B / (A - B);

	double u_a = (1/R_fg) + (1/R_gs) + (1/R_gg);
	R_1_Delta = R_fg + R_gs;
	R_2_Delta = R_1_Delta;
	double C = u_a * R_fg * R_gg;
	R_12_Delta = (C * C - R_fg * R_fg) / R_gg;

	/*LOG("x:          " << x);
	LOG("R_g:        " << R_g);
	LOG("R_ar:       " << R_ar);

	LOG("R_adv:      " << R_adv);
	LOG("R_con_a:    " << R_con_a);
	LOG("R_con_b:    " << R_con_b);
	LOG("R_gs:       " << R_gs);
	LOG("R_fg:       " << R_fg);
	LOG("R_gg:       " << R_gg);
	LOG("R_1_Delta:  " << R_1_Delta);
	LOG("R_12_Delta: " << R_12_Delta);*/
	return {R_1_Delta, R_2_Delta};
}

Greeks Configuration_U::set_greeks(Piping* piping)
{
	double factor = piping->get_fluid().get_c_vol() * piping->u_0 * piping->A_0;
	double beta_1 = 1. / (factor * R_1_Delta);
	double beta_12 = 1. / (factor * R_12_Delta);
	double gamma = sqrt(beta_1 * (beta_1 + 2 * beta_12));
	

	return Greeks(beta_1, beta_12, 0., gamma, (beta_1 + beta_12) / gamma);
}

///// 2U


Configuration_2U::Configuration_2U(Piping* _piping) : Configuration(_piping)  
{
	LOG("2U");
}

Greeks Configuration_2U::set_greeks(Piping* piping)
{
	return Greeks();
}

///// CX

Configuration_CX::Configuration_CX(Piping* _piping) : Configuration(_piping)  
{
	LOG("CX");
}

Greeks Configuration_CX::set_greeks(Piping* piping)
{
	return Greeks();
}

}
