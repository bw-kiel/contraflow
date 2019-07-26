#include "configuration.h"
#include "utilities.h"
#include "fluid.h"
#include "greeks.h"
#include <math.h>

#define _SQ2 1.41421356237
#define _2PI 6.283185307179586

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

	piping->Nu_0 = piping->fluid.Nusselt_pipe(piping->Re_0, L, piping->d_0_i);
	piping->Nu_1 = piping->Nu_0;
	
	LOG("u:  " << piping->u_0);
	LOG("Re: " << piping->Re_0);
	LOG("Nu: " << piping->Nu_0);
}

Resistances Configuration_U::set_resistances(double D, double lambda_g)
{
	double D2 = D * D;
	double d2 = piping->d_0_o * piping->d_0_o;
	double w2 = piping->w * piping->w;
	double l = _2PI * lambda_g;

	double x = log(sqrt(D2 + 2 * d2)/(2*piping->d_0_o)) / 
		log(D/(_SQ2*piping->d_0_o));
	double R_g = (1.601 - 0.888 * piping->w / D) * acosh((D2 + d2 - w2)/(2 * D * piping->d_0_o)) / l;
	double R_ar = acosh((2*w2 - d2)/d2) / l;

	R_adv = 1 / (piping->Nu_0 * piping->fluid.get_lambda() * M_PI);
	R_con_a = log(piping->d_0_o / piping->d_0_i) /(piping->lambda_0 * _2PI);
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

void Configuration_U::set_functions(double& f1, double& f2, double& f3, const double& gz, const Greeks& greeks)
{
	const double ch = cosh(gz);
	const double sh = sinh(gz);

	f1 = ch - sh * greeks.get_delta();
	f2 = sh * greeks.get_beta_12() / greeks.get_gamma(); 
	f3 = ch + sh * greeks.get_delta();
}

double Configuration_U::F4(const double &z, const double &a, const double &b, const Greeks& greeks)
{
	const double gamma = greeks.get_gamma();
	const double beta_1 = greeks.get_beta_1();
	const double beta_12 = greeks.get_beta_12();
	const double delta = greeks.get_delta();

        return  (sinh(gamma*(z-a)) - sinh(gamma*(z-b))) * beta_1 / gamma + \
				(cosh(gamma*(z-b)) - cosh(gamma*(z-a))) *
                                (beta_1 / gamma) * (delta + beta_12 / gamma);
}

double Configuration_U::F5(const double &z, const double &a, const double &b, const Greeks& greeks)
{

	const double gamma = greeks.get_gamma();
	const double beta_1 = greeks.get_beta_1();
	const double beta_12 = greeks.get_beta_12();
	const double delta = greeks.get_delta();

        return (sinh(gamma*(z-a)) - sinh(gamma*(z-b))) * beta_1 / gamma \ 
				- (cosh(gamma*(z-b)) - cosh(gamma*(z-a))) *
                                (beta_1 / gamma) * (delta + beta_12 / gamma);
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

void Configuration_2U::set_functions(double& f1, double& f2, double& f3, const double& gz, const Greeks& greeks)
{
	const double ch = cosh(gz);
	const double sh = sinh(gz);

	f1 = ch - sh * greeks.get_delta();
	f2 = sh * greeks.get_beta_12() / greeks.get_gamma(); 
	f3 = ch + sh * greeks.get_delta();
}

double Configuration_2U::F4(const double &z, const double &a, const double &b, const Greeks& greeks)
{

	return 0.;
}

double Configuration_2U::F5(const double &z, const double &a, const double &b, const Greeks& greeks)
{

	return 0.;
}

///// CX

Configuration_CX::Configuration_CX(Piping* _piping) : Configuration(_piping)  
{
	LOG("CX");
	piping->A_0 = piping->d_0_i * piping->d_0_i * M_PI / 4; 
	piping->A_1 = (piping->d_1_i * piping->d_1_i - piping->d_0_o * piping->d_0_o) * M_PI / 4;  
	// LOG("A_0: " << piping->A_0);
	// LOG("A_1: " << piping->A_1);
}

Greeks Configuration_CX::set_greeks(Piping* piping)
{
	double factor = piping->get_fluid().get_c_vol() * piping->u_0 * piping->A_0;
	double beta_1 = 1. / (factor * R_1_Delta);
	double beta_12 = 1. / (factor * R_12_Delta);

	double beta = - beta_1 / 2;
	double gamma = sqrt(beta_1 * (beta_1 * .25 + beta_12));

	return Greeks(beta_1, beta_12, beta, gamma, (beta_1*.5 + beta_12) / gamma);
}

void Configuration_CX::set_functions(double& f1, double& f2, double& f3, const double& gz, const Greeks& greeks)
{
	const double ch = cosh(gz);
	const double sh = sinh(gz);
	const double ebz = exp(greeks.get_beta() * gz / greeks.get_gamma());

	f1 = ebz * (ch - sh * greeks.get_delta());
	f2 = ebz * (sh * greeks.get_beta_12() / greeks.get_gamma()); 
	f3 = ebz * (ch + sh * greeks.get_delta());
}

double Configuration_CX::F4(const double &z, const double &a, const double &b, const Greeks& greeks)
{
	const double gamma = greeks.get_gamma();
	const double beta_1 = greeks.get_beta_1();
	const double beta_12 = greeks.get_beta_12();
	const double beta = greeks.get_beta();
	const double delta = greeks.get_delta();

	double ex = (exp(beta*(z-b)) - exp(beta*(z-a))) * beta_1 / (gamma*gamma - beta*beta);
	double co = (cosh(gamma*(z-b)) - cosh(gamma*(z-a))) * (gamma * delta + beta); 
	double si = (sinh(gamma*(z-a)) - sinh(gamma*(z-b))) * (gamma + delta*beta);
	
	return (co + si) * ex;
}

double Configuration_CX::F5(const double &z, const double &a, const double &b, const Greeks& greeks)
{

	const double gamma = greeks.get_gamma();
	const double beta_1 = greeks.get_beta_1();
	const double beta_12 = greeks.get_beta_12();
	const double beta = greeks.get_beta();
	const double delta = greeks.get_delta();

	double ex = (exp(beta*(z-b)) - exp(beta*(z-a))) * beta_1 *beta_12 / (gamma*gamma - beta*beta);
	double si = (sinh(gamma*(z-b)) - sinh(gamma*(z-a))) * beta / gamma; 
	double co = - (cosh(gamma*(z-b)) - cosh(gamma*(z-a)));
	return ex*(si+co);
}

void Configuration_CX::set_flow(double L)
{
	piping->u_0 = piping->Q / piping->A_0; 
	piping->u_1 = piping->Q / piping->A_1; 

	piping->Re_0 = piping->fluid.Reynolds(piping->u_0, piping->d_0_i);
	piping->Re_1 = piping->fluid.Reynolds(piping->u_1, piping->d_1_i - piping->d_0_o);

	piping->Nu_0 = piping->fluid.Nusselt_pipe(piping->Re_0, L, piping->d_0_i);
	piping->Nu_1 = piping->fluid.Nusselt_ring(piping->Re_1, L, piping->d_0_o, piping->d_1_i);
	
	/*LOG("u_0:  " << piping->u_0);
	LOG("u_1:  " << piping->u_1);
	LOG("Re_0: " << piping->Re_0);
	LOG("Re_1: " << piping->Re_1);
	LOG("Nu_0: " << piping->Nu_0);
	LOG("Nu_1: " << piping->Nu_1);
	*/
}
Resistances Configuration_CX::set_resistances(double D, double lambda_g)
{
	const double d = piping->d_1_i - piping->d_0_o;
	const double lnDd = log(D/(piping->d_1_o));
	const double lapi = piping->fluid.get_lambda() * M_PI;

	const double x = log(sqrt(D*D + piping->d_1_o * piping->d_1_o)/(_SQ2 * piping->d_1_o)) / lnDd;
	const double R_g = lnDd / (_2PI * lambda_g);

	R_adv_0 = 1 / (piping->Nu_0 * lapi);
	R_adv_1 = 1 / (piping->Nu_1 * lapi) * d / piping->d_0_o;
	R_adv_2 = 1 / (piping->Nu_1 * lapi) * d / piping->d_1_i;

	R_con_0_a = log(piping->d_0_o / piping->d_0_i) /(piping->lambda_0 * _2PI);
	R_con_1_a = log(piping->d_1_o / piping->d_1_i) /(piping->lambda_1 * _2PI);
	R_con_b = x * R_g;

	R_ff = R_adv_0 + R_adv_1 + R_con_0_a;
	R_fg = R_adv_2 + R_con_1_a + R_con_b;
	R_gs = (1 - x) * R_g;

	R_1_Delta = R_fg + R_gs;
	R_2_Delta = 0.;
	R_12_Delta = R_ff;

	/*LOG("x:          " << x);
	LOG("R_g:        " << R_g);

	LOG("R_adv_0:      " << R_adv_0);
	LOG("R_adv_1:      " << R_adv_1);
	LOG("R_adv_2:      " << R_adv_2);
	LOG("R_con_0_a:    " << R_con_0_a);
	LOG("R_con_1_a:    " << R_con_1_a);
	LOG("R_con_b:      " << R_con_b);
	LOG("R_ff:         " << R_ff);
	LOG("R_fg:         " << R_fg);
	LOG("R_gs:         " << R_gs);
	LOG("R_1_Delta:  " << R_1_Delta);
	LOG("R_2_Delta:  " << R_2_Delta);
	LOG("R_12_Delta: " << R_12_Delta);
	*/
	return {R_1_Delta, R_2_Delta};
}

}
