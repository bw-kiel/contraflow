#ifndef CONFIGURATION_H
#define CONFIGURATION_H

#include "piping.h"
#include "greeks.h"

namespace contra
{

class Configuration
{
public:
	Configuration(Piping* _piping) : piping(_piping) {}
	virtual Resistances set_resistances(double D, double lambda_g) = 0;
	virtual void set_flow(double L) = 0;
	virtual Greeks set_greeks(Piping* piping) = 0;

	double get_R_1_Delta() { return R_1_Delta; }
	double get_R_12_Delta() { return R_12_Delta; }

protected:
	Piping* piping;
	double R_1_Delta;
	double R_2_Delta;
	double R_12_Delta;

private:

};

class Configuration_U : public Configuration
{
public:
	Configuration_U(Piping* _piping);
	//void set_resistances_pipe(){}
	Resistances set_resistances(double D, double lambda_g);

	void set_flow(double L);
	Greeks set_greeks(Piping* piping);

private:
	double R_adv;	// advective resistance
	double R_con_a;
	double R_con_b;

	double R_gs;
	double R_fg;
	double R_gg;


};

class Configuration_2U : public Configuration
{
public:
	Configuration_2U(Piping* _piping);
	//void set_resistances_pipe(){}
	Resistances set_resistances(double D, double lambda_g)
	{
		//set_resistances_pipe();
	}

	void set_flow(double L) {}
	Greeks set_greeks(Piping* piping);

private:
	double R_gs;
	double R_con0_a, R_con_b;
	double R_fg;
	double R_gg1;
	double R_gg2;
	double R_adv;	// advective resistance
};

class Configuration_CX : public Configuration
{
public:
	Configuration_CX(Piping* _piping);

	//void set_resistances_pipe(){}
	//void set_resistances_ring(){}
	Resistances set_resistances(double D, double lambda_g) 
	{
		//set_resistances_pipe();
		//set_resistances_ring();
	}

	void set_flow(double L) {}
	Greeks set_greeks(Piping* piping);
private:
	double R_gs;
	double R_con0_a, R_con_b;
	double R_ff;
	double R_fg;
	double R_con1_a;

	double R_adv_0;	// advective resistance - pipe 0
	double R_adv_1;	// advective resistance - pipe 1
	double R_adv_2;	// advective resistance
};

}

#endif
