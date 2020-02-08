#pragma once

#ifndef _VALVE_MODELS_H
#define _VALVE_MODELS_H

#include <vector>
#include <string>

typedef struct simdata
{
	std::vector<double> t;
	std::vector<double> P;
	std::vector<double> x;
	std::vector<double> v;
	std::vector<double> a;
	std::vector<double> F_at;
} simdata;


enum friction_model {
	kano,
	karnopp,
	lugre,
	gms
};


class ValveModel {
private:
	std::vector<double> param_valve; // param valve has size 9 or 7, depending if the estimated parameters includes k and F_init or not
	std::vector<double> param_friction;

	double Ts = { 1e-3 }; // sampling time
	double dt = { 1e-5 }; // integration time
	double pos0[2] = { 0, -1 }; // initial position
	double t0 = { 0 }; // initial time
	double d0u0[2] = { -1, 0 }; // initial direction and input position the stem stopped, for Kano model only

	// input data for model
	std::vector<double> u;

	// Valve model
	friction_model model = kano;

	// Valve parameters
	double m{ 1.6 };
	double k{ 210490 };
	double S_a{ 0.0445 };
	double F_init{ 2550 };
	double x_min{ 0 };
	double x_max{ 29e-3 };
	double p_max{ 180000 };
	double p_min{ 41368 };
	double tau_ip{ 3.571428571428571 };

	// Friction parameters
	double S{ 25.287220681574890 };
	double J{ 1.296780547773071 };
	double D{ 11.494727614487891 };
	double x_minP{ 0 };
	double x_maxP{ 29e-3 };
	double F_c{ 700 };
	double F_s{ 780 };
	double F_v{ 125000 };
	double v_s{ 5e-4 };
	double sigma_0{ 26000000 };
	double sigma_1{ 2.039607805437114e+04 };
	double kappa_1{ 29250000 };
	double kappa_2{ 15600000 };
	double kappa_3{ 31200000 };
	double nu_1{ 2910.90841109400 };
	double nu_2{ 1455.45420554700 };
	double nu_3{ 4366.36261664100 };
	double alpha_1{ 0.75 };
	double alpha_2{ 0.2 };
	double alpha_3{ 0.05 };
	double C{ 20 };

	void sim_kano();
	void sim_karnopp();
	void sim_lugre();
	void sim_gms();

	
	void allocate_sim_data(int len_u);

	bool sim_data_initialized = false;

public:
	ValveModel() {};
	void set_valve_param_value(std::vector<double>); // set new values for valve parameters
	void set_friction_param_value(std::vector<double>); // set new values for friction parameters
	double get_param_friction(size_t i) { return param_friction[i]; };
	void set_d0u0(double* input);
	void set_input_data(std::vector<double>);
	void set_model(friction_model new_model);
	void valve_simulation();

	std::vector<double> OP2P_1order(std::vector<double>* OP, double tau, double Ts);
	std::vector<double> filter2orderZP(const std::vector<double>* data, double wn, double xi);
	std::vector<double> kalman_filter(const std::vector<double>* u, const std::vector<double>* y, double Rv, double Rw);

	friction_model get_model() { return model; };
	void clear_sim_data();

	simdata simulation_results;

	double get_Fc() { return F_c; };
	double get_Fs() { return F_s; };
	double get_alpha1() { return alpha_1; };
	double get_alpha2() { return alpha_2; };
	double get_alpha3() { return alpha_3; };

	double get_k() { return k; };
	double get_Finit() { return F_init; };
	double get_pmax() { return p_max; };
	double get_pmin() { return p_min; };
	double get_Sa() { return S_a; };
	double get_xmax() { return x_max; };
	double get_xmin() { return x_min; };
	double get_tauip() { return tau_ip; };

};
#endif //_VALVE_MODEL_CLASS_H




