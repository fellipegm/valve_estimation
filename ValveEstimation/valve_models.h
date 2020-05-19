#pragma once

#ifndef _VALVE_MODELS_H
#define _VALVE_MODELS_H

#include <vector>
#include <string>
#include <map>
#include "controller.h"


#define PARAM_VALVULA {1.6, 210490, 0.0445, 2550, 0, 0.029, 180000, 41368, 3.571428571428571}
#define PARAM_ATRITO_KANO {25.287220681574890, 1.296780547773071, 11.494727614487891, 0.0, 0.029}
#define PARAM_ATRITO_HE {12.8, 12.2, 11.49, 0.0, 0.029}
#define PARAM_ATRITO_CHOUDHURY {25.287220681574890, 1.296780547773071, 11.494727614487891, 0.0, 0.029}
#define PARAM_ATRITO_KARNOPP {700, 780, 125000, 5.0e-04}
#define PARAM_ATRITO_LUGRE {700, 780, 125000, 5.0e-04, 26000000, 2.039607805437114e+04}
#define PARAM_ATRITO_GMS {700, 780, 125000, 5.0e-04, 29250000, 15600000, 3120000, 2910.90841109400, 1455.45420554700, 4366.36261664100, 0.75, 0.2, 20}
//#define PARAM_ATRITO_GMS {210490, 2550, 709.496, 792.732, 201.992, 0.1, 5.54913e+07, 7.2306e+06, 1.13612e+08, 5e+08, 4.83728e+08, 1.16607e+08, 0.573763, 0.386197, 100}
#define PARAM_ATRITO_SGMS {700, 780, 125000, 5.0e-04, 29250000, 1060000, 0.73, 20}
#define PARAM_ATRITO_GMS1 {700, 780, 125000, 5.0e-04, 29250000, 1e5, 20}

#define GRAPHITE_HE {1.30683999e+01, 1.21924899e+01, 1.41805636e+01, 1.42446524e-03, 2.80129515e-02}
#define GRAPHITE_CHOUDHURY {2.64235916e+01, 8.61440961e-02, 1.45382759e+01, 1.44082678e-03, 2.79368992e-02}
#define GRAPHITE_KANO {2.64235916e+01, 8.61440961e-02, 1.45382758e+01, 1.44082680e-03, 2.79368992e-02}
#define GRAPHITE_KARNOPP {2.06941204e+05, 2.57139941e+03, 6.65423782e+02, 9.42393206e+02, 7.05010949e+04, 3.25785201e-03}
#define GRAPHITE_LUGRE {2.06916034e+05, 2.57207236e+03, 6.69865241e+02, 9.42121876e+02, 7.04812705e+04, 3.22122369e-03, 1.03805024e+09, 2.44167398e+05}
#define GRAPHITE_GMS {2.06769675e+05, 2.58092547e+03, 5.62344054e+02, 9.48677951e+02, 1.74349226e+04, 9.55237056e-02, 7.35429193e+08, 9.10410728e+08, 7.88656394e+08, 9.16629073e+04, 7.36314259e+02, 8.21418019e+05, 6.83524626e-01, 2.85912819e-01, 2.75461126e+02}
#define GRAPHITE_SGMS {2.06637366e+05, 2.58351795e+03, 9.47176012e+02, 9.47422389e+02, 1.82757730e+04, 9.93769567e-02, 1.01185945e+09, 1.03612472e+09, 5.11584968e-01, 1.06951925e+02}
#define GRAPHITE_GMS1 {2.07076853e+05, 2.57426611e+03, 9.43156660e+02, 9.49184643e+02, 1.81991278e+04, 1.00428913e-05, 1.03910639e+09, 1.23555418e+04, 1.08299711e+02}

#define TEFLON_HE {0.16518302, 0.7, 0.30079113, 0.00484964, 0.02393324}
#define TEFLON_CHOUDHURY {2.41594624e-04, 7.65700535e-02, 2.40000000e-17, 4.92463004e-03, 2.40762172e-02}
#define TEFLON_KANO {4.17542911e-03, 8.59583543e-02, 2.63460000e-16, 4.92512454e-03, 2.40767470e-02}
#define TEFLON_KARNOPP {3.89647893e+05, -1.49696401e+00,  3.77675259e-06,  1.31635130e-03, 2.96076838e+00,  6.01146031e-04}
#define TEFLON_LUGRE {3.89653322e+05, -1.49882171e+00,  4.14998423e-01,  5.14757900e-01, 6.81727594e+00,  3.43304880e-04,  3.08426159e+05,  6.50046385e+02}
#define TEFLON_GMS {3.89648815e+05, -1.49697665e+00, 700, 780, 125000, 5.0e-04, 29250000, 15600000, 3120000, 2910.90841109400, 1455.45420554700, 4366.36261664100, 0.75, 0.2, 20}
#define TEFLON_SGMS {3.89648815e+05, -1.49697665e+00, 700, 780, 125000, 5.0e-04, 29250000, 1060000, 0.73, 20}
#define TEFLON_GMS1 {3.89648815e+05, -1.49697665e+00,  6.05097442e-02,  7.75291618e-02, 3.14642711e+00,  8.33723581e-04,  1.75655618e+07,  5.00004020e+02, 1.00000012e+00}



typedef struct simdata
{
	std::vector<double> t;
	std::vector<double> OP;
	std::vector<double> P;
	std::vector<double> x;
	std::vector<double> v;
	std::vector<double> a;
	std::vector<double> F_at;
	std::vector<double> SP;
} simdata;


enum class friction_model {
	kano,
	karnopp,
	lugre,
	gms,
	sgms,
	gms1,
	he,
	choudhury
};


enum sim_type {
	ol,
	cl,
	h_cl
};


class ValveModel {
private:
	std::vector<double> param_valve { 1.6, 210490, }; // param valve has size 9 or 7, depending if the estimated parameters includes k and F_init or not
	std::vector<double> param_friction;

	double Ts = { 1e-3 }; // sampling time
	double dt = { 1e-5 }; // integration time
	std::vector<double> pos0{ 0, -1 }; // initial position
	double t0 = { 0 }; // initial time
	std::vector<double> d0u0{ -1, 0 }; // initial direction and input position the stem stopped, for Kano model only

	// Simulation type
	sim_type simulation_type = ol;

	// Variance for stem position sensor noise
	double std_noise_controller = 1e-20;

	// input data for ol model (diaphragm pressure) or cl model (SP)
	std::vector<double> u;

	// excitation signal for estimating the parameters in closed loop
	std::vector<double> exc_cl;

	// Valve model
	friction_model model = friction_model::kano;

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
	void sim_sgms();
	void sim_gms1();
	void sim_he();
	void sim_choudhury();


	void allocate_sim_data(int len_u);

	bool sim_data_initialized = false;

	double hydraulic_model(double SP, double x, int ct);

public:
	ValveModel();
	
	//ValveModel(const ValveModel& rhs); // copy constructor
	//void operator=(const ValveModel& rhs);

	Controller controller;
	Controller controller_hydraulic;

	void set_valve_param_value(std::vector<double>); // set new values for valve parameters
	void set_friction_param_value(std::vector<double>); // set new values for friction parameters
	void set_sampling_time(double input) { Ts = input; };
	void set_integration_time(double input) { dt = input; };
	void set_d0u0(std::vector<double> input);
	void set_pos0(std::vector<double> input) { pos0 = input; };
	void set_t0(double input) { t0 = input; };
	void set_input_data(std::vector<double>);
	void set_model(friction_model new_model);
	void set_simulation_type(sim_type input) { simulation_type = input; };
	void set_var_noise_controller(double input) { std_noise_controller = input; };


	std::vector<double> get_valve_param_vector(void) { return param_valve; };
	std::vector<double> get_friction_param_vector(void) { return param_friction; };
	double get_sampling_time(void) { return Ts; };
	double get_integration_time(void) { return dt; };
	std::vector<double> get_d0u0(void) { return d0u0; };
	std::vector<double> get_pos0(void) { return pos0; };
	double get_t0(void) { return t0; };
	std::vector<double> get_input_data(void) { return u; };
	friction_model get_model(void) { return model; };
	sim_type get_sim_type(void) { return simulation_type; };
	double get_var_noise_controller(void) { return std_noise_controller; };

	double get_param_friction(size_t i) { return param_friction[i]; };
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

	
	void valve_simulation();

	std::vector<double> OP2P_1order(std::vector<double>* OP);
	std::vector<double> filter2orderZP(const std::vector<double>* data, double wn, double xi);
	std::vector<double> kalman_filter(const std::vector<double>* u, const std::vector<double>* y, double Rv, double Rw);

	std::vector<double> Q_int, Q;

	void clear_sim_data();

	simdata simulation_results;

};
#endif //_VALVE_MODEL_H




