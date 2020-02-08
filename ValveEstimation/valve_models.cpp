

#include <iostream>
#include <math.h>
#include <vector>
#include <string>
#include <complex>

#include "valve_models.h"
#include "utils.h"

void upsample(double x_init, double x_end, int nSamp, double* return_data);


void ValveModel::valve_simulation() {
	switch (model) {
	case kano:
		sim_kano();
		break;
	case karnopp:
		sim_karnopp();
		break;
	case lugre:
		sim_lugre();
		break;
	case gms:
		sim_gms();
		break;
	default:
		std::cout << "Error: Model unindentified. (avalilable options: kano, karnopp, lugre, gms)" << std::endl;
		std::exit;
	}
}

void ValveModel::set_model(friction_model new_model) {
	model = new_model;
	clear_sim_data();
	sim_data_initialized = false;
}

void ValveModel::set_valve_param_value(std::vector<double> input_param_valve) {
	param_valve = input_param_valve;
	if (param_valve.size() == 9) {
		m = param_valve[0];
		k = param_valve[1];
		S_a = param_valve[2];
		F_init = param_valve[3];
		x_min = param_valve[4];
		x_max = param_valve[5];
		p_max = param_valve[6];
		p_min = param_valve[7];
		tau_ip = param_valve[8];
	}
	else if (param_valve.size() == 7) {
		m = param_valve[0];
		S_a = param_valve[1];
		x_min = param_valve[2];
		x_max = param_valve[3];
		p_max = param_valve[4];
		p_min = param_valve[5];
		tau_ip = param_valve[6];
	}
	else
		std::cout << "Error: Valve parameters size has to be 7 or 9" << std::endl;
}

void ValveModel::set_friction_param_value(std::vector<double> input_param_friction) {
	param_friction = input_param_friction;
	switch (model) {
	case kano:
		if (param_friction.size() == 5) {
			S = input_param_friction[0];
			J = input_param_friction[1];
			D = input_param_friction[2];
			x_minP = input_param_friction[3];
			x_maxP = input_param_friction[4];
		}
		else if (param_friction.size() == 3) {
			S = input_param_friction[0];
			J = input_param_friction[1];
			D = input_param_friction[2];
			// Use the same limits of the valve physical constraints
			x_minP = x_min;
			x_maxP = x_max;
		}
		break;
	case karnopp:
		if (param_friction.size() == 4) {
			F_c = input_param_friction[0];
			F_s = input_param_friction[1];
			F_v = input_param_friction[2];
			v_s = input_param_friction[3];
		}
		else if (param_friction.size() == 6) {
			k = input_param_friction[0];
			F_init = input_param_friction[1];
			F_c = input_param_friction[2];
			F_s = input_param_friction[3];
			F_v = input_param_friction[4];
			v_s = input_param_friction[5];
		}
		break;
	case lugre:
		if (param_friction.size() == 6) {
			F_c = input_param_friction[0];
			F_s = input_param_friction[1];
			F_v = input_param_friction[2];
			v_s = input_param_friction[3];
			sigma_0 = input_param_friction[4];
			sigma_1 = input_param_friction[5];
		}
		else if (param_friction.size() == 8) {
			k = input_param_friction[0];
			F_init = input_param_friction[1];
			F_c = input_param_friction[2];
			F_s = input_param_friction[3];
			F_v = input_param_friction[4];
			v_s = input_param_friction[5];
			sigma_0 = input_param_friction[6];
			sigma_1 = input_param_friction[7];
		}
		break;
	case gms:
		if (param_friction.size() == 13) {
			F_c = input_param_friction[0];
			F_s = input_param_friction[1];
			F_v = input_param_friction[2];
			v_s = input_param_friction[3];
			kappa_1 = input_param_friction[4];
			kappa_2 = input_param_friction[5];
			kappa_3 = input_param_friction[6];
			nu_1 = input_param_friction[7];
			nu_2 = input_param_friction[8];
			nu_3 = input_param_friction[9];
			alpha_1 = input_param_friction[10];
			alpha_2 = input_param_friction[11];
			alpha_3 = 1 - alpha_1 - alpha_2;
			C = input_param_friction[12];
			if (alpha_1 < 0 || alpha_1 > 1 || alpha_2 < 0 || alpha_2 > 1 || alpha_1 + alpha_2 > 1)
				std::cout << "Error: GMS alpha parameters error (it has to be between 0 an 1 and alpha_1 + alpha_2 < 1)" << std::endl;
		}
		else if (param_friction.size() == 14) {
			F_c = input_param_friction[0];
			F_s = input_param_friction[1];
			F_v = input_param_friction[2];
			v_s = input_param_friction[3];
			kappa_1 = input_param_friction[4];
			kappa_2 = input_param_friction[5];
			kappa_3 = input_param_friction[6];
			nu_1 = input_param_friction[7];
			nu_2 = input_param_friction[8];
			nu_3 = input_param_friction[9];
			alpha_1 = input_param_friction[10];
			alpha_2 = input_param_friction[11];
			alpha_3 = input_param_friction[12];
			C = input_param_friction[13];
		}
		else if (param_friction.size() == 15) {
			k = input_param_friction[0];
			F_init = input_param_friction[1];
			F_c = input_param_friction[2];
			F_s = input_param_friction[3];
			F_v = input_param_friction[4];
			v_s = input_param_friction[5];
			kappa_1 = input_param_friction[6];
			kappa_2 = input_param_friction[7];
			kappa_3 = input_param_friction[8];
			nu_1 = input_param_friction[9];
			nu_2 = input_param_friction[10];
			nu_3 = input_param_friction[11];
			alpha_1 = input_param_friction[12];
			alpha_2 = input_param_friction[13];
			alpha_3 = 1 - alpha_1 - alpha_2;
			C = input_param_friction[14];
			if (alpha_1 < 0 || alpha_1 > 1 || alpha_2 < 0 || alpha_2 > 1 || alpha_1 + alpha_2 > 1)
				std::cout << "Error: GMS alpha parameters error (it has to be between 0 an 1 and alpha_1 + alpha_2 < 1)" << std::endl;
		}
		else if (param_friction.size() == 16) {
			k = input_param_friction[0];
			F_init = input_param_friction[1];
			F_c = input_param_friction[2];
			F_s = input_param_friction[3];
			F_v = input_param_friction[4];
			v_s = input_param_friction[5];
			kappa_1 = input_param_friction[6];
			kappa_2 = input_param_friction[7];
			kappa_3 = input_param_friction[8];
			nu_1 = input_param_friction[9];
			nu_2 = input_param_friction[10];
			nu_3 = input_param_friction[11];
			alpha_1 = input_param_friction[12];
			alpha_2 = input_param_friction[13];
			alpha_3 = input_param_friction[14];
			C = input_param_friction[15];
		}
		break;
	default:
		std::cout << "Error: Valve parameters size has to be 4 or 6 (Karnopp), 3 or 5 (Kano), 6 or 8 (LuGre) or 13, 14, 15 or 16 (GMS)" << std::endl;
	}
}

void ValveModel::set_d0u0(double* input) {
	d0u0[0] = input[0];
	d0u0[1] = input[1];
}

void ValveModel::set_input_data(std::vector<double> u_input) {
	u = u_input;
}

void ValveModel::sim_kano() {

	double* du = new double[u.size()] {0};
	double* x = new double[u.size()] {0};
	double* input = new double[u.size()] {0};
	double x0 = pos0[0];

	for (int i = 0; i < u.size(); ++i) {
		input[i] = (u[i] - p_min) / (p_max - p_min) * 100;
		if (input[i] < 0)
			input[i] = 0;
		if (input[i] > 100)
			input[i] = 100;
	}

	for (int i = 0; i < u.size(); i++) {
		if (i == 0) {
			input[i] = 0;
			du[i] = 0;
		}
		else {
			input[i] = input[i] - D;
			du[i] = input[i] - input[i - 1];
		}
	}

#ifdef _OLD_KANO_CALC
	for (int i = 0; i < len_u; i++) {
		if (input[i] > 100) {
			input[i] = 100;
		}
		else if (input[i] < 0) {
			input[i] = 0;
		}
	}
#endif // _OLD_KANO_CALC

	int stp = 1;

	double u_s, d;
#ifndef _OLD_KANO_CALC
	u_s = d0u0[1] - D;
#else
	u_s = 0;
#endif
	d = d0u0[0];

#ifndef _OLD_KANO_CALC
	x[0] = (x0 - x_minP) / (x_maxP - x_minP) * ((100 - D) - (S - J) / 2);
#else
	x[0] = x0;
#endif

	for (int i = 1; i < u.size(); i++) {
		if (du[i] * du[i - 1] <= 0.0 && stp == 0) {
			u_s = input[i - 1];
			stp = 1;
		}

		if (stp == 0) {
			x[i] = input[i] - d / 2 * (S - J);
			stp = 0;
		}
		else if ((-d * (input[i] - u_s)) > S) {
			d = -d;
			x[i] = input[i] - d / 2 * (S - J);
			stp = 0;
		}
		else if ((d * (input[i] - u_s)) > J) {
			x[i] = input[i] - d / 2 * (S - J);
			stp = 0;
		}
		else {
			x[i] = x[i - 1];
		}
	}

	allocate_sim_data(u.size());

	for (int i = 0; i < u.size(); i++) {
		simulation_results.t[i] = i * Ts;

#ifndef _OLD_KANO_CALC
		simulation_results.x[i] = x[i] * (x_maxP - x_minP) / ((100 - D) - (S - J) / 2) + x_minP;
#else
		simulation_results.x[i] = x[i] * (x_maxP - x_minP)/100 + x_minP;
#endif

		simulation_results.P[i] = u[i];

		if (simulation_results.x[i] < x_minP) {
			simulation_results.x[i] = x_minP;
		}
		else if (simulation_results.x[i] > x_maxP) {
			simulation_results.x[i] = x_maxP;
		}

		if (i == 0) {
			simulation_results.v[i] = 0;
		}
		else {
			simulation_results.v[i] = (simulation_results.x[i] - simulation_results.x[i - 1]) / Ts;
		}

		simulation_results.a[i] = 0;
		simulation_results.F_at[i] = 0;
	}

	delete[] du;
	delete[] x; 
	delete[] input;
}

void ValveModel::sim_karnopp() {

	if (Ts < dt)
		throw "Error: simulation sampling time has to be lower than data sampling time";

	int nSamp = Ts / dt;

	// Valve model with Karnopp friction model
	// valve stem position initialization
	double* x = new double[nSamp + 2]{ 0 };
	x[nSamp + 1] = pos0[0];
	x[nSamp] = pos0[0];

	// diaphragm pressure initialization
	double* P_exc = new double[nSamp + 2]{ 0 };
	double* P_us = new double[nSamp] {0};

	// valve stem velocity initialization
	double* v = new double[nSamp + 2]{ 0 };

	// valve stem acceleration initialization
	double* a = new double[nSamp + 2]{ 0 };

	// friction force initialization
	double* F_at = new double[nSamp + 2]{ 0 };

	// resultant force initialization
	double* F_res = new double[nSamp + 2]{ 0 };


	allocate_sim_data(u.size());

	int stick = 1;
	double F_r{ 0 }, sig_F{ 0 }, sinal{ 0 }, P_exc_old{ 0 }, P_exc_old_old{ 0 };

	for (int i = 0; i < u.size(); i++) {

		if (i == u.size() - 1) {
			simulation_results.x[i] = x[nSamp + 1];
			simulation_results.v[i] = v[nSamp + 1];
			simulation_results.a[i] = a[nSamp + 1];
			simulation_results.P[i] = u[u.size() - 1];
			simulation_results.F_at[i] = F_at[nSamp + 1];
			simulation_results.t[i] = i * Ts + t0;
			break;
		}

		// reinitialize intermediate vectors
		upsample(u[i], u[i + 1], nSamp, &P_us[0]);
		P_exc_old = P_exc[nSamp + 1];
		P_exc_old_old = P_exc[nSamp];
		if (i == 0) {
			P_exc[0] = u[0];
			P_exc[1] = u[0];
		}
		else {
			P_exc[0] = P_exc_old_old;
			P_exc[1] = P_exc_old;
		}
		for (int ct = 0; ct < nSamp; ct++) {
			P_exc[ct + 2] = P_us[ct];
		}


		x[0] = x[nSamp];
		x[1] = x[nSamp + 1];
		v[0] = v[nSamp];
		v[1] = v[nSamp + 1];
		a[0] = a[nSamp];
		a[1] = a[nSamp + 1];
		F_at[0] = F_at[nSamp];
		F_at[1] = F_at[nSamp + 1];
		F_res[0] = F_res[nSamp];
		F_res[1] = F_res[nSamp + 1];

		simulation_results.x[i] = x[1];
		simulation_results.v[i] = v[1];
		simulation_results.a[i] = a[1];
		simulation_results.P[i] = P_exc[2];
		simulation_results.F_at[i] = F_at[1];
		simulation_results.t[i] = i * Ts + t0;

		for (int j = 2; j < nSamp + 2; j++) {
			if (v[j - 2] == 0 && !stick) {
				stick = 0;
			}
			else if (v[j - 2] * v[j - 1] <= 0) {
				stick = 1;
			}
			else {
				stick = 0;
			}
			if (stick) {
				F_r = S_a * P_exc[j - 1] - k * x[j - 1] - F_init;
				if (F_r > 0) {
					sig_F = 1;
				}
				else if (F_r < 0) {
					sig_F = -1;
				}
				else {
					sig_F = 0;
				}
				F_at[j] = sig_F * fmin(std::abs(F_r), F_s);
				F_res[j] = S_a * P_exc[j - 1] - k * x[j - 1] - F_init - F_at[j];
				v[j - 1] = 0;
				sinal = 0;
				if (std::abs(F_at[j]) >= F_s) {
					stick = 0;
				}
			}
			else {
				if (v[j - 1] > 0) {
					sinal = 1;
				}
				else if (v[j - 1] < 0) {
					sinal = -1;
				}
				else {
					sinal = 0;
				}
				F_at[j] = (F_c + (F_s - F_c) * exp(-pow(v[j - 1] / v_s, 2))) * sinal + F_v * v[j - 1];
				F_res[j] = S_a * P_exc[j - 1] - k * x[j - 1] - F_init - F_at[j];
			}

			a[j] = 1 / m * F_res[j];
			v[j] = dt / 2 * (a[j] + a[j - 1]) + v[j - 1];
			x[j] = dt / 2 * (v[j] + v[j - 1]) + x[j - 1];

			if (x[j] < x_min) {
				F_res[j] = 0;
				a[j] = 0;
				v[j] = 0;
				x[j] = x_min;
			}
			else if (x[j] > x_max) {
				F_res[j] = 0;
				a[j] = 0;
				v[j] = 0;
				x[j] = x_max;
			}
		}

	}
	delete[] x;
	delete[] P_exc;
	delete[] P_us;
	delete[] v;
	delete[] a;
	delete[] F_at;
	delete[] F_res;
}

void ValveModel::sim_lugre() {
	if (Ts < dt)
		throw "Error: simulation sampling time has to be lower than data sampling time";

	int nSamp = Ts / dt;

	// Valve model with Karnopp friction model
	// valve stem position initialization
	double* x = new double[nSamp + 2]{ 0 };
	x[nSamp + 1] = pos0[0];
	x[nSamp] = pos0[0];

	// diaphragm pressure initialization
	double* P_exc = new double[nSamp + 2]{ 0 };
	double* P_us = new double[nSamp] {0};

	// valve stem velocity initialization
	double* v = new double[nSamp + 2]{ 0 };

	// valve stem acceleration initialization
	double* a = new double[nSamp + 2]{ 0 };

	// friction force initialization
	double* F_at = new double[nSamp + 2]{ 0 };

	// resultant force initialization
	double* F_res = new double[nSamp + 2]{ 0 };

	// internal states vector
	double* z = new double[nSamp + 2]{ 0 };
	double z0;

	if (pos0[1] == 1)
		z0 = std::abs(S_a * u[0] - k * pos0[0] - F_init) / sigma_0;
	else if (pos0[1] == -1)
		z0 = -std::abs(S_a * u[0] - k * pos0[0] - F_init) / sigma_0;
	else
		z0 = 0.0;
	z[nSamp + 1] = z0;
	z[nSamp] = z0;

	allocate_sim_data(u.size());

	double P_exc_old{ 0 }, P_exc_old_old{ 0 };

	double g_v{ 0 }, dot_z{ 0 }, dot_z_ant{ 0 };

	for (int i = 0; i < u.size(); i++) {

		if (i == u.size() - 1) {
			simulation_results.x[i] = x[nSamp + 1];
			simulation_results.v[i] = v[nSamp + 1];
			simulation_results.a[i] = a[nSamp + 1];
			simulation_results.P[i] = u[u.size() - 1];
			simulation_results.F_at[i] = F_at[nSamp + 1];
			simulation_results.t[i] = i * Ts + t0;
			break;
		}

		// reinitialize intermediate vectors
		upsample(u[i], u[i + 1], nSamp, &P_us[0]);
		P_exc_old = P_exc[nSamp + 1];
		P_exc_old_old = P_exc[nSamp];
		if (i == 0) {
			P_exc[0] = u[0];
			P_exc[1] = u[0];
		}
		else {
			P_exc[0] = P_exc_old_old;
			P_exc[1] = P_exc_old;
		}
		for (int ct = 0; ct < nSamp; ct++) {
			P_exc[ct + 2] = P_us[ct];
		}

		x[0] = x[nSamp];
		x[1] = x[nSamp + 1];
		v[0] = v[nSamp];
		v[1] = v[nSamp + 1];
		a[0] = a[nSamp];
		a[1] = a[nSamp + 1];
		F_at[0] = F_at[nSamp];
		F_at[1] = F_at[nSamp + 1];
		F_res[0] = F_res[nSamp];
		F_res[1] = F_res[nSamp + 1];
		z[0] = z[nSamp];
		z[1] = z[nSamp + 1];

		simulation_results.x[i] = x[1];
		simulation_results.v[i] = v[1];
		simulation_results.a[i] = a[1];
		simulation_results.P[i] = P_exc[2];
		simulation_results.F_at[i] = F_at[1];
		simulation_results.t[i] = i * Ts + t0;


		for (int j = 2; j < nSamp + 2; j++) {
			g_v = 1 / sigma_0 * (F_c + (F_s - F_c) * exp(-pow(v[j - 1] / v_s, 2)));

			dot_z = v[j - 1] - std::abs(v[j - 1]) / g_v * z[j - 1];

			z[j] = dt / 2 * (dot_z + dot_z_ant) + z[j - 1];

			F_at[j] = sigma_0 * z[j] + sigma_1 * dot_z + F_v * v[j - 1];

			F_res[j] = S_a * P_exc[j - 1] - k * x[j - 1] - F_init - F_at[j];

			a[j] = 1 / m * F_res[j];
			v[j] = dt / 2 * (a[j] + a[j - 1]) + v[j - 1];
			x[j] = dt / 2 * (v[j] + v[j - 1]) + x[j - 1];

			if (x[j] < x_min) {
				F_res[j] = 0;
				a[j] = 0;
				v[j] = 0;
				x[j] = x_min;
			}
			else if (x[j] > x_max) {
				F_res[j] = 0;
				a[j] = 0;
				v[j] = 0;
				x[j] = x_max;
			}

			dot_z_ant = dot_z;
		}

	}
	delete[] x;
	delete[] P_exc;
	delete[] P_us;
	delete[] v;
	delete[] a;
	delete[] F_at;
	delete[] F_res;
	delete[] z;
}

void ValveModel::sim_gms() {
	if (Ts < dt)
		throw "Error: simulation sampling time has to be lower than data sampling time";

	int nSamp = Ts / dt;

	// Valve model with Karnopp friction model
	// valve stem position initialization
	double* x = new double[nSamp + 2] { 0 };
	x[nSamp + 1] = pos0[0];
	x[nSamp] = pos0[0];

	// diaphragm pressure initialization
	double* P_exc = new double[nSamp + 2] { 0 };
	double* P_us = new double[nSamp] { 0 };

	// valve stem velocity initialization
	double* v = new double[nSamp + 2]{ 0 };

	// valve stem acceleration initialization
	double* a = new double[nSamp + 2]{ 0 };

	// friction force initialization
	double* F_at = new double[nSamp + 2]{ 0 };

	// resultant force initialization
	double* F_res = new double[nSamp + 2]{ 0 };

	// internal states vector
	double* z_1 = new double[nSamp + 2]{ 0 };
	double* z_2 = new double[nSamp + 2]{ 0 };
	double* z_3 = new double[nSamp + 2]{ 0 };
	double z_1_0{ 0 }, z_2_0{ 0 }, z_3_0{ 0 };


	if (pos0[1] == 1) {
		z_1_0 = std::abs(S_a * u[0] - k * pos0[0] - F_init) * alpha_1 / kappa_1;
		z_2_0 = std::abs(S_a * u[0] - k * pos0[0] - F_init) * alpha_2 / kappa_2;
		z_3_0 = std::abs(S_a * u[0] - k * pos0[0] - F_init) * alpha_3 / kappa_3;
	}
	else if (pos0[1] == -1) {
		z_1_0 = -std::abs(S_a * u[0] - k * pos0[0] - F_init) * alpha_1 / kappa_1;
		z_2_0 = -std::abs(S_a * u[0] - k * pos0[0] - F_init) * alpha_2 / kappa_2;
		z_3_0 = -std::abs(S_a * u[0] - k * pos0[0] - F_init) * alpha_3 / kappa_3;
	}
	else {
		z_1_0 = 0;
		z_2_0 = 0;
		z_3_0 = 0;
	}
	z_1[nSamp + 1] = z_1_0;
	z_1[nSamp] = z_1_0;
	z_2[nSamp + 1] = z_2_0;
	z_2[nSamp] = z_2_0;
	z_3[nSamp + 1] = z_3_0;
	z_3[nSamp] = z_3_0;

	allocate_sim_data(u.size());

	double P_exc_old{ 0 }, P_exc_old_old{ 0 };

	double F_at_1{ 0 }, F_at_1_ant{ z_1_0 * kappa_1 * alpha_1 }, F_at_2{ 0 }, F_at_2_ant{ z_2_0 * kappa_2 * alpha_2 }, F_at_3{ 0 }, F_at_3_ant{ z_3_0 * kappa_3 * alpha_3 };
	double dot_z_1{ 0 }, dot_z_1_ant{ 0 }, dot_z_2{ 0 }, dot_z_2_ant{ 0 }, dot_z_3{ 0 }, dot_z_3_ant{ 0 };
	double sinal{ 0 }, s_v{ 0 };
	int stick_1{ 1 }, stick_2{ 1 }, stick_3{ 1 };

	for (int i = 0; i < u.size(); i++) {

		if (i == u.size() - 1) {
			simulation_results.x[i] = x[nSamp + 1];
			simulation_results.v[i] = v[nSamp + 1];
			simulation_results.a[i] = a[nSamp + 1];
			simulation_results.P[i] = u[u.size() - 1];
			simulation_results.F_at[i] = F_at[nSamp + 1];
			simulation_results.t[i] = i * Ts + t0;
			break;
		}

		// reinitialize intermediate vectors
		upsample(u[i], u[i + 1], nSamp, &P_us[0]);
		P_exc_old = P_exc[nSamp + 1];
		P_exc_old_old = P_exc[nSamp];
		if (i == 0) {
			P_exc[0] = u[0];
			P_exc[1] = u[0];
		}
		else {
			P_exc[0] = P_exc_old_old;
			P_exc[1] = P_exc_old;
		}
		for (int ct = 0; ct < nSamp; ct++) {
			P_exc[ct + 2] = P_us[ct];
		}

		x[0] = x[nSamp];
		x[1] = x[nSamp + 1];
		v[0] = v[nSamp];
		v[1] = v[nSamp + 1];
		a[0] = a[nSamp];
		a[1] = a[nSamp + 1];
		F_at[0] = F_at[nSamp];
		F_at[1] = F_at[nSamp + 1];
		F_res[0] = F_res[nSamp];
		F_res[1] = F_res[nSamp + 1];
		z_1[0] = z_1[nSamp];
		z_1[1] = z_1[nSamp + 1];
		z_2[0] = z_2[nSamp];
		z_2[1] = z_2[nSamp + 1];
		z_3[0] = z_3[nSamp];
		z_3[1] = z_3[nSamp + 1];

		simulation_results.x[i] = x[1];
		simulation_results.v[i] = v[1];
		simulation_results.a[i] = a[1];
		simulation_results.P[i] = P_exc[2];
		simulation_results.F_at[i] = F_at[1];
		simulation_results.t[i] = i * Ts + t0;


		for (int j = 2; j < nSamp + 2; j++) {
			if (v[j - 1] > 0)
				sinal = 1;
			else if (v[j - 1] < 0)
				sinal = -1;
			else
				sinal = 0;

			if (v[j - 2] == 0.0 && stick_1 == 0)
				stick_1 = 0;
			else if (v[j - 2] * v[j - 1] <= 0)
				stick_1 = 1;

			if (v[j - 2] == 0.0 && stick_2 == 0)
				stick_2 = 0;
			else if (v[j - 2] * v[j - 1] <= 0)
				stick_2 = 1;

			if (v[j - 2] == 0.0 && stick_3 == 0)
				stick_3 = 0;
			else if (v[j - 2] * v[j - 1] <= 0)
				stick_3 = 1;

			s_v = (F_c + (F_s - F_c) * exp(-pow(v[j - 1] / v_s, 2)));

			// First element
			if (stick_1) {
				dot_z_1 = v[j - 1];
				z_1[j] = dt / 2 * (v[j - 1] + v[j - 2]) + z_1[j - 1];
			}
			else {
				dot_z_1 = C * alpha_1 / kappa_1 * (sinal - F_at_1_ant / alpha_1 / s_v);
				z_1[j] = dt / 2 * (dot_z_1 + dot_z_1_ant) + z_1[j - 1];
			}
			F_at_1 = kappa_1 * z_1[j] + nu_1 * dot_z_1;
			if (std::abs(F_at_1) > alpha_1 * s_v)
				stick_1 = 0;

			// Second element
			if (stick_2) {
				dot_z_2 = v[j - 1];
				z_2[j] = dt / 2 * (v[j - 1] + v[j - 2]) + z_2[j - 1];
			}
			else {
				dot_z_2 = C * alpha_2 / kappa_2 * (sinal - F_at_2_ant / alpha_2 / s_v);
				z_2[j] = dt / 2 * (dot_z_2 + dot_z_2_ant) + z_2[j - 1];
			}
			F_at_2 = kappa_2 * z_2[j] + nu_2 * dot_z_2;
			if (std::abs(F_at_2) > alpha_2 * s_v)
				stick_2 = 0;

			// Third element
			if (stick_3) {
				dot_z_3 = v[j - 1];
				z_3[j] = dt / 2 * (v[j - 1] + v[j - 2]) + z_3[j - 1];
			}
			else {
				dot_z_3 = C * alpha_3 / kappa_3 * (sinal - F_at_3_ant / alpha_3 / s_v);
				z_3[j] = dt / 2 * (dot_z_3 + dot_z_3_ant) + z_3[j - 1];
			}
			F_at_3 = kappa_3 * z_3[j] + nu_3 * dot_z_3;
			if (std::abs(F_at_3) > alpha_3 * s_v)
				stick_3 = 0;

			F_at[j] = F_at_1 + F_at_2 + F_at_3 + F_v * v[j - 1];

			F_res[j] = S_a * P_exc[j - 1] - k * x[j - 1] - F_init - F_at[j];

			a[j] = 1 / m * F_res[j];
			v[j] = dt / 2 * (a[j] + a[j - 1]) + v[j - 1];
			x[j] = dt / 2 * (v[j] + v[j - 1]) + x[j - 1];

			if (x[j] < x_min) {
				F_res[j] = 0;
				a[j] = 0;
				v[j] = 0;
				x[j] = x_min;
			}
			else if (x[j] > x_max) {
				F_res[j] = 0;
				a[j] = 0;
				v[j] = 0;
				x[j] = x_max;
			}

			dot_z_1_ant = dot_z_1;
			dot_z_2_ant = dot_z_2;
			dot_z_3_ant = dot_z_3;
			F_at_1_ant = F_at_1;
			F_at_2_ant = F_at_2;
			F_at_3_ant = F_at_3;
		}

	}
	delete[] x;
	delete[] P_exc;
	delete[] P_us;
	delete[] v;
	delete[] a;
	delete[] F_at;
	delete[] F_res;
	delete[] z_1;
	delete[] z_2;
	delete[] z_3;
}

void ValveModel::clear_sim_data() {
	simulation_results.t.clear();
	simulation_results.P.clear();
	simulation_results.x.clear();
	simulation_results.v.clear();
	simulation_results.a.clear();
	simulation_results.F_at.clear();
	simulation_results.t.shrink_to_fit();
	simulation_results.P.shrink_to_fit();
	simulation_results.x.shrink_to_fit();
	simulation_results.v.shrink_to_fit();
	simulation_results.a.shrink_to_fit();
	simulation_results.F_at.shrink_to_fit();
	sim_data_initialized = false;
}

void ValveModel::allocate_sim_data(int len_u) {
	if (!sim_data_initialized) {
		simulation_results.t.reserve(len_u);
		simulation_results.P.reserve(len_u);
		simulation_results.x.reserve(len_u);
		simulation_results.v.reserve(len_u);
		simulation_results.a.reserve(len_u);
		simulation_results.F_at.reserve(len_u);
		for (int i = 0; i < len_u; i++) {
			simulation_results.t.push_back(0.0);
			simulation_results.P.push_back(0.0);
			simulation_results.x.push_back(0.0);
			simulation_results.v.push_back(0.0);
			simulation_results.a.push_back(0.0);
			simulation_results.F_at.push_back(0.0);
		}
		sim_data_initialized = true;
	}
}

std::vector<double> ValveModel::OP2P_1order(std::vector<double>* OP, double tau, double Ts) {
	std::vector<double> P = filter1order(OP, tau, Ts);
	for (int i = 0; i < P.size(); ++i) {
		P[i] = P[i] * (p_max - p_min) / 100 + p_min;
		if (P[i] < p_min)
			P[i] = p_min;
		else if (P[i] > p_max)
			P[i] = p_max;
	}
	return P;
}


std::vector<double> ValveModel::filter2orderZP(const std::vector<double>* data, double wn, double xi) {
	
	double a0 = 1 + 2 * xi * wn * Ts + pow(Ts, 2) * pow(wn, 2);
	double b0 = (pow(Ts, 2) * pow(wn, 2))/a0;
	double a1 = (2 + 2 * xi * wn * Ts)/a0;
	double a2 = -1 / a0;

	std::vector<double> filt_data((*data).size(), 0.0);
	filt_data[0] = ((*data)[0]);
	filt_data[1] = ((*data)[1]);

	for (int i = 2; i < data->size(); ++i)
		filt_data[i] = b0 * (*data)[i - 1] + a1 * filt_data[i - 1] + a2 * filt_data[i - 2];

	std::vector<double> inv_data(filt_data.size(), 0.0), filt2(filt_data.size(), 0.0);
	for (int i = 0; i < filt_data.size(); ++i)
		inv_data[i] = filt_data[filt_data.size()-1-i];

	filt2[0] = filt_data[0];
	filt2[1] = filt_data[1];

	for (int i = 2; i < data->size(); ++i)
		filt2[i] = b0 * inv_data[i - 1] + a1 * filt2[i - 1] + a2 * filt2[i - 2];

	std::vector<double> retdata(filt2.size(), 0.0);
	for (int i = 0; i < filt2.size(); ++i)
		retdata[i] = filt2[filt2.size() - 1 - i];

	return retdata;
}

std::vector<double> ValveModel::kalman_filter(const std::vector<double>* u, const std::vector<double>* y, double Rv, double Rw) {
	
	std::vector<double> P_norm(y->size(), 0.0);
	// Normalize the diaphragm pressure
	for (int i = 0; i < y->size(); ++i)
		P_norm[i] = (*y)[i] - p_min;

	double a = exp(-1 / tau_ip * Ts);
	double b = (p_max - p_min) / 100 * Ts / tau_ip;
	double c = 1;
	
	std::vector<double> M(u->size(), 0.0), P(u->size(), 0.0), Lc(u->size(), 0.0), hat_x(u->size(), 0.0), bar_x(u->size(), 0.0);
	M[0] = 1.0;
	P[0] = 1.0;
	Lc[0] = 0.0;
	hat_x[0] = P_norm[0];
	bar_x[0] = P_norm[0];

	for (int i = 1; i < u->size(); ++i) {
		P[i] = M[i - 1] - M[i - 1] * c / (c * M[i - 1] * c + Rv) * c * M[i - 1];
		Lc[i] = P[i] * c / Rv;
		hat_x[i] = bar_x[i - 1] + Lc[i] * (P_norm[i - 1] - c * bar_x[i - 1]);
		bar_x[i] = a * hat_x[i] + b * (*u)[i - 1];
		M[i] = a * P[i] * a + Rw;
	}

	std::vector<double> retdata(y->size(), 0.0);
	for (int i = 0; i < y->size(); ++i)
		retdata[i] = hat_x[i] + p_min;

	return retdata;
}


void upsample(double x_init, double x_end, int nSamp, double* return_data) {
	for (int i = 0; i < nSamp; i++) {
		return_data[i] = x_init + (x_end - x_init) * i / nSamp;
	}
}