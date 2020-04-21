#include <iostream>
#include <string>
#include <cstring>
#include <ctime>
#include <vector>

#include "estimator.h"
#include "valve_models.h"
#include "csv_utils.h"
#include "utils.h"

#ifdef _SET_WIN7
	#define WINVER 0x0601
	#define _WIN32_WINNT 0x0601
#endif


std::map<std::string, std::vector<double>> model_param_map{
		{"kano", PARAM_ATRITO_KANO},
		{"karnopp", PARAM_ATRITO_KARNOPP},
		{"lugre", PARAM_ATRITO_LUGRE},
		{"gms", PARAM_ATRITO_GMS},
		{"sgms", PARAM_ATRITO_SGMS},
		{"gms1", PARAM_ATRITO_GMS1},
		{"he", PARAM_ATRITO_HE},
		{"choudhury", PARAM_ATRITO_CHOUDHURY}
};

std::map<std::string, std::vector<double>> est_frict_graphite{
		{"kano", GRAPHITE_KANO},
		{"karnopp", GRAPHITE_KARNOPP},
		{"lugre", GRAPHITE_LUGRE},
		{"gms1", GRAPHITE_GMS1},
		{"gms", GRAPHITE_GMS},
		{"sgms", GRAPHITE_SGMS},
		{"he", GRAPHITE_HE},
		{"choudhury", GRAPHITE_CHOUDHURY}
};

std::map<std::string, std::vector<double>> est_frict_teflon{
		{"kano", TEFLON_KANO},
		{"karnopp", TEFLON_KARNOPP},
		{"lugre", TEFLON_LUGRE},
		{"gms1", TEFLON_GMS1},
		{"gms", TEFLON_GMS},
		{"sgms", TEFLON_SGMS},
		{"he", TEFLON_HE},
		{"choudhury", TEFLON_CHOUDHURY}
};

std::map<std::string, friction_model> fric_map{
		{"kano", friction_model::kano},
		{"karnopp", friction_model::karnopp},
		{"lugre", friction_model::lugre},
		{"gms", friction_model::gms},
		{"sgms", friction_model::sgms},
		{"gms1", friction_model::gms1},
		{"he", friction_model::he},
		{"choudhury", friction_model::choudhury}
};


int main(int argc, char** argv) {

	std::string save_dir = "";
	std::string type = "";
	std::string excitation = "";
	std::vector<std::string> models{};
	bool sim_noise = false;
	bool estimate_k_finit = true;
	int n_tests = 1;
	std::string valve_est = "";
	std::string load_file = "";


	int parse_res = parse_arguments(argc, argv, save_dir, type, excitation, models, sim_noise, estimate_k_finit, n_tests, valve_est, load_file);

	if (parse_res) {
		return EXIT_FAILURE;
	}


	if (type.compare("simulation") == 0) {
		// Valve data and simulation configuration
		double w_n = 0.28;
		double S_exc = 25.7;
		double v_max = 3e-3 * 100 / 29e-3;
		double v_min = 1e-6 * 100 / 29e-3;

		std::vector<double> input;
		if (excitation.compare("sinusoidal") == 0) {

			std::vector<std::vector<double>> exc = exc_vel_senoidal(v_min, v_max,
				2 * M_PI / w_n / 20, (2 * (100 - 1.1 * S_exc)) / ((v_max - v_min) / 2 + v_min), S_exc, 1 / w_n, 5, 1e-3);
			input = exc[1];
		}
		else if (excitation.compare("aleatory") == 0) {
			std::vector<std::vector<double>> exc = exc_vel_aleatoria(v_min * 10, v_max, S_exc, 700, 1e-3);
			input = exc[1];
		}

#ifdef _SIM_SAVE_OP
		// Save the OP value
		write_matrix(save_dir + "OP.csv", &input, {"OP"});
#endif

		clock_t t;
		double time_taken;
		ValveModel valve = ValveModel();
		std::vector<double> P;

		for (auto model : models) {
			// Simulate the model
			ValveModel valve = ValveModel();
			valve.set_valve_param_value(PARAM_VALVULA);
			valve.set_model(fric_map[model]);
			valve.set_friction_param_value(model_param_map[model]);
			P = valve.OP2P_1order(&input);
			valve.set_input_data(P);
			t = clock();
			valve.valve_simulation();
			t = clock() - t;
			time_taken = ((double)t) / CLOCKS_PER_SEC;
			std::cout << model << " exec time: " << time_taken << std::endl;
			std::string filename = save_dir + "simulation_" + model + ".sim";
			write_simulation(filename, valve.simulation_results);
		}
	}
	if (type.compare("estimation") == 0) {
		// Valve data and configuration of estimation
		double w_n = 0.28;
		double S_exc = 25.7;
		double v_max = 3e-3 * 100 / 29e-3;
		double v_min = 1e-6 * 100 / 29e-3;
		double std_k;
		double SNR;
		if (sim_noise)
			SNR = 25.0;
		else
			SNR = 1e10;
		if (estimate_k_finit)
			std_k = 210490 / 10;
		else
			std_k = 0;

		for (int i = 0; i < n_tests; ++i) {
			std::vector<double> input;
			if (excitation.compare("sinusoidal") == 0) {
				std::vector<std::vector<double>> exc = exc_vel_senoidal(v_min, v_max,
					2 * M_PI / w_n / 20, (2 * (100 - 1.1 * S_exc)) / ((v_max - v_min) / 2 + v_min), S_exc, 1 / w_n, 5, 1e-3);
				input = exc[1];
			}
			else if (excitation.compare("aleatory") == 0) {
				std::vector<std::vector<double>> exc = exc_vel_aleatoria(v_min * 10, v_max, S_exc, 700, 1e-3);
				input = exc[1];
			}

			// Initialize input data
			ValveModel valve_init_data = ValveModel();
			valve_init_data.set_valve_param_value(PARAM_VALVULA);
			std::vector<double> P = valve_init_data.OP2P_1order(&input);

			double Rv, time_taken;
			clock_t t;

			// Add noise to input data
			std::vector<double> P_noise = simulateNoise(P, SNR);
			Rv = pow(mean_vec(P_noise) / pow(10.0, SNR / 10.0), 2.0);
			// Filter input data
			std::vector<double> P_filt;
			if (sim_noise) {
				P_filt = valve_init_data.kalman_filter(&input, &P_noise, Rv, Rv / 50);
				P_filt = valve_init_data.filter2orderZP(&P_filt, 1 / valve_init_data.get_tauip() * 100.0, 0.9);
			}
			else {
				P_filt = P;
			}

			// Estimate for each model
			for (auto model : models) {
				// Simulate the valve
				ValveModel valve_sim = ValveModel();
				valve_sim.set_model(fric_map[model]);
				valve_sim.set_valve_param_value(PARAM_VALVULA);
				valve_sim.set_friction_param_value(model_param_map[model]);
				valve_sim.set_input_data(P);
				valve_sim.set_simulation_type(ol);
				valve_sim.valve_simulation();

				// Add noise and filter the model output
				std::vector<double> x_noise = simulateNoise(valve_sim.simulation_results.x, SNR);
				std::vector<double> x_filt;
				if (sim_noise) {
					x_filt = valve_sim.filter2orderZP(&x_noise, 1 / valve_init_data.get_tauip() * 100.0, 0.9);
				}
				else {
					x_filt = valve_sim.simulation_results.x;
				}

				// Initialize the model estimation
				Estimator estimator = Estimator();
				estimator.valve.set_model(fric_map[model]);
				estimator.valve.set_simulation_type(ol);
				estimator.valve.set_valve_param_value(PARAM_VALVULA);
				estimator.valve.set_input_data(P_filt);
				estimator.set_des_data(x_filt);
				estimator.calc_lbub(S_exc, std_k, estimate_k_finit);
				t = clock();
				estimator_output results = estimator.run_estimator();
				t = clock() - t;
				time_taken = ((double)t) / CLOCKS_PER_SEC;

				// Evaluate the minimum residual
				valve_sim.set_input_data(P_filt);
				valve_sim.valve_simulation();
				double min_residual{ 0 };
				for (int i = 0; i < valve_sim.simulation_results.x.size(); i ++)
					min_residual += pow(valve_sim.simulation_results.x[i] - x_filt[i], 2.0);
				min_residual = min_residual / ((double)x_filt.size());

				write_estimation(save_dir, results, model, time_taken, min_residual);
			}
		}
		std::cout << "Estimation is over" << std::endl;
		std::cin;
	}
	else if (type.compare("real-estimation") == 0) {

		bool simulate_output = false;

		csv_real_data dados = get_real_data(load_file);

		double m, k, std_k, S_a, F_init, x_min, x_max, x_minP, x_maxP, p_min, p_max, tau_ip, Rv, S0;
		if (valve_est.compare("graphite") == 0) {
			m = 1.6;
			k = 2.046708857254743e+05;
			std_k = 8.682259138949736e+03;
			S_a = 445e-4;
			F_init = 2.619434530232170e+03;
			x_min = 0.0014;
			x_max = 0.02858;
			p_min = 4.223107229033772e+04;
			p_max = 100 * 1.640909524103048e+03 + p_min;
			tau_ip = 0.933;
			Rv = 7.722764547684506e+04;
			S0 = 23.718774311954462;
		}
		else if (valve_est.compare("teflon") == 0) {
			m = 1.6;
			k = 3.839115597423488e+05;
			std_k = 1.936244943421518e+02;
			S_a = 445e-4;
			F_init = -1.575751572954500;
			x_min = 5.0e-03;
			x_max = 0.0284;
			p_min = 4.341768573042717e+04;
			p_max = 100 * 1.671912917332373e+03 + p_min;
			tau_ip = 0.979;
			Rv = 3.845874036452569e+04;
			S0 = 0.516156822550698;
		}

		std::vector<double> param_valvula{ m, k, S_a, F_init, x_min, x_max, p_max, p_min, tau_ip };

		ValveModel valve_init_data = ValveModel();
		valve_init_data.set_model(friction_model::kano);
		valve_init_data.set_valve_param_value(param_valvula);
		std::vector<double> P_filt = valve_init_data.kalman_filter(&dados.OP, &dados.P, Rv, Rv / 50);
		P_filt = valve_init_data.filter2orderZP(&P_filt, 1 / valve_init_data.get_tauip() * 100.0, 0.9);
		std::vector<double> x_filt = valve_init_data.filter2orderZP(&dados.x, 1 / valve_init_data.get_tauip() * 100.0, 0.9);

		std::vector<double> sliced_x = std::vector<double>(x_filt.begin(), x_filt.begin() + 5);
		double initial_pos = mean_vec(sliced_x);

		clock_t t;
		double time_taken;
		for (int i = 0; i < n_tests; ++i) {
			for (auto model : models) {
				// Initialize the model estimation
				Estimator estimator = Estimator();
				estimator.valve.set_sampling_time(5e-3);
				estimator.valve.set_integration_time(5e-6);
				estimator.valve.set_pos0({ initial_pos, -1 });
				estimator.valve.set_valve_param_value(param_valvula);
				estimator.valve.set_model(fric_map[model]);
				estimator.valve.set_simulation_type(ol);
				estimator.valve.set_input_data(P_filt);
				estimator.set_des_data(x_filt);
				estimator.calc_lbub(S0, std_k, estimate_k_finit);
				t = clock();
				estimator_output results = estimator.run_estimator();
				t = clock() - t;
				time_taken = ((double)t) / CLOCKS_PER_SEC;

				write_estimation(save_dir, results, model, time_taken, 0.0);

				if (simulate_output) {
					estimator.valve.set_friction_param_value(results.parameters[0]);
					estimator.valve.valve_simulation();
					std::string filename = save_dir + "simulation_" + model + "_" + valve_est + ".csv";
					write_simulation(filename, estimator.valve.simulation_results);
				}
			}
		}

		std::cout << "Estimation is over" << std::endl;
		std::cin;
	}
	else if (type.compare("real-simulation") == 0) {

		csv_real_data dados = get_real_data(load_file);

		double m, k, std_k, S_a, F_init, x_min, x_max, x_minP, x_maxP, p_min, p_max, tau_ip, Rv, S0;
		std::map<std::string, std::vector<double>> est_fric_param;
		if (valve_est.compare("graphite") == 0) {
			est_fric_param = est_frict_graphite;
			m = 1.6;
			k = 2.046708857254743e+05;
			std_k = 8.682259138949736e+03;
			S_a = 445e-4;
			F_init = 2.619434530232170e+03;
			x_min = 0.0014;
			x_max = 0.02858;
			p_min = 4.223107229033772e+04;
			p_max = 100 * 1.640909524103048e+03 + p_min;
			tau_ip = 0.933;
			Rv = 7.722764547684506e+04;
			S0 = 23.718774311954462;
		}
		else {
			est_fric_param = est_frict_teflon;
			m = 1.6;
			k = 3.839115597423488e+05;
			std_k = 1.936244943421518e+02;
			S_a = 445e-4;
			F_init = -1.575751572954500;
			x_min = 5.0e-03;
			x_max = 0.0284;
			p_min = 4.341768573042717e+04;
			p_max = 100 * 1.671912917332373e+03 + p_min;
			tau_ip = 0.979;
			Rv = 3.845874036452569e+04;
			S0 = 0.516156822550698;
		}

		std::vector<double> param_valvula{ m, k, S_a, F_init, x_min, x_max, p_max, p_min, tau_ip };

		ValveModel valve_init_data = ValveModel();
		valve_init_data.set_model(friction_model::kano);
		valve_init_data.set_valve_param_value(param_valvula);
		std::vector<double> P_filt = valve_init_data.kalman_filter(&dados.OP, &dados.P, Rv, Rv / 50);
		P_filt = valve_init_data.filter2orderZP(&P_filt, 1 / valve_init_data.get_tauip() * 100.0, 0.9);
		std::vector<double> x_filt = valve_init_data.filter2orderZP(&dados.x, 1 / valve_init_data.get_tauip() * 100.0, 0.9);

		std::vector<double> sliced_x = std::vector<double>(x_filt.begin(), x_filt.begin() + 5);
		double initial_pos = mean_vec(sliced_x);

		std::vector<std::vector<double>> save_filter{ dados.P, dados.x, P_filt, x_filt };
		write_matrix(save_dir + "filtered_data_" + valve_est + "_" + excitation + ".csv", save_filter, {"P", "x", "P_filt", "x_filt"});

		// Initialize the valve model
		for (auto model : models) {
			ValveModel valve = ValveModel();
			valve.set_sampling_time(5e-3);
			valve.set_integration_time(5e-6);
			valve.set_pos0({ initial_pos, -1 });
			valve.set_valve_param_value(param_valvula);
			valve.set_model(fric_map[model]);
			valve.set_simulation_type(ol);
			valve.set_input_data(P_filt);
			valve.set_friction_param_value(est_fric_param[model]);
			valve.valve_simulation();

			std::string filename = save_dir + "simulation_" + model + "_" + valve_est + "_" + excitation + ".sim";
			write_simulation(filename, valve.simulation_results);
		}

		std::cout << "Real simulation is over" << std::endl;
		std::cin;
	}
	else if (type.compare("cl-estimation") == 0) {
		//double w_n = 0.28;
		//double S_exc = 25.7;
		//double v_max = 3e-3 * 100 / 29e-3;
		//double v_min = 1e-6 * 100 / 29e-3;

		//int n_tests = 20;

		//for (int i = 0; i < n_tests; ++i) {
		//	std::string save_dir("D:\\Drive\\Projetos\\Doutorado2C\\test_data\\");

		//	double stdNoise = 50.0 / pow(10.0, 25.0 / 10.0);

		//	std::vector<std::vector<double>> SP = exc_SP_cl_simulation(60, 0.5, 120, 1200, 1e-3);
		//	std::vector<std::vector<double>> exc = exc_vel_aleatoria(v_min * 10, v_max, S_exc, 3000, 1e-3);

		//	clock_t t;
		//	double time_taken;
		//	ValveModel valve_init_data = ValveModel();

		//	//	//--------------------------------------------------------------------------------------------
		//	//	// Kano model
		//	//	ValveModel valve_kano = ValveModel();
		//	//	valve_kano.set_model(kano);
		//	//	valve_kano.set_valve_param_value(PARAM_VALVULA);
		//	//	valve_kano.set_friction_param_value(PARAM_ATRITO_KANO);
		//	//	valve_kano.set_simulation_type(h_cl);
		//	//	valve_kano.set_var_noise_controller(stdNoise);
		//	//	valve_kano.set_input_data(SP[1]);
		//	//	valve_kano.controller.set_controller_parameters({ 0.1, 0, 0, 0, 100, -100 });
		//	//	valve_kano.controller.set_excitation(exc[1]);
		//	//	valve_kano.controller.set_estimation(true);
		//	//	valve_kano.valve_simulation();
		//	//	// Cut data to start at the beginning of the excitation
		//	//	std::vector<double> P_kano_noise = simulateNoise(valve_kano.simulation_results.P, 25);
		//	//	procDataCL kano_data = preProcessCLdata(valve_kano.simulation_results.t, valve_kano.simulation_results.OP,
		//	//		P_kano_noise, valve_kano.simulation_results.x,
		//	//		valve_kano.controller.t_exc);
		//	//	//procDataCL kano_data_test = preProcessCLdata(valve_kano.simulation_results.t, valve_kano.simulation_results.OP,
		//	//	//	valve_kano.simulation_results.P, valve_kano.simulation_results.x,
		//	//	//	valve_kano.controller.t_exc);
		//	//	// Add noise to the pressure data
		//	//	valve_init_data.set_model(kano);
		//	//	valve_init_data.set_valve_param_value(PARAM_VALVULA);
		//	//	// Filter data
		//	//	std::vector<double> P_kano_filt = valve_init_data.kalman_filter(&kano_data.OP, &kano_data.P, kano_data.Rv, kano_data.Rv);
		//	//	kano_data.P = valve_init_data.filter2orderZP(&P_kano_filt, 1 / valve_init_data.get_tauip() * 100, 0.9);
		//	//	kano_data.x = valve_init_data.filter2orderZP(&kano_data.x, 1 / valve_init_data.get_tauip() * 100, 0.9);
		//	//	//std::vector<std::vector<double>> save_data{ kano_data_test.t, kano_data_test.P, kano_data.t, kano_data.P };
		//	//	//write_matrix(save_dir.append("test_pressure.csv"), &save_data);
		//	//	//double residual{ 0 };
		//	//	//for (int i = 0; i < kano_data_test.P.size(); i++) {
		//	//	//	residual += pow(kano_data_test.P[i] - kano_data.P[i], 2);
		//	//	//}
		//	//	// Initialize the Kano model estimation
		//	//	Estimator estimator_kano = Estimator();
		//	//	estimator_kano.valve.set_model(kano);
		//	//	estimator_kano.valve.set_simulation_type(ol);
		//	//	estimator_kano.valve.set_valve_param_value(PARAM_VALVULA);
		//	//	estimator_kano.valve.set_input_data(kano_data.P);
		//	//	estimator_kano.valve.set_d0u0({ kano_data.d0, kano_data.u0 });
		//	//	estimator_kano.valve.set_pos0({ kano_data.x0, kano_data.d0 });
		//	//	estimator_kano.set_des_data(kano_data.x);
		//	//	estimator_kano.calc_lbub(S_exc, 0);
		//	//	t = clock();
		//	//	estimator_output results_kano = estimator_kano.run_estimator();
		//	//	t = clock() - t;
		//	//	time_taken = ((double)t) / CLOCKS_PER_SEC;
		//	//	write_estimation(save_dir, results_kano, &estimator_kano, time_taken);
		//}

		//std::cout << "Estimation is over" << std::endl;
		//std::cin;
	}
	else if (type.compare("cl-simulation") == 0) {
		//std::string experimento = "sinusoidal_velocity"; // loaded_data, sinusoidal_velocity, aleatory_velocity

		//// Valve data
		//double w_n = 0.28;
		//double S_exc = 25.7;
		//double v_max = 3e-3 * 100 / 29e-3;
		//double v_min = 1e-6 * 100 / 29e-3;

		//double stdNoise = 50.0 / pow(10.0, 25.0 / 10.0);

		//std::vector<std::vector<double>> SP = exc_SP_cl_simulation(60, 0.25, 120, 1800, 1e-3);
		//std::vector<std::vector<double>> exc = exc_vel_aleatoria(v_min * 10, v_max, S_exc, 3000, 1e-3);

		////// Kano model cl simulation
		//ValveModel valve = ValveModel();
		////valve.set_model(kano);
		////valve.set_valve_param_value(PARAM_VALVULA);
		////valve.set_friction_param_value(PARAM_ATRITO_KANO);
		////valve.set_simulation_type(h_cl);
		////valve.set_var_noise_controller(stdNoise);
		////valve.set_input_data(SP[1]);
		////valve.controller.set_controller_parameters({0.1, 0, 0, 0, 100, -100});
		////valve.controller.set_excitation(exc[1]);
		////valve.controller.set_estimation(true);
		////valve.valve_simulation();
		////write_simulation("D:\\Drive\\Projetos\\Doutorado2C\\test_data\\simulation_kano_cl.csv", valve.simulation_results);
	}

	return 0;
}