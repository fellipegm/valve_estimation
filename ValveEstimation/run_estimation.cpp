#include <iostream>
#include <string>
#include <cstring>
#include <ctime>
#include <vector>
#include <direct.h>
#include <windows.h>
#include <tchar.h>


#include "estimator.h"
#include "valve_models.h"
#include "csv_utils.h"
#include "utils.h"


#define PARAM_VALVULA {1.6, 210490, 0.0445, 2550, 0, 0.029, 180000, 41368, 3.571428571428571}
#define PARAM_ATRITO_KANO {25.287220681574890, 1.296780547773071, 11.494727614487891, 0.0, 0.029}
#define PARAM_ATRITO_HE {12.8, 12.2, 11.49, 0.0, 0.029}
#define PARAM_ATRITO_CHOUDHURY {25.287220681574890, 1.296780547773071, 11.494727614487891, 0.0, 0.029}
#define PARAM_ATRITO_KARNOPP {700, 780, 125000, 5.0e-04}
#define PARAM_ATRITO_LUGRE {700, 780, 125000, 5.0e-04, 26000000, 2.039607805437114e+04}
#define PARAM_ATRITO_GMS {700, 780, 125000, 5.0e-04, 29250000, 15600000, 3120000, 2910.90841109400, 1455.45420554700, 4366.36261664100, 0.75, 0.2, 20}
#define PARAM_ATRITO_SGMS {700, 780, 125000, 5.0e-04, 29250000, 1060000, 0.73, 20}
#define PARAM_ATRITO_GMS1 {700, 780, 125000, 5.0e-04, 29250000, 1e5, 20}


#define TS 1e-3


/*
TODO: Verificar nome do computador para usar pastas
TODO: Simulação em malha fechada (gerar dados para estimação)
TODO: Funções para filtrar dados para estimação em malha fechada
TODO: Melhorar a ihm para rodar as otimizações
*/

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

int main() {

	// Get hostname
	TCHAR  comp_name[MAX_COMPUTERNAME_LENGTH + 1];
	DWORD  bufCharCount = MAX_COMPUTERNAME_LENGTH + 1;
	std::string save_dir;
	if (GetComputerName(comp_name, &bufCharCount)) {
		if (_tcscmp(comp_name, "MMI_WIN7-PC") == 0)
			save_dir = "D:\\Drive\\Projetos\\test_data\\";
		else if (_tcscmp(comp_name, "LABSERV") == 0)
			save_dir = "D:\\Users\\fmarques\\Results\\";
		else
			save_dir = "D:\\Drive\\Projetos\\Doutorado2C\\test_data\\";
	}

	std::string type = "estimation"; // estimation, simulation, cl simulation or cl estimation, real estimation, real cl estimation

	if (type.compare("estimation") == 0) {
		// Valve data and configuration of estimation
		double w_n = 0.28;
		double S_exc = 25.7;
		double v_max = 3e-3 * 100 / 29e-3;
		double v_min = 1e-6 * 100 / 29e-3;
		bool sim_noise = false;
		bool estimate_k_finit = true;
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

		std::string experimento = "sinusoidal_velocity";
		// std::vector<std::string> models = {"kano", "he", "choudhury", "karnopp", "lugre", "sgms"};
		std::vector<std::string> models = { "kano", "he", "choudhury", "karnopp", "lugre", "gms1" };
		int n_tests = 20;

		
		for (int i = 0; i < n_tests; ++i) {
			std::vector<double> input;
			if (experimento.compare("sinusoidal_velocity") == 0) {
				std::vector<std::vector<double>> exc = exc_vel_senoidal(v_min, v_max,
					2 * M_PI / w_n / 20, (2 * (100 - 1.1 * S_exc)) / ((v_max - v_min) / 2 + v_min), S_exc, 1 / w_n, 5, 1e-3);
				input = exc[1];
			}
			else if (experimento.compare("aleatory_velocity") == 0) {
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
				estimator.calc_lbub(S_exc, std_k);
				t = clock();
				estimator_output results = estimator.run_estimator();
				t = clock() - t;
				time_taken = ((double)t) / CLOCKS_PER_SEC;
				write_estimation(save_dir, results, model, time_taken);
			}
		}
		std::cout << "Estimation is over" << std::endl;
		std::cin;
	}
	else if (type.compare("cl estimation") == 0) {
		double w_n = 0.28;
		double S_exc = 25.7;
		double v_max = 3e-3 * 100 / 29e-3;
		double v_min = 1e-6 * 100 / 29e-3;

		int n_tests = 20;

		for (int i = 0; i < n_tests; ++i) {
			std::string save_dir("D:\\Drive\\Projetos\\Doutorado2C\\test_data\\");

			double stdNoise = 50.0 / pow(10.0, 25.0 / 10.0);

			std::vector<std::vector<double>> SP = exc_SP_cl_simulation(60, 0.5, 120, 1200, 1e-3);
			std::vector<std::vector<double>> exc = exc_vel_aleatoria(v_min * 10, v_max, S_exc, 3000, 1e-3);

			clock_t t;
			double time_taken;
			ValveModel valve_init_data = ValveModel();

			//	//--------------------------------------------------------------------------------------------
			//	// Kano model
			//	ValveModel valve_kano = ValveModel();
			//	valve_kano.set_model(kano);
			//	valve_kano.set_valve_param_value(PARAM_VALVULA);
			//	valve_kano.set_friction_param_value(PARAM_ATRITO_KANO);
			//	valve_kano.set_simulation_type(h_cl);
			//	valve_kano.set_var_noise_controller(stdNoise);
			//	valve_kano.set_input_data(SP[1]);
			//	valve_kano.controller.set_controller_parameters({ 0.1, 0, 0, 0, 100, -100 });
			//	valve_kano.controller.set_excitation(exc[1]);
			//	valve_kano.controller.set_estimation(true);
			//	valve_kano.valve_simulation();
			//	// Cut data to start at the beginning of the excitation
			//	std::vector<double> P_kano_noise = simulateNoise(valve_kano.simulation_results.P, 25);
			//	procDataCL kano_data = preProcessCLdata(valve_kano.simulation_results.t, valve_kano.simulation_results.OP,
			//		P_kano_noise, valve_kano.simulation_results.x,
			//		valve_kano.controller.t_exc);
			//	//procDataCL kano_data_test = preProcessCLdata(valve_kano.simulation_results.t, valve_kano.simulation_results.OP,
			//	//	valve_kano.simulation_results.P, valve_kano.simulation_results.x,
			//	//	valve_kano.controller.t_exc);
			//	// Add noise to the pressure data
			//	valve_init_data.set_model(kano);
			//	valve_init_data.set_valve_param_value(PARAM_VALVULA);
			//	// Filter data
			//	std::vector<double> P_kano_filt = valve_init_data.kalman_filter(&kano_data.OP, &kano_data.P, kano_data.Rv, kano_data.Rv);
			//	kano_data.P = valve_init_data.filter2orderZP(&P_kano_filt, 1 / valve_init_data.get_tauip() * 100, 0.9);
			//	kano_data.x = valve_init_data.filter2orderZP(&kano_data.x, 1 / valve_init_data.get_tauip() * 100, 0.9);
			//	//std::vector<std::vector<double>> save_data{ kano_data_test.t, kano_data_test.P, kano_data.t, kano_data.P };
			//	//write_matrix(save_dir.append("test_pressure.csv"), &save_data);
			//	//double residual{ 0 };
			//	//for (int i = 0; i < kano_data_test.P.size(); i++) {
			//	//	residual += pow(kano_data_test.P[i] - kano_data.P[i], 2);
			//	//}
			//	// Initialize the Kano model estimation
			//	Estimator estimator_kano = Estimator();
			//	estimator_kano.valve.set_model(kano);
			//	estimator_kano.valve.set_simulation_type(ol);
			//	estimator_kano.valve.set_valve_param_value(PARAM_VALVULA);
			//	estimator_kano.valve.set_input_data(kano_data.P);
			//	estimator_kano.valve.set_d0u0({ kano_data.d0, kano_data.u0 });
			//	estimator_kano.valve.set_pos0({ kano_data.x0, kano_data.d0 });
			//	estimator_kano.set_des_data(kano_data.x);
			//	estimator_kano.calc_lbub(S_exc, 0);
			//	t = clock();
			//	estimator_output results_kano = estimator_kano.run_estimator();
			//	t = clock() - t;
			//	time_taken = ((double)t) / CLOCKS_PER_SEC;
			//	write_estimation(save_dir, results_kano, &estimator_kano, time_taken);
		}

		std::cout << "Estimation is over" << std::endl;
		std::cin;
	}
	//-------------------------------------------------------------------------------------------------------------
	// OL SIMULATION
	else if (type.compare("simulation") == 0) {

		std::string experimento = "aleatory_velocity"; // loaded_data, sinusoidal_velocity, aleatory_velocity

		// Valve data and simulation configuration
		double w_n = 0.28;
		double S_exc = 25.7;
		double v_max = 3e-3 * 100 / 29e-3;
		double v_min = 1e-6 * 100 / 29e-3;

		// std::vector<std::string> models = { "kano", "he", "choudhury", "karnopp", "lugre", "gms", "sgms" };
		std::vector<std::string> models = { "kano", "he", "choudhury", "karnopp", "lugre", "gms", "sgms", "gms1" };

		std::vector<double> input;
		if (experimento.compare("loaded_data") == 0) {
			std::string file("D:\\Drive\\Projetos\\Doutorado2C\\test_data\\data_estimation.csv");
			csv_data dados;
			dados = get_data(file);
			input = dados.P;
		}
		else if (experimento.compare("sinusoidal_velocity") == 0) {

			std::vector<std::vector<double>> exc = exc_vel_senoidal(v_min, v_max,
				2 * M_PI / w_n / 20, (2 * (100 - 1.1 * S_exc)) / ((v_max - v_min) / 2 + v_min), S_exc, 1 / w_n, 5, 1e-3);
			input = exc[1];
		}
		else if (experimento.compare("aleatory_velocity") == 0) {
			std::vector<std::vector<double>> exc = exc_vel_aleatoria(v_min*10, v_max, S_exc, 700, 1e-3);
			input = exc[1];
		}

		// Save the OP value
		// write_vector("D:\\Drive\\Projetos\\Doutorado2C\\test_data\\OP.csv", &input);

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
			if (experimento.compare("loaded_data") == 0)
				valve.set_input_data(input);
			else {
				P = valve.OP2P_1order(&input, 1 / w_n, 1e-3);
				valve.set_input_data(P);
			}
			t = clock();
			valve.valve_simulation();
			t = clock() - t;
			time_taken = ((double)t) / CLOCKS_PER_SEC;
			std::cout << model << " exec time: " << time_taken << std::endl;
			std::string filename = save_dir + "simulation_" + model + ".csv";
			write_simulation(filename, valve.simulation_results);
		}
	}
	else if (type.compare("cl simulation") == 0) {
		std::string experimento = "sinusoidal_velocity"; // loaded_data, sinusoidal_velocity, aleatory_velocity

		// Valve data
		double w_n = 0.28;
		double S_exc = 25.7;
		double v_max = 3e-3 * 100 / 29e-3;
		double v_min = 1e-6 * 100 / 29e-3;

		double stdNoise = 50.0 / pow(10.0, 25.0 / 10.0);

		std::vector<std::vector<double>> SP = exc_SP_cl_simulation(60, 0.25, 120, 1800, 1e-3);
		std::vector<std::vector<double>> exc = exc_vel_aleatoria(v_min * 10, v_max, S_exc, 3000, 1e-3);

		//// Kano model cl simulation
		ValveModel valve = ValveModel();
		//valve.set_model(kano);
		//valve.set_valve_param_value(PARAM_VALVULA);
		//valve.set_friction_param_value(PARAM_ATRITO_KANO);
		//valve.set_simulation_type(h_cl);
		//valve.set_var_noise_controller(stdNoise);
		//valve.set_input_data(SP[1]);
		//valve.controller.set_controller_parameters({0.1, 0, 0, 0, 100, -100});
		//valve.controller.set_excitation(exc[1]);
		//valve.controller.set_estimation(true);
		//valve.valve_simulation();
		//write_simulation("D:\\Drive\\Projetos\\Doutorado2C\\test_data\\simulation_kano_cl.csv", valve.simulation_results);
	}
	else if (type.compare("real estimation") == 0) {

		std::string save_dir("D:\\Drive\\Projetos\\Doutorado2C\\test_data\\");

		std::string valvula{"graphite"};
		std::string test{"sinusoidal"};
		int n_tests = 20;

		std::string open_file = save_dir;
		open_file.append(valvula).append("_data_ol_").append(test).append(".csv");
		csv_real_data dados = get_real_data(open_file);
		
		double m, k, std_k, S_a, F_init, x_min, x_max, x_minP, x_maxP, p_min, p_max, tau_ip, Rv, S0;
		if (valvula.compare("graphite") == 0) {
			m = 1.6;
			k = 2.046708857254743e+05; 
			std_k = 8.682259138949736e+03;
			S_a = 445e-4;
			F_init = 2.619434530232170e+03;
			x_min = 0.0014;
			x_max = 0.02858;
			x_minP = x_min;
			x_maxP = 0.0284;
			p_min = 4.223107229033772e+04;
			p_max = 100 * 1.640909524103048e+03 + p_min;
			tau_ip = 0.933;
			Rv = 7.722764547684506e+04;
			S0 = 23.718774311954462;
		}
		else if (valvula.compare("teflon") == 0) {
			m = 1.6;
			k = 3.839115597423488e+05; 
			std_k = 1.936244943421518e+02;
			S_a = 445e-4;
			F_init = -1.575751572954500;
			x_min = 5.075e-03;
			x_max = 0.02858;
			x_minP = x_min;
			x_maxP = 0.0241;
			p_min = 4.341768573042717e+04;
			p_max = 100 * 1.671912917332373e+03 + p_min;
			tau_ip = 0.979;
			Rv = 3.845874036452569e+04;
			S0 = 0.516156822550698;
		}

		std::vector<double> param_valvula{m, k, S_a, F_init, x_min, x_max, p_max, p_min, tau_ip};

		ValveModel valve_init_data = ValveModel();
		valve_init_data.set_model(friction_model::kano);
		valve_init_data.set_valve_param_value(param_valvula);
		std::vector<double> P_filt = valve_init_data.kalman_filter(&dados.OP, &dados.P, Rv, Rv);
		P_filt = valve_init_data.filter2orderZP(&P_filt, 1 / valve_init_data.get_tauip() * 100.0, 0.9);
		std::vector<double> x_filt = valve_init_data.filter2orderZP(&dados.x, 1 / valve_init_data.get_tauip() * 100.0, 0.9);

		clock_t t;
		double time_taken;
		for (int i = 0; i < n_tests; ++i) {
			//////--------------------------------------------------------------------------------------------
			////// Kano model
			////Estimator estimator_kano = Estimator();
			////estimator_kano.valve.set_model(kano);
			////estimator_kano.valve.set_simulation_type(ol);
			////estimator_kano.valve.set_sampling_time(5e-3);
			////estimator_kano.valve.set_valve_param_value(param_valvula);
			////estimator_kano.valve.set_input_data(P_filt);
			////estimator_kano.set_des_data(x_filt);
			////estimator_kano.calc_lbub(S0, k/3);
			////t = clock();
			////estimator_output results_kano = estimator_kano.run_estimator();
			////t = clock() - t;
			////time_taken = ((double)t) / CLOCKS_PER_SEC;
			////write_estimation(save_dir, results_kano, &estimator_kano, time_taken);

		}

		std::cout << "Estimation is over" << std::endl;
		std::cin;
	}
	return 0;
}