#include <iostream>
#include <string>
#include <cstring>
#include <ctime>
#include <vector>
#include <direct.h>

#include "estimator.h"
#include "valve_models.h"
#include "csv_utils.h"
#include "utils.h"


#define PARAM_VALVULA {1.6, 210490, 0.0445, 2550, 0, 0.029, 180000, 41368, 3.571428571428571}
#define PARAM_ATRITO_KANO {25.287220681574890, 1.296780547773071, 11.494727614487891, 0.0, 0.029}
#define PARAM_ATRITO_HE {12.6, 12.2, 11.49, 0.0, 0.029}
#define PARAM_ATRITO_CHOUDHURY {25.287220681574890, 1.296780547773071, 11.494727614487891, 0.0, 0.029}
#define PARAM_ATRITO_KARNOPP {700, 780, 125000, 5.0e-04}
#define PARAM_ATRITO_LUGRE {700, 780, 125000, 5.0e-04, 26000000, 2.039607805437114e+04}
#define PARAM_ATRITO_GMS {700, 780, 125000, 5.0e-04, 29250000, 15600000, 31200000, 2910.90841109400, 1455.45420554700, 4366.36261664100, 0.75, 0.2, 20}

#define TS 1e-3


/*
TODO: Verificar nome do computador para usar pastas
TODO: Simulação em malha fechada (gerar dados para estimação)
TODO: Funções para filtrar dados para estimação em malha fechada
TODO: melhorar a ihm para rodar as otimizações
*/


int main() {
	std::string type = "cl estimation"; // estimation, simulation, cl simulation or cl estimation

	if (type.compare("estimation") == 0) {
		//std::string file("D:\\Doc2C\\test_data\\data_estimation.csv");
		std::string estimation_data("D:\\Drive\\Projetos\\Doutorado2C\\test_data\\data_estimation_sin.csv");
		std::string save_dir("D:\\Drive\\Projetos\\Doutorado2C\\test_data\\");

		csv_data dados;

		double S0 = 25.7;
		int n_tests = 20;
		bool simulate_noise = true;

		for (int i = 0; i < n_tests; ++i) {
			if (simulate_noise) {
				dados = get_data(estimation_data);
				std::vector<double> P_noise = simulateNoise(dados.P, 25);
				std::vector<double> x_kano_noise = simulateNoise(dados.x_kano, 25);
				std::vector<double> x_karnopp_noise = simulateNoise(dados.x_karnopp, 25);
				std::vector<double> x_lugre_noise = simulateNoise(dados.x_lugre, 25);
				std::vector<double> x_gms_noise = simulateNoise(dados.x_gms, 25);
				std::vector<double> x_he_noise = simulateNoise(dados.x_he, 25);
				std::vector<double> x_choudhury_noise = simulateNoise(dados.x_choudhury, 25);

				ValveModel valve_init_data = ValveModel();
				valve_init_data.set_model(kano);
				valve_init_data.set_valve_param_value(PARAM_VALVULA);

				double var = pow(mean_vec(dados.P) / pow(10, 25 / 10), 2);
				std::vector<double> P_filt = valve_init_data.kalman_filter(&dados.OP, &P_noise, var, var/50);
				dados.P = valve_init_data.filter2orderZP(&P_filt, 1/valve_init_data.get_tauip()*100, 0.9);
				
				dados.x_kano = valve_init_data.filter2orderZP(&x_kano_noise, 1 / valve_init_data.get_tauip() * 100, 0.9);
				dados.x_karnopp = valve_init_data.filter2orderZP(&x_karnopp_noise, 1 / valve_init_data.get_tauip() * 100, 0.9);
				dados.x_lugre = valve_init_data.filter2orderZP(&x_lugre_noise, 1 / valve_init_data.get_tauip() * 100, 0.9);
				dados.x_gms = valve_init_data.filter2orderZP(&x_gms_noise, 1 / valve_init_data.get_tauip() * 100, 0.9);
				dados.x_he = valve_init_data.filter2orderZP(&x_he_noise, 1 / valve_init_data.get_tauip() * 100, 0.9);
				dados.x_choudhury = valve_init_data.filter2orderZP(&x_choudhury_noise, 1 / valve_init_data.get_tauip() * 100, 0.9);
			}

			//// Initialize the Kano model estimation
			//Estimator estimator_kano = Estimator();
			//estimator_kano.valve.set_model(kano);
			//estimator_kano.valve.set_valve_param_value(PARAM_VALVULA);
			//estimator_kano.valve.set_input_data(dados.P);
			//estimator_kano.set_des_data(dados.x_kano);
			//estimator_kano.calc_lbub(S0, 0);
			//clock_t t;
			//t = clock();
			//estimator_output results_kano = estimator_kano.run_estimator();
			//t = clock() - t;
			//double time_taken = ((double)t) / CLOCKS_PER_SEC;
			//write_estimation(save_dir, results_kano, &estimator_kano, time_taken);


			//// Initialize the Karnopp model estimation
			//Estimator estimator_karnopp = Estimator();
			//estimator_karnopp.valve.set_model(karnopp);
			//estimator_karnopp.valve.set_valve_param_value(PARAM_VALVULA);
			//estimator_karnopp.valve.set_input_data(dados.P);
			//estimator_karnopp.set_des_data(dados.x_karnopp);
			//estimator_karnopp.calc_lbub(S0, 0);
			//t = clock();
			//estimator_output results_karnopp = estimator_karnopp.run_estimator();
			//t = clock() - t;
			//time_taken = ((double)t) / CLOCKS_PER_SEC;
			//write_estimation(save_dir, results_karnopp, &estimator_karnopp, time_taken);


			//// Initialize the LuGre model estimation
			//Estimator estimator_lugre = Estimator();
			//estimator_lugre.valve.set_model(lugre);
			//estimator_lugre.valve.set_valve_param_value(PARAM_VALVULA);
			//estimator_lugre.valve.set_input_data(dados.P);
			//estimator_lugre.set_des_data(dados.x_lugre);
			//estimator_lugre.calc_lbub(S0, 0);
			//t = clock();
			//estimator_output results_lugre = estimator_lugre.run_estimator();
			//t = clock() - t;
			//time_taken = ((double)t) / CLOCKS_PER_SEC;
			//write_estimation(save_dir, results_lugre, &estimator_lugre, time_taken);


			// Initialize the He model estimation
			Estimator estimator_he = Estimator();
			estimator_he.valve.set_model(he);
			estimator_he.valve.set_valve_param_value(PARAM_VALVULA);
			estimator_he.valve.set_input_data(dados.P);
			estimator_he.set_des_data(dados.x_he);
			estimator_he.calc_lbub(S0, 0);
			clock_t t = clock();
			estimator_output results_he = estimator_he.run_estimator();
			t = clock() - t;
			double time_taken = ((double)t) / CLOCKS_PER_SEC;
			write_estimation(save_dir, results_he, &estimator_he, time_taken);


			// Initialize the Choudhury model estimation
			Estimator estimator_choudhury = Estimator();
			estimator_choudhury.valve.set_model(choudhury);
			estimator_choudhury.valve.set_valve_param_value(PARAM_VALVULA);
			estimator_choudhury.valve.set_input_data(dados.P);
			estimator_choudhury.set_des_data(dados.x_choudhury);
			estimator_choudhury.calc_lbub(S0, 0);
			t = clock();
			estimator_output results_choudhury = estimator_choudhury.run_estimator();
			t = clock() - t;
			time_taken = ((double)t) / CLOCKS_PER_SEC;
			write_estimation(save_dir, results_choudhury, &estimator_choudhury, time_taken);
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

			double stdNoise = 50 / pow(10.0, 25.0 / 10.0);

			std::vector<std::vector<double>> SP = exc_SP_cl_simulation(50, 0.25, 120, 1800, 1e-3);
			std::vector<std::vector<double>> exc = exc_vel_aleatoria(v_min * 10, v_max, S_exc, 300, 1e-3);

			//--------------------------------------------------------------------------------------------
			// Kano model
			ValveModel valve_kano = ValveModel();
			valve_kano.set_model(kano);
			valve_kano.set_valve_param_value(PARAM_VALVULA);
			valve_kano.set_friction_param_value(PARAM_ATRITO_KANO);
			valve_kano.set_simulation_type(h_cl);
			valve_kano.set_var_noise_controller(stdNoise);
			valve_kano.set_input_data(SP[1]);
			valve_kano.controller.set_controller_parameters({ 0.1, 0, 0, 0, 100, -100 });
			valve_kano.controller.set_excitation(exc[1]);
			valve_kano.controller.set_estimation(true);
			valve_kano.valve_simulation();
			// Cut data to start at the beginning of the excitation
			std::vector<double> P_kano_noise = simulateNoise(valve_kano.simulation_results.P, 25.0);
			procDataCL kano_data = preProcessCLdata(valve_kano.simulation_results.t, valve_kano.simulation_results.OP,
				P_kano_noise, valve_kano.simulation_results.x,
				valve_kano.controller.t_exc);
			// std::vector<std::vector<double>> save_data{ kano_data.t, kano_data.OP, kano_data.P, kano_data.x };
			// write_matrix(save_dir.append("test_kano.csv"), &save_data);

			// Add noise to the pressure data
			ValveModel valve_init_data = ValveModel();
			valve_init_data.set_model(kano);
			valve_init_data.set_valve_param_value(PARAM_VALVULA);
			// Filter data
			std::vector<double> P_kano_filt = valve_init_data.kalman_filter(&kano_data.OP, &kano_data.P, kano_data.Rv, kano_data.Rv / 50);
			kano_data.P = valve_init_data.filter2orderZP(&P_kano_filt, 1 / valve_init_data.get_tauip() * 100, 0.9);
			kano_data.x = valve_init_data.filter2orderZP(&kano_data.x, 1 / valve_init_data.get_tauip() * 100, 0.9);
			// Initialize the Kano model estimation
			Estimator estimator_kano = Estimator();
			estimator_kano.valve.set_model(kano);
			estimator_kano.valve.set_simulation_type(ol);
			estimator_kano.valve.set_valve_param_value(PARAM_VALVULA);
			estimator_kano.valve.set_input_data(kano_data.P);
			estimator_kano.valve.set_d0u0({ kano_data.d0, kano_data.u0 });
			estimator_kano.valve.set_pos0({ kano_data.x0, kano_data.d0 });
			estimator_kano.set_des_data(kano_data.x);
			estimator_kano.calc_lbub(S_exc, 0);
			clock_t t;
			t = clock();
			estimator_output results_kano = estimator_kano.run_estimator();
			t = clock() - t;
			double time_taken = ((double)t) / CLOCKS_PER_SEC;
			write_estimation(save_dir, results_kano, &estimator_kano, time_taken);
			

			//--------------------------------------------------------------------------------------------
			// He model
			ValveModel valve_he = ValveModel();
			valve_he.set_model(he);
			valve_he.set_valve_param_value(PARAM_VALVULA);
			valve_he.set_friction_param_value(PARAM_ATRITO_HE);
			valve_he.set_simulation_type(h_cl);
			valve_he.set_var_noise_controller(stdNoise);
			valve_he.set_input_data(SP[1]);
			valve_he.controller.set_controller_parameters({ 0.1, 0, 0, 0, 100, -100 });
			valve_he.controller.set_excitation(exc[1]);
			valve_he.controller.set_estimation(true);
			valve_he.valve_simulation();
			// Cut data to start at the beginning of the excitation
			std::vector<double> P_he_noise = simulateNoise(valve_he.simulation_results.P, 25.0);
			procDataCL he_data = preProcessCLdata(valve_he.simulation_results.t, valve_he.simulation_results.OP,
				P_he_noise, valve_he.simulation_results.x,
				valve_he.controller.t_exc);

			// Add noise to the pressure data
			valve_init_data.set_model(he);
			valve_init_data.set_valve_param_value(PARAM_VALVULA);
			// Filter data
			std::vector<double> P_he_filt = valve_init_data.kalman_filter(&he_data.OP, &he_data.P, he_data.Rv, he_data.Rv / 50);
			he_data.P = valve_init_data.filter2orderZP(&P_he_filt, 1 / valve_init_data.get_tauip() * 100, 0.9);
			he_data.x = valve_init_data.filter2orderZP(&he_data.x, 1 / valve_init_data.get_tauip() * 100, 0.9);
			// Initialize the He model estimation
			Estimator estimator_he = Estimator();
			estimator_he.valve.set_model(he);
			estimator_he.valve.set_simulation_type(ol);
			estimator_he.valve.set_valve_param_value(PARAM_VALVULA);
			estimator_he.valve.set_input_data(he_data.P);
			estimator_he.valve.set_d0u0({ he_data.d0, he_data.u0 });
			estimator_he.valve.set_pos0({ he_data.x0, he_data.d0 });
			estimator_he.set_des_data(he_data.x);
			estimator_he.calc_lbub(S_exc, 0);
			t = clock();
			estimator_output results_he = estimator_he.run_estimator();
			t = clock() - t;
			time_taken = ((double)t) / CLOCKS_PER_SEC;
			write_estimation(save_dir, results_he, &estimator_he, time_taken);


			//--------------------------------------------------------------------------------------------
			// Choudhury model
			ValveModel valve_choudhury = ValveModel();
			valve_choudhury.set_model(choudhury);
			valve_choudhury.set_valve_param_value(PARAM_VALVULA);
			valve_choudhury.set_friction_param_value(PARAM_ATRITO_CHOUDHURY);
			valve_choudhury.set_simulation_type(h_cl);
			valve_choudhury.set_var_noise_controller(stdNoise);
			valve_choudhury.set_input_data(SP[1]);
			valve_choudhury.controller.set_controller_parameters({ 0.1, 0, 0, 0, 100, -100 });
			valve_choudhury.controller.set_excitation(exc[1]);
			valve_choudhury.controller.set_estimation(true);
			valve_choudhury.valve_simulation();
			// Cut data to start at the beginning of the excitation
			std::vector<double> P_choudhury_noise = simulateNoise(valve_choudhury.simulation_results.P, 25.0);
			procDataCL choudhury_data = preProcessCLdata(valve_choudhury.simulation_results.t, valve_choudhury.simulation_results.OP,
				P_choudhury_noise, valve_choudhury.simulation_results.x,
				valve_choudhury.controller.t_exc);

			// Add noise to the pressure data
			valve_init_data.set_model(choudhury);
			valve_init_data.set_valve_param_value(PARAM_VALVULA);
			// Filter data
			std::vector<double> P_choudhury_filt = valve_init_data.kalman_filter(&choudhury_data.OP, &choudhury_data.P, choudhury_data.Rv, choudhury_data.Rv / 50);
			choudhury_data.P = valve_init_data.filter2orderZP(&P_choudhury_filt, 1 / valve_init_data.get_tauip() * 100, 0.9);
			choudhury_data.x = valve_init_data.filter2orderZP(&choudhury_data.x, 1 / valve_init_data.get_tauip() * 100, 0.9);
			// Initialize the Choudhury model estimation
			Estimator estimator_choudhury = Estimator();
			estimator_choudhury.valve.set_model(choudhury);
			estimator_choudhury.valve.set_simulation_type(ol);
			estimator_choudhury.valve.set_valve_param_value(PARAM_VALVULA);
			estimator_choudhury.valve.set_input_data(choudhury_data.P);
			estimator_choudhury.valve.set_d0u0({ choudhury_data.d0, choudhury_data.u0 });
			estimator_choudhury.valve.set_pos0({ choudhury_data.x0, choudhury_data.d0 });
			estimator_choudhury.set_des_data(choudhury_data.x);
			estimator_choudhury.calc_lbub(S_exc, 0);
			t = clock();
			estimator_output results_choudhury = estimator_choudhury.run_estimator();
			t = clock() - t;
			time_taken = ((double)t) / CLOCKS_PER_SEC;
			write_estimation(save_dir, results_choudhury, &estimator_choudhury, time_taken);


			////--------------------------------------------------------------------------------------------
			//// Karnopp model
			//ValveModel valve_karnopp = ValveModel();
			//valve_karnopp.set_model(karnopp);
			//valve_karnopp.set_valve_param_value(PARAM_VALVULA);
			//valve_karnopp.set_friction_param_value(PARAM_ATRITO_KARNOPP);
			//valve_karnopp.set_simulation_type(h_cl);
			//valve_karnopp.set_var_noise_controller(stdNoise);
			//valve_karnopp.set_input_data(SP[1]);
			//valve_karnopp.controller.set_controller_parameters({ 0.1, 0, 0, 0, 100, -100 });
			//valve_karnopp.controller.set_excitation(exc[1]);
			//valve_karnopp.controller.set_estimation(true);
			//valve_karnopp.valve_simulation();
			//// Cut data to start at the beginning of the excitation
			//std::vector<double> P_karnopp_noise = simulateNoise(valve_karnopp.simulation_results.P, 25.0);
			//procDataCL karnopp_data = preProcessCLdata(valve_karnopp.simulation_results.t, valve_karnopp.simulation_results.OP,
			//	P_karnopp_noise, valve_karnopp.simulation_results.x,
			//	valve_karnopp.controller.t_exc);

			//// Add noise to the pressure data
			//valve_init_data.set_model(karnopp);
			//valve_init_data.set_valve_param_value(PARAM_VALVULA);
			//// Filter data
			//std::vector<double> P_karnopp_filt = valve_init_data.kalman_filter(&karnopp_data.OP, &karnopp_data.P, karnopp_data.Rv, karnopp_data.Rv / 50);
			//karnopp_data.P = valve_init_data.filter2orderZP(&P_karnopp_filt, 1 / valve_init_data.get_tauip() * 100, 0.9);
			//karnopp_data.x = valve_init_data.filter2orderZP(&karnopp_data.x, 1 / valve_init_data.get_tauip() * 100, 0.9);
			//// Initialize the Karnopp model estimation
			//Estimator estimator_karnopp = Estimator();
			//estimator_karnopp.valve.set_model(karnopp);
			//estimator_karnopp.valve.set_simulation_type(ol);
			//estimator_karnopp.valve.set_valve_param_value(PARAM_VALVULA);
			//estimator_karnopp.valve.set_input_data(karnopp_data.P);
			//estimator_karnopp.valve.set_d0u0({ karnopp_data.d0, karnopp_data.u0 });
			//estimator_karnopp.valve.set_pos0({ karnopp_data.x0, karnopp_data.d0 });
			//estimator_karnopp.set_des_data(karnopp_data.x);
			//estimator_karnopp.calc_lbub(S_exc, 0);
			//t = clock();
			//estimator_output results_karnopp = estimator_karnopp.run_estimator();
			//t = clock() - t;
			//time_taken = ((double)t) / CLOCKS_PER_SEC;
			//write_estimation(save_dir, results_karnopp, &estimator_karnopp, time_taken);
			//

			////--------------------------------------------------------------------------------------------
			//// LuGre model
			//ValveModel valve_lugre = ValveModel();
			//valve_lugre.set_model(lugre);
			//valve_lugre.set_valve_param_value(PARAM_VALVULA);
			//valve_lugre.set_friction_param_value(PARAM_ATRITO_LUGRE);
			//valve_lugre.set_simulation_type(h_cl);
			//valve_lugre.set_var_noise_controller(stdNoise);
			//valve_lugre.set_input_data(SP[1]);
			//valve_lugre.controller.set_controller_parameters({ 0.1, 0, 0, 0, 100, -100 });
			//valve_lugre.controller.set_excitation(exc[1]);
			//valve_lugre.controller.set_estimation(true);
			//valve_lugre.valve_simulation();
			//// Cut data to start at the beginning of the excitation
			//std::vector<double> P_lugre_noise = simulateNoise(valve_lugre.simulation_results.P, 25.0);
			//procDataCL lugre_data = preProcessCLdata(valve_lugre.simulation_results.t, valve_lugre.simulation_results.OP,
			//	P_lugre_noise, valve_lugre.simulation_results.x,
			//	valve_lugre.controller.t_exc);

			//// Add noise to the pressure data
			//valve_init_data.set_model(lugre);
			//valve_init_data.set_valve_param_value(PARAM_VALVULA);
			//// Filter data
			//std::vector<double> P_lugre_filt = valve_init_data.kalman_filter(&lugre_data.OP, &lugre_data.P, lugre_data.Rv, lugre_data.Rv / 50);
			//lugre_data.P = valve_init_data.filter2orderZP(&P_lugre_filt, 1 / valve_init_data.get_tauip() * 100, 0.9);
			//lugre_data.x = valve_init_data.filter2orderZP(&lugre_data.x, 1 / valve_init_data.get_tauip() * 100, 0.9);
			//// Initialize the LuGre model estimation
			//Estimator estimator_lugre = Estimator();
			//estimator_lugre.valve.set_model(lugre);
			//estimator_lugre.valve.set_simulation_type(ol);
			//estimator_lugre.valve.set_valve_param_value(PARAM_VALVULA);
			//estimator_lugre.valve.set_input_data(lugre_data.P);
			//estimator_lugre.valve.set_d0u0({ lugre_data.d0, lugre_data.u0 });
			//estimator_lugre.valve.set_pos0({ lugre_data.x0, lugre_data.d0 });
			//estimator_lugre.set_des_data(lugre_data.x);
			//estimator_lugre.calc_lbub(S_exc, 0);
			//t = clock();
			//estimator_output results_lugre = estimator_lugre.run_estimator();
			//t = clock() - t;
			//time_taken = ((double)t) / CLOCKS_PER_SEC;
			//write_estimation(save_dir, results_lugre, &estimator_lugre, time_taken);


			////--------------------------------------------------------------------------------------------
			//// GMS model
			//ValveModel valve_gms = ValveModel();
			//valve_gms.set_model(gms);
			//valve_gms.set_valve_param_value(PARAM_VALVULA);
			//valve_gms.set_friction_param_value(PARAM_ATRITO_GMS);
			//valve_gms.set_simulation_type(h_cl);
			//valve_gms.set_var_noise_controller(stdNoise);
			//valve_gms.set_input_data(SP[1]);
			//valve_gms.controller.set_controller_parameters({ 0.1, 0, 0, 0, 100, -100 });
			//valve_gms.controller.set_excitation(exc[1]);
			//valve_gms.controller.set_estimation(true);
			//valve_gms.valve_simulation();
			//// Cut data to start at the beginning of the excitation
			//std::vector<double> P_gms_noise = simulateNoise(valve_gms.simulation_results.P, 25.0);
			//procDataCL gms_data = preProcessCLdata(valve_gms.simulation_results.t, valve_gms.simulation_results.OP,
			//	P_gms_noise, valve_gms.simulation_results.x,
			//	valve_gms.controller.t_exc);

			//// Add noise to the pressure data
			//valve_init_data.set_model(gms);
			//valve_init_data.set_valve_param_value(PARAM_VALVULA);
			//// Filter data
			//std::vector<double> P_gms_filt = valve_init_data.kalman_filter(&gms_data.OP, &gms_data.P, gms_data.Rv, gms_data.Rv / 50);
			//gms_data.P = valve_init_data.filter2orderZP(&P_gms_filt, 1 / valve_init_data.get_tauip() * 100, 0.9);
			//gms_data.x = valve_init_data.filter2orderZP(&gms_data.x, 1 / valve_init_data.get_tauip() * 100, 0.9);
			//// Initialize the GMS model estimation
			//Estimator estimator_gms = Estimator();
			//estimator_gms.valve.set_model(gms);
			//estimator_gms.valve.set_simulation_type(ol);
			//estimator_gms.valve.set_valve_param_value(PARAM_VALVULA);
			//estimator_gms.valve.set_input_data(gms_data.P);
			//estimator_gms.valve.set_d0u0({ gms_data.d0, gms_data.u0 });
			//estimator_gms.valve.set_pos0({ gms_data.x0, gms_data.d0 });
			//estimator_gms.set_des_data(gms_data.x);
			//estimator_gms.calc_lbub(S_exc, 0);
			//t = clock();
			//estimator_output results_gms = estimator_gms.run_estimator();
			//t = clock() - t;
			//time_taken = ((double)t) / CLOCKS_PER_SEC;
			//write_estimation(save_dir, results_gms, &estimator_gms, time_taken);
		}

		std::cout << "Estimation is over" << std::endl;
		std::cin;
	}
	//-------------------------------------------------------------------------------------------------------------
	// OL SIMULATION
	else if (type.compare("simulation") == 0) {

		std::string experimento = "sinusoidal_velocity"; // loaded_data, sinusoidal_velocity, aleatory_velocity

		// Valve data
		double w_n = 0.28;
		double S_exc = 25.7;
		double v_max = 3e-3 * 100 / 29e-3;
		double v_min = 1e-6 * 100 / 29e-3;

		std::vector<double> input;
		if (experimento.compare("loaded_data") == 0) {
			//std::string file("D:\\Doc2C\\test_data\\dados.csv");
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
		else {
			std::vector<std::vector<double>> exc = exc_vel_aleatoria(v_min*10, v_max, S_exc, 50, 1e-3);
			input = exc[1];
		}

		// Save the OP value
		 write_vector("D:\\Drive\\Projetos\\Doutorado2C\\test_data\\OP.csv", &input);

		// Initialize valve model
		ValveModel valve = ValveModel();
		std::vector<double> P;
		if (experimento.compare("loaded_data") == 0)
			valve.set_input_data(input);
		else {
			P = valve.OP2P_1order(&input, 1 / w_n, 1e-3);
			valve.set_input_data(P);
		}
		valve.set_valve_param_value(PARAM_VALVULA);

		//std::vector<double> P_filt = valve.kalman_filter(&input, &P, 1, 1/50);
		//std::vector<double> P_filt = valve.filter2orderZP(&P, w_n*100, 0.9);
		//std::vector<std::vector<double>> writer;
		//writer.push_back(P);
		//writer.push_back(P_filt);
		//write_matrix("D:\\Drive\\Projetos\\Doutorado2C\\test_data\\test_filter.csv", &writer);

		//// SIMULATE KANO MODEL
		//valve.set_model(kano);
		//std::vector<double> param_at PARAM_ATRITO_KANO;
		//valve.set_friction_param_value(param_at);
		//clock_t t;
		//t = clock();
		//valve.valve_simulation();
		//t = clock() - t;
		//double time_taken = ((double)t) / CLOCKS_PER_SEC;
		//printf("Kano exec time: %e\n", time_taken);
		//write_simulation("D:\\Drive\\Projetos\\Doutorado2C\\test_data\\simulation_kano.csv", valve.simulation_results);


		// SIMULATE HE MODEL
		valve.set_model(he);
		valve.set_friction_param_value(PARAM_ATRITO_HE);
		clock_t t;
		t = clock();
		valve.valve_simulation();
		t = clock() - t;
		double time_taken = ((double) t) / CLOCKS_PER_SEC;
		printf("He exec time: %e\n", time_taken);
		write_simulation("D:\\Drive\\Projetos\\Doutorado2C\\test_data\\simulation_he.csv", valve.simulation_results);


		// SIMULATE CHOUDHURY MODEL
		valve.set_model(choudhury);
		valve.set_friction_param_value(PARAM_ATRITO_CHOUDHURY);
		t = clock();
		valve.valve_simulation();
		t = clock() - t;
		time_taken = ((double)t) / CLOCKS_PER_SEC;
		printf("Choudhury exec time: %e\n", time_taken);
		write_simulation("D:\\Drive\\Projetos\\Doutorado2C\\test_data\\simulation_choudhury.csv", valve.simulation_results);


		//// SIMULATE KARNOPP MODEL
		//valve.set_model(karnopp);
		//param_at = PARAM_ATRITO_KARNOPP;
		//valve.set_friction_param_value(param_at);
		//t = clock();
		//valve.valve_simulation();
		//t = clock() - t;
		//time_taken = ((double)t) / CLOCKS_PER_SEC;
		//printf("Karnopp exec time: %e\n", time_taken);
		//write_simulation("D:\\Drive\\Projetos\\Doutorado2C\\test_data\\simulation_karnopp.csv", valve.simulation_results);


		//// SIMULATE LUGRE MODEL
		//valve.set_model(lugre);
		//param_at = PARAM_ATRITO_LUGRE;
		//valve.set_friction_param_value(param_at);
		//t = clock();
		//valve.valve_simulation();
		//t = clock() - t;
		//time_taken = ((double)t) / CLOCKS_PER_SEC;
		//printf("LuGre exec time: %e\n", time_taken);
		//write_simulation("D:\\Drive\\Projetos\\Doutorado2C\\test_data\\simulation_lugre.csv", valve.simulation_results);


		//// SIMULATE GMS MODEL
		//valve.set_model(gms);
		//param_at = PARAM_ATRITO_GMS;
		//valve.set_friction_param_value(param_at);
		//t = clock();
		//valve.valve_simulation();
		//t = clock() - t;
		//time_taken = ((double)t) / CLOCKS_PER_SEC;
		//printf("GMS exec time: %e\n", time_taken);
		//write_simulation("D:\\Drive\\Projetos\\Doutorado2C\\test_data\\simulation_gms.csv", valve.simulation_results);

		std::cout << "End of simulations" << std::endl;
		std::cin;
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
		std::vector<std::vector<double>> exc = exc_vel_aleatoria(v_min * 10, v_max, S_exc, 300, 1e-3);

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


		// Choudhury model cl simulation
		valve.set_model(choudhury);
		valve.set_valve_param_value(PARAM_VALVULA);
		valve.set_friction_param_value(PARAM_ATRITO_CHOUDHURY);
		valve.set_simulation_type(h_cl);
		valve.set_var_noise_controller(stdNoise);
		valve.set_input_data(SP[1]);
		valve.controller.set_controller_parameters({ 0.1, 0, 0, 0, 100, -100 });
		valve.controller.set_excitation(exc[1]);
		valve.controller.set_estimation(true);
		valve.valve_simulation();
		write_simulation("D:\\Drive\\Projetos\\Doutorado2C\\test_data\\simulation_choudhury_cl.csv", valve.simulation_results);


		//// Karnopp model cl simulation
		//valve.set_model(karnopp);
		//valve.set_valve_param_value(PARAM_VALVULA);
		//valve.set_friction_param_value(PARAM_ATRITO_KARNOPP);
		//valve.set_simulation_type(h_cl);
		//valve.set_var_noise_controller(stdNoise);
		//valve.set_input_data(SP[1]);
		//valve.controller.set_controller_parameters({ 0.1, 0, 0, 0, 100, -100 });
		//valve.controller.set_excitation(exc[1]);
		//valve.controller.set_estimation(true);
		//valve.valve_simulation();
		//write_simulation("D:\\Drive\\Projetos\\Doutorado2C\\test_data\\simulation_karnopp_cl.csv", valve.simulation_results);


		//// LuGre model cl simulation
		//valve.set_model(lugre);
		//valve.set_valve_param_value(PARAM_VALVULA);
		//valve.set_friction_param_value(PARAM_ATRITO_LUGRE);
		//valve.set_simulation_type(h_cl);
		//valve.set_var_noise_controller(stdNoise);
		//valve.set_input_data(SP[1]);
		//valve.controller.set_controller_parameters({ 0.1, 0, 0, 0, 100, -100 });
		//valve.controller.set_excitation(exc[1]);
		//valve.controller.set_estimation(true);
		//valve.valve_simulation();
		//write_simulation("D:\\Drive\\Projetos\\Doutorado2C\\test_data\\simulation_lugre_cl.csv", valve.simulation_results);


		//// GMS model cl simulation
		//valve.set_model(gms);
		//valve.set_valve_param_value(PARAM_VALVULA);
		//valve.set_friction_param_value(PARAM_ATRITO_GMS);
		//valve.set_simulation_type(h_cl);
		//valve.set_var_noise_controller(stdNoise);
		//valve.set_input_data(SP[1]);
		//valve.controller.set_controller_parameters({ 0.1, 0, 0, 0, 100, -100 });
		//valve.controller.set_excitation(exc[1]);
		//valve.controller.set_estimation(true);
		//valve.valve_simulation();
		//write_simulation("D:\\Drive\\Projetos\\Doutorado2C\\test_data\\simulation_gms_cl.csv", valve.simulation_results);
	}
	return 0;
}