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

	
	std::string type = "estimation"; // estimation or simulation


	if (type.compare("estimation") == 0) {
		//std::string file("D:\\Doc2C\\test_data\\data_estimation.csv");
		std::string estimation_data("D:\\Drive\\Projetos\\Doutorado2C\\test_data\\data_estimation_sin.csv");
		std::string save_dir("D:\\Drive\\Projetos\\Doutorado2C\\test_data\\");

		csv_data dados;

		double S0 = 25.7;
		int n_tests = 10;
		bool simulate_noise = true;

		for (int i = 0; i < n_tests; ++i) {
			if (simulate_noise) {
				dados = get_data(estimation_data);
				std::vector<double> P_noise = simulateNoise(dados.P, 25);
				std::vector<double> x_kano_noise = simulateNoise(dados.x_kano, 25);
				std::vector<double> x_karnopp_noise = simulateNoise(dados.x_karnopp, 25);
				std::vector<double> x_lugre_noise = simulateNoise(dados.x_lugre, 25);
				std::vector<double> x_gms_noise = simulateNoise(dados.x_gms, 25);

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
			}

			// Initialize the Kano model estimation
			Estimator estimator_kano = Estimator();
			estimator_kano.valve.set_model(kano);
			estimator_kano.valve.set_valve_param_value(PARAM_VALVULA);
			estimator_kano.valve.set_input_data(dados.P);
			estimator_kano.set_des_data(dados.x_kano);
			estimator_kano.calc_lbub(S0, 0);
			clock_t t;
			t = clock();
			estimator_output results_kano = estimator_kano.run_estimator();
			t = clock() - t;
			double time_taken = ((double)t) / CLOCKS_PER_SEC;
			write_estimation(save_dir, results_kano, &estimator_kano, time_taken);


			// Initialize the Karnopp model estimation
			Estimator estimator_karnopp = Estimator();
			estimator_karnopp.valve.set_model(karnopp);
			estimator_karnopp.valve.set_valve_param_value(PARAM_VALVULA);
			estimator_karnopp.valve.set_input_data(dados.P);
			estimator_karnopp.set_des_data(dados.x_karnopp);
			estimator_karnopp.calc_lbub(S0, 0);
			t = clock();
			estimator_output results_karnopp = estimator_karnopp.run_estimator();
			t = clock() - t;
			time_taken = ((double)t) / CLOCKS_PER_SEC;
			write_estimation(save_dir, results_karnopp, &estimator_karnopp, time_taken);


			// Initialize the LuGre model estimation
			Estimator estimator_lugre = Estimator();
			estimator_lugre.valve.set_model(lugre);
			estimator_lugre.valve.set_valve_param_value(PARAM_VALVULA);
			estimator_lugre.valve.set_input_data(dados.P);
			estimator_lugre.set_des_data(dados.x_lugre);
			estimator_lugre.calc_lbub(S0, 0);
			t = clock();
			estimator_output results_lugre = estimator_lugre.run_estimator();
			t = clock() - t;
			time_taken = ((double)t) / CLOCKS_PER_SEC;
			write_estimation(save_dir, results_lugre, &estimator_lugre, time_taken);


			// Initialize the GMS model estimation
			Estimator estimator_gms = Estimator();
			estimator_gms.valve.set_model(gms);
			estimator_gms.valve.set_valve_param_value(PARAM_VALVULA);
			estimator_gms.valve.set_input_data(dados.P);
			estimator_gms.set_des_data(dados.x_gms);
			estimator_gms.calc_lbub(S0, 0);
			t = clock();
			estimator_output results_gms = estimator_gms.run_estimator();
			t = clock() - t;
			time_taken = ((double)t) / CLOCKS_PER_SEC;
			write_estimation(save_dir, results_gms, &estimator_gms, time_taken);
		}

		std::cout << "Estimation is over" << std::endl;
		std::cin;
	}
	//-------------------------------------------------------------------------------------------------------------
	// SIMULATION
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

		// SIMULATE KANO MODEL
		valve.set_model(kano);
		std::vector<double> param_at PARAM_ATRITO_KANO;
		valve.set_friction_param_value(param_at);
		clock_t t;
		t = clock();
		valve.valve_simulation();
		t = clock() - t;
		double time_taken = ((double)t) / CLOCKS_PER_SEC;
		printf("Kano exec time: %e\n", time_taken);
		write_simulation("D:\\Drive\\Projetos\\Doutorado2C\\test_data\\simulation_kano.csv", valve.simulation_results);


		// SIMULATE KARNOPP MODEL
		valve.set_model(karnopp);
		param_at = PARAM_ATRITO_KARNOPP;
		valve.set_friction_param_value(param_at);
		t = clock();
		valve.valve_simulation();
		t = clock() - t;
		time_taken = ((double)t) / CLOCKS_PER_SEC;
		printf("Karnopp exec time: %e\n", time_taken);
		write_simulation("D:\\Drive\\Projetos\\Doutorado2C\\test_data\\simulation_karnopp.csv", valve.simulation_results);


		// SIMULATE LUGRE MODEL
		valve.set_model(lugre);
		param_at = PARAM_ATRITO_LUGRE;
		valve.set_friction_param_value(param_at);
		t = clock();
		valve.valve_simulation();
		t = clock() - t;
		time_taken = ((double)t) / CLOCKS_PER_SEC;
		printf("LuGre exec time: %e\n", time_taken);
		write_simulation("D:\\Drive\\Projetos\\Doutorado2C\\test_data\\simulation_lugre.csv", valve.simulation_results);


		// SIMULATE GMS MODEL
		valve.set_model(gms);
		param_at = PARAM_ATRITO_GMS;
		valve.set_friction_param_value(param_at);
		t = clock();
		valve.valve_simulation();
		t = clock() - t;
		time_taken = ((double)t) / CLOCKS_PER_SEC;
		printf("GMS exec time: %e\n", time_taken);
		write_simulation("D:\\Drive\\Projetos\\Doutorado2C\\test_data\\simulation_gms.csv", valve.simulation_results);

		std::cout << "End of simulations" << std::endl;
		std::cin;
	}

	return 0;
}