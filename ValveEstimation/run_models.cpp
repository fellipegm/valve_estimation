#include <iostream>
#include <string>
#include <cstring>
#include <ctime>
#include <vector>
#include <direct.h>

#include "valve_models.h"
#include "csv_utils.h"
#include "utils.h"


#define PARAM_VALVULA {1.6, 210490, 0.0445, 2550, 0, 0.029, 180000, 41368, 3.571428571428571}
#define PARAM_ATRITO_KANO {25.287220681574890, 1.296780547773071, 11.494727614487891, 0.0, 0.029}
#define PARAM_ATRITO_KARNOPP {700, 780, 125000, 5.0e-04}
#define PARAM_ATRITO_LUGRE {700, 780, 125000, 5.0e-04, 26000000, 2.039607805437114e+04}
#define PARAM_ATRITO_GMS {700, 780, 125000, 5.0e-04, 29250000, 15600000, 31200000, 2910.90841109400, 1455.45420554700, 4366.36261664100, 0.75, 0.2, 20}

#define TS 1e-3

int main() {

	std::string experimento = "aleatory_velocity"; // loaded_data, sinusoidal_velocity, aleatory_velocities

	// Valve data
	double w_n = 0.28;
	double S_exc = 25.7;
	double v_max = 3e-3 * 100 / 29e-3;
	double v_min = 1e-5 * 100 / 29e-3;

	std::vector<double> input;
	if (experimento.compare("loaded_data") == 0) {
		std::string file("D:\\Doc2C\\test_data\\dados.csv");

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
		std::vector<std::vector<double>> exc = exc_vel_aleatoria(v_min, v_max, S_exc, 50, 1e-3);
		input = exc[1];
	}
	

	// Initialize valve model
	ValveModel* valve;
	valve = ValveModel::New();
	if (experimento.compare("loaded_data") == 0)
		valve->set_input_data(input);
	else {
		std::vector<double> P = valve->OP2P_1order(&input, 1 / w_n, 1e-3);
		valve->set_input_data(P);
	}
	std::vector<double> param_valv PARAM_VALVULA;
	valve->set_valve_param_value(param_valv);
	
	// SIMULATE KANO MODEL
	valve->set_model(kano);
	std::vector<double> param_at PARAM_ATRITO_KANO;
	valve->set_friction_param_value(param_at);
	clock_t t;
	t = clock();
	valve->valve_simulation();
	t = clock() - t;
	double time_taken = ((double)t) / CLOCKS_PER_SEC;
	printf("Kano exec time: %e\n", time_taken);
	write_simulation("D:\\Doc2C\\test_data\\simulation_kano.csv", valve->simulation_results);


	// SIMULATE KARNOPP MODEL
	valve->set_model(karnopp);
	param_at = PARAM_ATRITO_KARNOPP;
	valve->set_friction_param_value(param_at);
	t = clock();
	valve->valve_simulation();
	t = clock() - t;
	time_taken = ((double)t) / CLOCKS_PER_SEC;
	printf("Karnopp exec time: %e\n", time_taken);
	write_simulation("D:\\Doc2C\\test_data\\simulation_karnopp.csv", valve->simulation_results);


	// SIMULATE LUGRE MODEL
	valve->set_model(lugre);
	param_at = PARAM_ATRITO_LUGRE;
	valve->set_friction_param_value(param_at);
	t = clock();
	valve->valve_simulation();
	t = clock() - t;
	time_taken = ((double)t) / CLOCKS_PER_SEC;
	printf("LuGre exec time: %e\n", time_taken);
	write_simulation("D:\\Doc2C\\test_data\\simulation_lugre.csv", valve->simulation_results);
	

	// SIMULATE GMS MODEL
	valve->set_model(gms);
	param_at = PARAM_ATRITO_GMS;
	valve->set_friction_param_value(param_at);
	t = clock();
	valve->valve_simulation();
	t = clock() - t;
	time_taken = ((double)t) / CLOCKS_PER_SEC;
	printf("GMS exec time: %e\n", time_taken);
	write_simulation("D:\\Doc2C\\test_data\\simulation_gms.csv", valve->simulation_results);

	delete valve;

	std::cout << "End of simulations" << std::endl;
	std::cin;

	return 0;
}