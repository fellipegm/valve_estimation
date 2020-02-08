

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include "csv_utils.h"
#include "valve_models.h"
#include "estimator.h"
#include <ctime>


csv_data get_data(const std::string filename) {

	std::fstream fin;
	fin.open(filename, std::ios::in);

	std::string line, word;
	std::vector<std::vector<std::string>> data;
	std::vector<std::string> row;
	while (fin >> line) {
		row.clear();
		std::stringstream ss(line);
		while (getline(ss, word, ','))
			row.push_back(word);
		data.push_back(row);
	}

	csv_data loaded_data;
	for (auto row : data) {
		if (row.size() == 7) {
			loaded_data.t.push_back(stod(row[0]));
			loaded_data.OP.push_back(stod(row[1]));
			loaded_data.P.push_back(stod(row[2]));
			loaded_data.x_kano.push_back(stod(row[3]));
			loaded_data.x_karnopp.push_back(stod(row[4]));
			loaded_data.x_lugre.push_back(stod(row[5]));
			loaded_data.x_gms.push_back(stod(row[6]));
		}
	}

	fin.close();

	return loaded_data;
}

void write_simulation(std::string filename, simdata simulacao) {
	std::ofstream fout(filename);
	
	fout << std::fixed << std::setprecision(20);
	for (int i = 0; i < simulacao.t.size(); ++i) {
		fout << simulacao.t[i] << "," << simulacao.P[i] << "," << simulacao.x[i] << "," << simulacao.v[i] << "," << simulacao.a[i] << "," <<
			simulacao.F_at[i] << "\n";
	}
	
	fout.close();
}


void write_estimation(std::string dir, estimator_output data, Estimator* estimator, double run_time) {
	// Get valve model
	std::string model = "";
	switch (estimator->valve.get_model()){
		case kano:
			model = "kano";
			break;
		case karnopp:
			model = "karnopp";
			break;
		case lugre:
			model = "lugre";
			break;
		case gms:
			model = "gms";
			break;
		default:
			break;
	}
	std::time_t t = std::time(0);
	std::tm* now = std::localtime(&t);
	std::stringstream filename;
	filename << std::to_string(now->tm_mday) << "-" << std::to_string(now->tm_mon + 1) << "-" << std::to_string(now->tm_year + 1900) << '_' <<
		now->tm_hour << '-' << now->tm_min << '-' << now->tm_sec << '_' << 
		model;

	std::ofstream fout(dir.append(filename.str()).append(".out"));
	fout << std::fixed << std::setprecision(20);

	fout << "N Func. Evals:," << data.nFuncEvals << "\n" <<
		"residual:," << data.residuals[0] << "\n" <<
		"time:," << run_time << "\n";
	for (int i = 0; i < data.parameters[0].size(); i++) {
		fout << data.parameters[0][i];
		if (i < data.parameters[0].size() - 1)
			fout << ",";
		else
			fout << "\n";
	}
	fout.close();
}


void write_vector(std::string filename, std::vector<double>* u) {
	std::ofstream fout(filename);
	fout << std::fixed << std::setprecision(20);
	for (auto data : *u) {
		fout << data << "\n";
	}
	fout.close();
}


void write_matrix(std::string filename, std::vector<std::vector<double>>* u) {
	std::ofstream fout(filename);
	fout << std::fixed << std::setprecision(20);
	for (int i = 0; i < (*u)[0].size(); ++i) {
		for (int j = 0; j < (*u).size(); ++j) {
			fout << (*u)[j][i];
			if (j == (*u).size() - 1)
				fout << "\n";
			else
				fout << ",";
		}
	}
	fout.close();
}