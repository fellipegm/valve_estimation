

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
		if (row.size() == 9) {
			loaded_data.t.push_back(stod(row[0]));
			loaded_data.OP.push_back(stod(row[1]));
			loaded_data.P.push_back(stod(row[2]));
			loaded_data.x_kano.push_back(stod(row[3]));
			loaded_data.x_karnopp.push_back(stod(row[4]));
			loaded_data.x_lugre.push_back(stod(row[5]));
			loaded_data.x_gms.push_back(stod(row[6]));
			loaded_data.x_he.push_back(stod(row[7]));
			loaded_data.x_choudhury.push_back(stod(row[8]));
		}
	}

	fin.close();

	return loaded_data;
}

void write_simulation(std::string filename, simdata simulacao) {
	std::ofstream fout(filename);
	
	fout << std::fixed << std::setprecision(20);
	for (int i = 0; i < simulacao.t.size(); ++i) {
		fout << simulacao.t[i] << "," << simulacao.OP[i] << "," << simulacao.P[i] << "," << simulacao.x[i] << "," << simulacao.v[i] << 
			"," << simulacao.a[i] << "," << simulacao.F_at[i] << "," << simulacao.SP[i] << "\n";
	}
	
	fout.close();
}


void write_estimation(std::string dir, estimator_output data, std::string model, double run_time) {

	std::time_t now = std::time(0);
	struct std::tm newtime; 
	localtime_s(&newtime, &now);
	std::stringstream filename;
	filename << std::to_string(newtime.tm_mday) << "-" << std::to_string(newtime.tm_mon + 1) << "-" << std::to_string(newtime.tm_year + 1900) << '_' <<
		newtime.tm_hour << '-' << newtime.tm_min << '-' << newtime.tm_sec << '_' <<
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