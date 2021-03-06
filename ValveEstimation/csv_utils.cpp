

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include "csv_utils.h"
#include "valve_models.h"
#include "estimator.h"
#include <ctime>
#include <algorithm>


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


csv_real_data get_real_data(const std::string filename) {

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

	csv_real_data loaded_data;
	for (auto row : data) {
		if (row.size() == 3) {
			loaded_data.OP.push_back(stod(row[0]));
			loaded_data.P.push_back(stod(row[1]));
			loaded_data.x.push_back(stod(row[2]));
		}
	}

	fin.close();

	return loaded_data;
}



void writer_variable(std::ofstream& file, const std::vector<double>& variable) {
	for (double data : variable)
		file.write((char*)&data, sizeof(double));
}

void write_header_sim(std::ofstream& file, unsigned int n_variables, unsigned int size_vector, std::vector<std::string> variables) {
	// file identifier (.sim 0x73, 0x69, 0x6D, 0x00, 0x00)
	char* identification = new char[5]{ 0x73, 0x69, 0x6D, 0x00, 0x00 };
	file.write((&(*identification)), sizeof(identification[0]) * 5);

	file.write((char*)&n_variables, sizeof(unsigned int));
	file.write((char*)&size_vector, sizeof(size_vector));

	for (std::string variable : variables) {
		char tag[20]{};
		for (int i = 0; i < variable.size() && i < 20; i++)
			tag[i] = variable.c_str()[i];
		file.write((&(*tag)), sizeof(tag[0]) * 20);
	}
}

void write_simulation(std::string filename, const simdata& simulacao) {
	
#ifdef _WRITE_SIMULATION_CSV
	
	std::ofstream fout(filename);
	fout << std::fixed << std::setprecision(20);
	for (int i = 0; i < simulacao.t.size(); ++i) {
		fout << simulacao.t[i] << "," << simulacao.OP[i] << "," << simulacao.P[i] << "," << simulacao.x[i] << "," << simulacao.v[i] <<
			"," << simulacao.a[i] << "," << simulacao.F_at[i] << "," << simulacao.SP[i] << "\n";
	}
	fout.close()

#else

	std::ofstream fout(filename, std::ios::out | std::ios::binary);
	
	std::vector<std::string> variables{ "t", "OP", "P", "x", "v", "a", "F_at", "SP" };
	write_header_sim(fout, 8, simulacao.t.size(), variables);
	writer_variable(fout, simulacao.t);
	writer_variable(fout, simulacao.OP);
	writer_variable(fout, simulacao.P);
	writer_variable(fout, simulacao.x);
	writer_variable(fout, simulacao.v);
	writer_variable(fout, simulacao.a);
	writer_variable(fout, simulacao.F_at);
	writer_variable(fout, simulacao.SP);

	fout.close();

#endif
}




void write_estimation(std::string dir, estimator_output data, std::string model, double run_time, double min_residual) {

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
		"time:," << run_time << "\n" <<
		"original param. residual:," << min_residual << "\n";
	for (int i = 0; i < data.parameters[0].size(); i++) {
		fout << data.parameters[0][i];
		if (i < data.parameters[0].size() - 1)
			fout << ",";
		else
			fout << "\n";
	}
	fout.close();
}


void write_matrix(std::string filename, const std::vector<std::vector<double>>& u, std::vector<std::string> variables) {

#ifdef _WRITE_SIMULATION_CSV

	std::ofstream fout(filename);
	fout << std::fixed << std::setprecision(20);
	size_t minsize = { 1000000000000 };
	for (int i = 0; i < (*u).size(); ++i) {
		minsize = std::min((*u)[i].size(), minsize);
	}
	for (int i = 0; i < minsize; ++i) {
		for (int j = 0; j < (*u).size(); ++j) {
			fout << (*u)[j][i];
			if (j == (*u).size() - 1)
				fout << "\n";
			else
				fout << ",";
}
	}
	fout.close();

#else

	std::ofstream fout(filename, std::ios::out | std::ios::binary);

	write_header_sim(fout, u.size(), u[0].size(), variables);

	for (auto vector : u)
		writer_variable(fout, vector);

	fout.close();

#endif
}