#pragma once

#ifndef _OPEN_CSV_H_
#define _OPEN_CSV_H_

#include <vector>
#include <iostream>
#include "valve_models.h"
#include "estimator.h"

typedef struct csv_data
{
	std::vector<double> t;
	std::vector<double> P;
	std::vector<double> OP;
	std::vector<double> x_kano;
	std::vector<double> x_karnopp;
	std::vector<double> x_lugre;
	std::vector<double> x_gms;
	std::vector<double> x_he;
	std::vector<double> x_choudhury;
} csv_data;

typedef struct csv_real_data
{
	std::vector<double> OP;
	std::vector<double> P;
	std::vector<double> x;
} csv_real_data;

csv_data get_data(const std::string);
csv_real_data get_real_data(const std::string filename);
void write_simulation(std::string, simdata);
void write_estimation(std::string dir, estimator_output data, std::string model, double run_time, double min_residual);
void write_vector(std::string filename, std::vector<double>* u);
void write_matrix(std::string filename, std::vector<std::vector<double>>* u);

#endif // _OPEN_CSV_H_