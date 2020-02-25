#pragma once


#ifndef _UTILS_H
#define _UTILS_H

#include <vector>
#include <algorithm>
#include <numeric>


# define M_PI           3.14159265358979323846


typedef struct procDataCL
{
	std::vector<double> t;
	std::vector<double> OP;
	std::vector<double> P;
	std::vector<double> x;
	double x0{ 0.0 };
	double u0{ 0.0 };
	double d0{ 1 };
	double Rv{ 0.0 };
} procDataCL;


std::vector<size_t> sort_indexes(const std::vector<double>& v);

std::vector<double> sum_vectors(const std::vector<double>&, const std::vector<double>&);

std::vector<double> subtract_vectors(const std::vector<double>&, const std::vector<double>&);
std::vector<double> subtract_vect_const(const std::vector<double>& v1, const double c2);
double min_vec(const std::vector<double>& v1);
double max_vec(const std::vector<double>& v1);
std::vector<double> abs_vec(const std::vector<double>& v1);

std::vector<double> vector_const_multiply(const std::vector<double>&, const double);

std::vector<std::vector<double>> sum_matrices(const std::vector<std::vector<double>>&, const std::vector<std::vector<double>>&);

std::vector<std::vector<double>> subtract_matrices(const std::vector<std::vector<double>>&, const std::vector<std::vector<double>>&);

std::vector<std::vector<double>> matrix_const_multiply(const std::vector<std::vector<double>>&, const double);

std::vector<std::vector<double>> matrix_element_multiply(const std::vector<std::vector<double>>&, const std::vector<std::vector<double>>&);



template<typename T>
class Logspace {
private:
	T curValue, base;

public:
	Logspace(T first, T base) : curValue(first), base(base) {}

	T operator()() {
		T retval = curValue;
		curValue *= base;
		return retval;
	}
};

std::vector<double> pyLogspace(double, double, int, double);
std::vector<double> pyLinspace(double start, double stop, int num);

std::vector<std::vector<double>> exc_vel_senoidal(double, double, double, double, double, double, int, double);
std::vector<std::vector<double>> exc_vel_aleatoria(double v_min, double v_max, double S, int n, double dt);
std::vector<std::vector<double>> exc_SP_cl_simulation(double SPm, double varSP, double Tsp, double Tsim, double dt);

procDataCL preProcessCLdata(const std::vector<double>& t, const std::vector<double>& OP,
	const std::vector<double>& P, const std::vector<double>& x, double t_exc);

std::vector<double> filter1order(std::vector<double>* data, double tau, double Ts);


template<typename T>
T variance_vec(const std::vector<T>& vec)
{
	size_t sz = vec.size();
	if (sz == 1)
		return 0.0;

	// Calculate the mean
	T mean = std::accumulate(vec.begin(), vec.end(), 0.0) / sz;

	// Now calculate the variance
	auto variance_func = [&mean, &sz](T accumulator, const T& val)
	{
		return accumulator + ((val - mean) * (val - mean) / (sz - 1));
	};

	return std::accumulate(vec.begin(), vec.end(), 0.0, variance_func);
}

template<typename T>
T std_vec(const std::vector<T>& vec)
{
	return sqrt(variance_vec(vec));
}

template<typename T>
T mean_vec(const std::vector<T>& vec)
{
	return std::accumulate(vec.begin(), vec.end(), 0.0)/std::size(vec);
}

std::vector<double> simulateNoise(const std::vector<double>& data, double snr);
std::vector<double> simulateNoise(const std::vector<double>& data, double snr, double mean);

template<typename T>
T signal_fnc(const T value)
{
	if (value > 0)
		return 1;
	else if (value < 0)
		return -1;
	else
		return 0;
}

#endif