#pragma once


#ifndef _UTILS_H
#define _UTILS_H

#include <vector>
#include <algorithm>
#include <numeric>


# define M_PI           3.14159265358979323846


std::vector<size_t> sort_indexes(const std::vector<double>& v);

std::vector<double> sum_vectors(const std::vector<double>&, const std::vector<double>&);

std::vector<double> subtract_vectors(const std::vector<double>&, const std::vector<double>&);

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
T mean_vec(const std::vector<T>& vec)
{
	return std::accumulate(vec.begin(), vec.end(), 0.0)/std::size(vec);
}

std::vector<double> simulateNoise(const std::vector<double>& data, double snr);

#endif