

#include "utils.h"
#include <algorithm>
#include <numeric>
#include <math.h>
#include <random>
#include <chrono>
#include <float.h>


std::vector<size_t> sort_indexes(const std::vector<double>& v) {
	// initialize original index locations
	std::vector<size_t> idx(v.size());
	std::iota(idx.begin(), idx.end(), 0);

	// sort indexes based on comparing values in v
	sort(idx.begin(), idx.end(),
		[&v](size_t i1, size_t i2) {return v[i1] < v[i2]; });

	return idx;
}


std::vector<double> sum_vectors(const std::vector<double>& v1, const std::vector<double>& v2) {
	std::vector<double> results(v1.size());
	for (int i = 0; i < v1.size(); i++)
		results[i] = v1[i] + v2[i];
	return results;
}

std::vector<double> subtract_vectors(const std::vector<double>& v1, const std::vector<double>& v2) {
	std::vector<double> results(v1.size());
	for (int i = 0; i < v1.size(); i++)
		results[i] = v1[i] - v2[i];
	return results;
}

std::vector<double> subtract_vect_const(const std::vector<double>& v1, const double c2) {
	std::vector<double> results(v1.size());
	for (int i = 0; i < v1.size(); i++)
		results[i] = v1[i] - c2;
	return results;
}

double min_vec(const std::vector<double>& v1) {
	double result = DBL_MAX;
	for (int i = 0; i < v1.size(); i++)
		result = std::min(result, v1[i]);
	return result;
}

double max_vec(const std::vector<double>& v1) {
	double result = -DBL_MAX;
	for (int i = 0; i < v1.size(); i++)
		result = std::max(result, v1[i]);
	return result;
}

std::vector<double> abs_vec(const std::vector<double>& v1) {
	std::vector<double> results(v1.size());
	for (int i = 0; i < v1.size(); i++)
		results[i] = std::abs(v1[i]);
	return results;
}


std::vector<double> vector_const_multiply(const std::vector<double>& v1, const double constant) {
	std::vector<double> results(v1.size());
	for (int i = 0; i < v1.size(); i++)
		results[i] = v1[i] * constant;
	return results;
}


std::vector<std::vector<double>> sum_matrices(const std::vector<std::vector<double>>& m1, const std::vector<std::vector<double>>& m2) {
	std::vector<std::vector<double>> results(m1);
	for (int i = 0; i < m1.size(); i++)
		for (int j = 0; j < m1[0].size(); j++)
			results[i][j] = m1[i][j] + m2[i][j];
	return results;
}

std::vector<std::vector<double>> subtract_matrices(const std::vector<std::vector<double>>& m1, const std::vector<std::vector<double>>& m2) {
	std::vector<std::vector<double>> results(m1);
	for (int i = 0; i < m1.size(); i++)
		for (int j = 0; j < m1[0].size(); j++)
			results[i][j] = m1[i][j] - m2[i][j];
	return results;
}

std::vector<std::vector<double>> matrix_const_multiply(const std::vector<std::vector<double>>& m1, const double constant) {
	std::vector<std::vector<double>> results(m1);
	for (int i = 0; i < m1.size(); i++)
		for (int j = 0; j < m1[0].size(); j++)
			results[i][j] = m1[i][j] * constant;
	return results;
}


std::vector<std::vector<double>> matrix_element_multiply(const std::vector<std::vector<double>>& m1, const std::vector<std::vector<double>>& m2) {
	std::vector<std::vector<double>> results(m1);
	for (int i = 0; i < m1.size(); i++)
		for (int j = 0; j < m1[0].size(); j++)
			results[i][j] = m1[i][j] * m2[i][j];
	return results;
}


std::vector<std::vector<double>> exc_vel_senoidal(double v_min, double v_max, double T_min, double T_max, double S, double tau, int n, double dt) {
	double exc_max = 100;

	// Periods and velocities
	std::vector<double> T, v;
	T = pyLogspace(log10(T_min), log10(T_max), n, 10);
	v = pyLogspace(log10(v_min*10), log10(v_max), n*2, 10);

	double t_init;
	std::vector<double> v_atual, t, exc;
	t.push_back(0);
	exc.push_back(0);
	for (int i = 0; i < T.size(); ++i) {
		v_atual.push_back(0.0);
		while (exc.back() < 0.9 * S) {
			exc.push_back(exc.back() + 2*dt);
			t.push_back(t.back() + dt);
		}
		while (exc.back() < 1.1 * S) {
			exc.push_back(exc.back() + v[i*2] * dt);
			t.push_back(t.back() + dt);
		}
		for (int ct = 0; ct < int(5 * tau / dt); ++ct) {
			exc.push_back(1.1*S);
			t.push_back(t.back() + dt);
		}

		t_init = t.back();
		while (exc.back() < exc_max) {
			t.push_back(t.back() + dt);
			v_atual.push_back( ((v_max - v_min) / 2 + v_min - (v_max - v_min) / 2 * cos(2 * M_PI / T[i] * (t.back() - t_init))) );
			exc.push_back(exc.back() + v_atual.back() * dt);
		}

		for (int ct = 0; ct < int(5 * tau / dt); ++ct) {
			exc.push_back(exc_max);
			t.push_back(t.back() + dt);
		}
		while (exc.back() > 100-0.9*S) {
			exc.push_back(exc.back() - 2 * dt);
			t.push_back(t.back() + dt);
		}
		while (exc.back() >100-1.1*S) {
			exc.push_back(exc.back() - v[i*2+1] * dt);
			t.push_back(t.back() + dt);
		}
		while (exc.back() > 0) {
			exc.push_back(exc.back() -10 * dt);
			t.push_back(t.back() + dt);
		}
		for (int ct = 0; ct < int(5 * tau / dt); ++ct) {
			exc.push_back(0);
			t.push_back(t.back() + dt);
		}
	}

	std::vector<std::vector<double>> exc_signal;
	exc_signal.push_back(t);
	exc_signal.push_back(exc);

	return exc_signal;
}


std::vector<std::vector<double>> exc_vel_aleatoria(double v_min, double v_max, double S, int n, double dt) {
	double exc_max = 100;

	// Initialize velocities, amplitudes and directions
	std::vector<double> tam, v, dir;

	std::random_device rd;
	std::default_random_engine gen_uniform(rd());
	std::uniform_real_distribution<double> rand_uniform(0.0, 1.0);

	for (int i = 0; i < n; i++) {
		v.push_back(rand_uniform(gen_uniform)*(v_max - v_min) + v_min);
		tam.push_back((rand_uniform(gen_uniform) * 1.5 + 0.5)*S);
		if (rand_uniform(gen_uniform) > 0.5)
			dir.push_back(1);
		else
			dir.push_back(-1);
	}
	
	std::vector<double> t, exc;
	t.push_back(0.0);
	exc.push_back(0.0);

	for (int ct = 0; ct < int(10 / dt); ++ct) {
		exc.push_back(0);
		t.push_back(t.back() + dt);
	}
	for (int ct = 0; ct < n; ct++) {
		for (int i = 0; i < int(tam[ct] / v[ct] / dt); ++i) {
			exc.push_back(exc.back() + dir[ct] * v[ct] * dt);
			t.push_back(t.back() + dt);
			if (exc.back() > 100) {
				exc.back() = 99.9;
				dir[ct] = -1;
			}
			if (exc.back() < 0) {
				exc.back() = 0.01;
				dir[ct] = 1;
			}
		}
	}
	std::vector<std::vector<double>> exc_signal;
	exc_signal.push_back(t);
	exc_signal.push_back(exc);

	return exc_signal;
}


std::vector<std::vector<double>> exc_SP_cl_simulation(double SPm, double varSP, double Tsp, double Tsim, double dt) {

	int seed = std::chrono::system_clock::now().time_since_epoch().count() + 115499;
	std::default_random_engine gen_normal(seed);
	std::normal_distribution<double> randn(SPm, varSP);

	std::vector<double> t(int(Tsim / dt), 0.0), SP(int(Tsim / dt), 0.0);
	double SP_i;
	SP_i = randn(gen_normal);
	for (int i = 0; i < int(Tsim / dt); ++i) {
		if (std::remainder(double(i) * dt, Tsp) == 0)
			SP_i = randn(gen_normal);
		t[i] = dt * i;
		SP[i] = SP_i;
	}

	return {t, SP};
}

std::vector<double> pyLogspace(double start, double stop, int num = 50, double base = 10) {
	double realStart = pow(base, start);
	double realBase = pow(base, (stop - start) / num);

	std::vector<double> retval;
	retval.reserve(num);
	std::generate_n(std::back_inserter(retval), num, Logspace<double>(realStart, realBase));
	return retval;
}

std::vector<double> pyLinspace(double start, double stop, int num = 50) {
	std::vector<double> retval;
	double delta = (stop - start)/(double(num-1));
	retval.reserve(num);

	retval.push_back(start);
	for (int i = 0; i < num-1; ++i) {
		retval.push_back(retval.back()+delta);
	}
	return retval;
}


std::vector<double> filter1order(std::vector<double>* data, double tau, double Ts) {
	std::vector<double> filt_data;
	filt_data.reserve(data->size());
	filt_data.push_back(data->front());
	for (int i = 1; i < data->size(); ++i) {
		filt_data.push_back((Ts * data->at(i) + tau * filt_data.back()) / (Ts + tau));
	}
	return filt_data;
}


std::vector<double> simulateNoise(const std::vector<double>& data, double snr) {
	int seed = std::chrono::system_clock::now().time_since_epoch().count() + 83;
	std::default_random_engine gen_normal(seed);
	std::normal_distribution<double> randn(0.0, 1.0);

	double mean = mean_vec(data);
	double sigNoiseR = pow(10.0, snr / 10.0);
	double stanDev = mean / sigNoiseR;

	std::vector<double> retval(data.size(), 0.0);
	for (int i = 0; i < data.size(); ++i)
		retval[i] = data[i] + randn(gen_normal) * stanDev;
	return retval;
}

std::vector<double> simulateNoise(const std::vector<double>& data, double snr, double mean) {
	int seed = std::chrono::system_clock::now().time_since_epoch().count() + 83;
	std::default_random_engine gen_normal(seed);
	std::normal_distribution<double> randn(0.0, 1.0);

	double sigNoiseR = pow(10.0, snr / 10.0);
	double stanDev = mean / sigNoiseR;

	std::vector<double> retval(data.size(), 0.0);
	for (int i = 0; i < data.size(); ++i)
		retval[i] = data[i] + randn(gen_normal) * stanDev;
	return retval;
}


procDataCL preProcessCLdata(const std::vector<double>& t, const std::vector<double>& OP, 
							const std::vector<double>& P, const std::vector<double>& x, double t_exc) {
	procDataCL retval;
	
	double Ts = t[1] - t[0];
	// Find the time which the valve is set to automatic again - 1s
	int ind = int( (t_exc - 1) /  Ts );

	int sz_retval = t.size() - ind;
	retval.t.reserve(sz_retval);
	retval.OP.reserve(sz_retval);
	retval.P.reserve(sz_retval);
	retval.x.reserve(sz_retval);
	for (int i = 0; i < sz_retval; ++i) {
		retval.t.push_back(0.0);
		retval.OP.push_back(0.0);
		retval.P.push_back(0.0);
		retval.x.push_back(0.0);
	}

	int k = 0;
	for (int i = ind; i < t.size(); ++i) {
		retval.t[k] = k * Ts;
		retval.OP[k] = OP[i];
		retval.P[k] = P[i];
		retval.x[k] = x[i];
		k += 1;
	}

	int ind_stats = int((t_exc - 2) / Ts);
	std::vector<double> x_stats, P_stats;
	for (int i = ind_stats; i > ind_stats - 2.0 / Ts; --i) {
		x_stats.push_back(x[i]);
		P_stats.push_back(P[i]);
	}
	double sDevx = std_vec(x_stats);
	double x_stp = mean_vec(x_stats);
	double sDevP = std_vec(P_stats);

	std::vector<double> x_int;
	int stp_data;
	for (int i = t_exc / Ts; i > 0; --i) {
		x_int.clear();
		double end_int = (i - 0.3 / Ts < 0) ? 0 : i - 0.3 / Ts;
		for (int k = i; k > end_int; --k) {
			x_int.push_back(x[k]);
		}
		if (min_vec(abs_vec(subtract_vect_const(x_int, x_stp))) > 2 * sDevx) {
			stp_data = i;
			break;
		}
	}
	retval.x0 = x_stp;
	retval.d0 = 1;
	retval.u0 = P[stp_data];
	retval.Rv = sDevP * sDevP;
	return retval;
}
