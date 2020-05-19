

#include <iostream>
#include <string>
#include <math.h>
#include <chrono>

#include "valve_models.h"
#include "estimator.h"
#include "utils.h"
#include <algorithm>
#include <random>
#include <omp.h>

Estimator::Estimator() {
	valve = ValveModel();
}

void Estimator::calc_lbub(double S0, double std_k, bool ident_k_finit) {
	double Fc_Fs_u2 = (valve.get_pmax() * valve.get_Sa() - valve.get_pmin() * valve.get_Sa()) * S0 / 100 / 2;
	double Fc_Fs_2pc = (valve.get_pmax() * valve.get_Sa() - valve.get_pmin() * valve.get_Sa()) / 100;
	double Fc_ub, Fc_lb, Fs_ub, Fs_lb, Fv_lb, Fv_ub, S_lb, S_ub, J_lb, J_ub, D_lb, D_ub, vs_lb, vs_ub;

	double k_lb, k_ub, finit_lb, finit_ub;
	double threshold_error_k = 0.1;
	if (std::abs((valve.get_k() - 2 * std_k) / valve.get_k() - 1) > threshold_error_k) {
		double error_k = std::max(0.1, std::abs((valve.get_k() - 2 * std_k) / valve.get_k() - 1));
		k_lb = valve.get_k() - valve.get_k() * error_k;
		k_ub = valve.get_k() + valve.get_k() * error_k;
		finit_lb = valve.get_Finit() - valve.get_Finit() * error_k;
		finit_ub = valve.get_Finit() + valve.get_Finit() * error_k;
	}
	else {
		k_lb = valve.get_k() - valve.get_k() * threshold_error_k;
		k_ub = valve.get_k() + valve.get_k() * threshold_error_k;
		finit_lb = valve.get_Finit() - valve.get_Finit() * threshold_error_k;
		finit_ub = valve.get_Finit() + valve.get_Finit() * threshold_error_k;
	}
	if (valve.get_Finit() < 0) {
		double aux_invert = finit_lb;
		finit_lb = finit_ub;
		finit_ub = aux_invert;
	}

	if (S0 > 2) {
		// Bounds for every first principles model
		Fc_ub = Fc_Fs_u2 * 1.3;
		Fc_lb = Fc_Fs_u2 * 0.49;
		Fs_ub = Fc_Fs_u2 * 1.3;
		Fs_lb = Fc_Fs_u2 * 0.7;
		Fv_lb = Fc_ub / 100 / ((valve.get_xmax() - 0) * 2 * M_PI / valve.get_tauip());

		// Kano bounds
		S_ub = S0 * 1.3;
		S_lb = S0 * 0.7;
		J_ub = S0 * 0.4;
	}
	else {
		Fc_ub = Fc_Fs_2pc * 1.3;
		Fc_lb = 0;
		Fs_ub = Fc_Fs_2pc * 1.3;
		Fs_lb = 0;
		Fv_lb = 0;

		// Kano bounds
		S_ub = 2;
		S_lb = 0;
		J_ub = 2 * 0.4;
	}

	// Kano bounds
	J_lb = 0;
	D_ub = 30;
	D_lb = 0;

	// Upper bound for viscous friction coefficient
	Fv_ub = Fc_ub * 100 / ((valve.get_xmax() - 0) * 2 * M_PI / valve.get_tauip());
	// Bounds for Stribeck velocity
	vs_ub = 0.1;
	vs_lb = 1e-5;

	//Bounds for LuGre
	double sigma0_ub = Fc_Fs_u2 * 1.3 * 1e6;
	double sigma0_lb = Fc_Fs_u2 * 0.7 * 1e3;
	double sigma1_ub = (sigma0_ub + sigma0_lb)/2;
	double sigma1_lb = 500;

	// Bounds for GMS
	double kappa_ub = Fc_Fs_u2 * 1.3 * 1e6;
	double kappa_lb = Fc_Fs_u2 * 0.7 * 1e3;
	double nu_ub = (kappa_ub + kappa_lb) / 2;
	double nu_lb = 500;
	double alpha_ub = 1;
	double alpha_lb = 0;
	double C_ub = 1e3;
	double C_lb = 1;

	switch (valve.get_model()) {
	case friction_model::kano:
		if (ident_k_finit) {
			double range = (valve.get_xmax() - valve.get_xmin()) / 2;
			set_lb_ub({ S_lb, J_lb, D_lb, std::max(0.0, valve.get_xmin() - range), std::max(0.0, valve.get_xmax() - range) },
				{ S_ub, J_ub, D_ub, valve.get_xmin() + range, valve.get_xmax() + range });
		}
		else
			set_lb_ub({ S_lb, J_lb, D_lb }, { S_ub, J_ub, D_ub });
		break;
	case friction_model::choudhury:
		if (ident_k_finit) {
			double range = (valve.get_xmax() - valve.get_xmin()) / 2;
			set_lb_ub({ S_lb, J_lb, D_lb, std::max(0.0, valve.get_xmin() - range), std::max(0.0, valve.get_xmax() - range) },
				{ S_ub, J_ub, D_ub, valve.get_xmin() + range, valve.get_xmax() + range });
		}
		else
			set_lb_ub({ S_lb, J_lb, D_lb }, { S_ub, J_ub, D_ub });
		break;
	case friction_model::he:
		if (ident_k_finit) {
			double range = (valve.get_xmax() - valve.get_xmin()) / 2;
			set_lb_ub({ 0.35*S0, 0.0, D_lb, std::max(0.0, valve.get_xmin() - range), std::max(0.0, valve.get_xmax() - range) },
				{ 0.845*S0, 0.845*S0, D_ub, valve.get_xmin() + range, valve.get_xmax() + range });
		}
		else
			set_lb_ub({ 0.35 * S0, 0.0, D_lb }, { 0.845 * S0, 0.845 * S0, D_ub });
		break;
	case friction_model::karnopp:
		if (ident_k_finit) {
			set_lb_ub({ k_lb, finit_lb, Fc_lb, Fs_lb, Fv_lb, vs_lb },
				{ k_ub, finit_ub, Fc_ub, Fs_ub, Fv_ub, vs_ub });
		}
		else
			set_lb_ub({ Fc_lb, Fs_lb, Fv_lb, vs_lb }, { Fc_ub, Fs_ub, Fv_ub, vs_ub });
		break;
	case friction_model::lugre:
		if (ident_k_finit) {
			set_lb_ub({ k_lb, finit_lb, Fc_lb, Fs_lb, Fv_lb, vs_lb, sigma0_lb, sigma1_lb },
				{ k_ub, finit_ub, Fc_ub, Fs_ub, Fv_ub, vs_ub, sigma0_ub, sigma1_ub });
		}
		else
			set_lb_ub({ Fc_lb, Fs_lb, Fv_lb, vs_lb, sigma0_lb, sigma1_lb }, { Fc_ub, Fs_ub, Fv_ub, vs_ub, sigma0_ub, sigma1_ub });
		break;
	case friction_model::gms:
		if (ident_k_finit) {
			set_lb_ub({ k_lb, finit_lb, Fc_lb, Fs_lb, Fv_lb, vs_lb, kappa_lb, kappa_lb, kappa_lb, nu_lb, nu_lb, nu_lb, alpha_lb, alpha_lb, C_lb },
				{ k_ub, finit_ub, Fc_ub, Fs_ub, Fv_ub, vs_ub, kappa_ub, kappa_ub, kappa_ub, nu_ub, nu_ub, nu_ub, alpha_ub, alpha_ub, C_ub });
		}
		else
			set_lb_ub({ Fc_lb, Fs_lb, Fv_lb, vs_lb, kappa_lb, kappa_lb, kappa_lb, nu_lb, nu_lb, nu_lb, alpha_lb, alpha_lb, C_lb },
				{ Fc_ub, Fs_ub, Fv_ub, vs_ub, kappa_ub, kappa_ub, kappa_ub, nu_ub, nu_ub, nu_ub, alpha_ub, alpha_ub, C_ub });
		break;
	case friction_model::sgms:
		if (ident_k_finit) {
			set_lb_ub({ k_lb, finit_lb, Fc_lb, Fs_lb, Fv_lb, vs_lb, kappa_lb, kappa_lb, alpha_lb, C_lb },
				{ k_ub, finit_ub, Fc_ub, Fs_ub, Fv_ub, vs_ub, kappa_ub, kappa_ub, alpha_ub, C_ub });
		}
		else
			set_lb_ub({ Fc_lb, Fs_lb, Fv_lb, vs_lb, kappa_lb, kappa_lb, alpha_lb, C_lb },
				{ Fc_ub, Fs_ub, Fv_ub, vs_ub, kappa_ub, kappa_ub, alpha_ub, C_ub });
		break;
	case friction_model::gms1:
		if (ident_k_finit) {
			set_lb_ub({ k_lb, finit_lb, Fc_lb, Fs_lb, Fv_lb, vs_lb, kappa_lb, nu_lb, C_lb },
				{ k_ub, finit_ub, Fc_ub, Fs_ub, Fv_ub, vs_ub, kappa_ub, nu_ub, C_ub });
		}
		else
			set_lb_ub({ Fc_lb, Fs_lb, Fv_lb, vs_lb, kappa_lb, nu_lb, C_lb },
				{ Fc_ub, Fs_ub, Fv_ub, vs_ub, kappa_ub, nu_ub, C_ub });
		break;
	default:
		break;
	}
}

void Estimator::set_des_data(std::vector<double> input_data) {
	des_data = input_data;

	mindecrease = pow(1e-7, 2) / 50;
}

double Estimator::residual_calc(ValveModel* model) {

	if (model->simulation_results.x.size() != des_data.size()) {
		std::cout << "Error: Simulated data and desired data does not have the same length" << std::endl;
		std::exit;
	}

	double residual{ 0 };

	for (int i = 0; i < des_data.size(); i++) {
		residual += std::pow(model->simulation_results.x[i] - des_data[i], 2);
	}
	residual = residual / ((double)des_data.size());

	double penalization{ 0.0 };
	if (pen_search_space) {
		for (size_t i = 0; i < lb.size(); i++)
			penalization += std::max(0.0, lb[i] - model->get_param_friction(i)) * 1e5 + std::max(0.0, model->get_param_friction(i) - ub[i]) * 1e5;
		if (model->get_model() != friction_model::kano && model->get_model() != friction_model::choudhury)
			penalization += std::max(0.0, model->get_Fc() - model->get_Fs()) * 1e5;
		if (model->get_model() == friction_model::gms) {
			double alpha_error = std::abs(model->get_alpha1() + model->get_alpha2() + model->get_alpha3() - 1);
			if (std::abs(alpha_error) > 1e-3)
				penalization += std::abs(alpha_error) * 1e5;
			if (model->get_alpha1() < 0.0 || model->get_alpha1() > 1.0)
				penalization += std::abs(model->get_alpha1()) * 1e5;
			if (model->get_alpha2() < 0.0 || model->get_alpha2() > 1.0)
				penalization += std::abs(model->get_alpha2()) * 1e5;
			if (model->get_alpha3() < 0.0 || model->get_alpha3() > 1.0)
				penalization += std::abs(model->get_alpha3()) * 1e5;
		}
		if (model->get_model() == friction_model::sgms) {
			if (model->get_alpha1() < 0 || model->get_alpha1() > 1)
				penalization += std::abs(model->get_alpha1()) * 1e5;
		}
	}
	return (isnan(residual + penalization) ? 1e20 : residual + penalization);
}

estimator_output Estimator::run_estimator() {
	// Run a grid search between the parameters bounds
	estimator_output grid_search = initial_map();

	// Initialize data for DE optimization
	DE_pop = lb.size() * 10;

	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine gen_uniform(seed);
	std::uniform_real_distribution<double> rand_uniform(0.0, 1.0);

	seed = std::chrono::system_clock::now().time_since_epoch().count() + 2458;
	std::default_random_engine gen_normal(seed);
	std::normal_distribution<double> randn(0.0, 1.0);

	std::vector<double> aux_pop_vector;
	std::vector<std::vector<double>> init_pop;
	for (int i = 0; i < DE_pop; i++) {
		if (i < 10 && i < grid_search.parameters.size())
			init_pop.push_back(grid_search.parameters[i]);
		else if (i < 20 && i < grid_search.parameters.size()) {
			aux_pop_vector.clear();
			for (int j = 0; j < lb.size(); j++) {
				aux_pop_vector.push_back(randn(gen_normal) * (ub[j] - lb[j]) / 100 + grid_search.parameters[i][j]);
				if (aux_pop_vector.back() < lb[j])
					aux_pop_vector.back() = lb[j];
				else if (aux_pop_vector.back() > ub[j])
					aux_pop_vector.back() = ub[j];
			}
			init_pop.push_back(aux_pop_vector);
		}
		else {
			aux_pop_vector.clear();
			for (int j = 0; j < lb.size(); j++)
				aux_pop_vector.push_back(rand_uniform(gen_uniform) * (ub[j] - lb[j]) + lb[j]);
			init_pop.push_back(aux_pop_vector);
		}
		if (valve.get_model() == friction_model::gms) {
			int offset = 0;
			if (init_pop[i].size() >= 15)
				offset = 2;
			while (init_pop[i][10 + offset] + init_pop[i][11 + offset] > 1) {
				init_pop[i][10 + offset] = rand_uniform(gen_uniform);
				init_pop[i][11 + offset] = rand_uniform(gen_uniform);
			}
		}
	}
	// Run a DE optimization
	estimator_output de_results = dif_evolution(init_pop);

	// Run a Simplex optimization
	estimator_output simplex_results = simplex(de_results.parameters[0]);

	simplex_results.nFuncEvals = grid_search.nFuncEvals + de_results.nFuncEvals + simplex_results.nFuncEvals;
	return simplex_results;

}


estimator_output Estimator::simplex(std::vector<double> param0) {
	std::cout << "Simplex Optimization \n\n";
	// Parameters dimension
	int dim = param0.size();

	// Initialize optimization parameters
	double rho{ 1 }, chi{ 2 }, psi{ 0.5 }, sigma{ 0.5 };

	// Maximum number of iterations
	int itermax = 1500;

	// Number of function evals
	int nfeval = 0;

	// Initialize vertices and residues
	std::vector<std::vector<double>> v;
	v.push_back(param0);
	std::vector<double> residual_v;

	// Generate vertices
	std::vector<double> aux_vert_init;
	for (int i = 0; i < dim; i++) {
		if (std::abs(param0[i]) > 1e-10) {
			aux_vert_init = param0;
			aux_vert_init[i] = param0[i] + (ub[i] - lb[i]) * 0.05; // five percent vertice creation
			v.push_back(aux_vert_init);
			aux_vert_init.clear();
		}
		else {
			aux_vert_init = param0;
			aux_vert_init[i] = (ub[i] - lb[i]) * 0.00025 + 0.00025; // zero term delta
			v.push_back(aux_vert_init);
			aux_vert_init.clear();
		}
	}

	// Evaluate the vertices
	for (int i = 0; i < v.size(); ++i)
		residual_v.push_back(0.0);

	ValveModel* models = new ValveModel[v.size()];
#pragma omp parallel for
	for (int i = 0; i < v.size(); i++) {
		models[i] = ValveModel();
		models[i] = valve;
		models[i].set_friction_param_value(v[i]);
		models[i].valve_simulation();
		residual_v[i] = residual_calc(&models[i]);
	}
	nfeval += dim + 1;

	std::vector<size_t> idx = sort_indexes(residual_v);
	std::vector<std::vector<double>> sorted_v;
	sort(residual_v.begin(), residual_v.end());
	for (int i = 0; i < idx.size(); i++) {
		sorted_v.push_back(v[idx[i]]);
	}
	v = sorted_v;

	delete[] models;

	// Simplex optimization algorithm
	int iteration = 0;
	std::cout << "Iteration: " << iteration << std::endl;
	std::cout << "Residual: " << residual_v[0] << std::endl;
	std::cout << "Procedure: " << "Initialization" << std::endl;
	for (int i = 0; i < v[0].size(); ++i)
		std::cout << "Best_par[" << i << "]: " << v[0][i] << std::endl;
	std::cout << std::endl;

	std::vector<double> hist_residual;
	hist_residual.push_back(residual_v[0]);
	std::vector<std::vector<double>> aux_v_contr, xcont, x;
	std::vector<double> xbar, xr, xe, xc, xcc, residual_x;
	for (int i = 0; i < dim + 4; ++i)
		residual_x.push_back(0.0);

	// Array of valve models to use parallel threads
	ValveModel* models_x = new ValveModel[dim + 4];
	for (int i = 0; i < dim + 4; i++) {
		x.push_back(v[0]);
		models_x[i] = ValveModel();
		models_x[i] = valve;
	}

	std::vector<double> fxcont;
	std::string procedure;
	while (iteration < itermax) {
		if (iteration > int(double(itermax) * 0.2) + 1) {
			if (hist_residual[iteration - int(double(itermax) * 0.1)] - hist_residual[iteration - 1] < mindecrease)
				break;
		}

		// Average of the best n points
		xbar.clear();
		for (int i = 0; i < dim; ++i) {
			double mean_aux = 0;
			for (int j = 0; j < dim; ++j)
				mean_aux += v[j][i] / double(dim);
			xbar.push_back(mean_aux);
		}
		xr = subtract_vectors(vector_const_multiply(xbar, 1 + rho), vector_const_multiply(v.back(), rho));
		xe = subtract_vectors(vector_const_multiply(xbar, 1 + rho * chi), vector_const_multiply(v.back(), rho * chi));
		xc = subtract_vectors(vector_const_multiply(xbar, 1 + psi * rho), vector_const_multiply(v.back(), psi * rho));
		xcc = subtract_vectors(vector_const_multiply(xbar, 1 - psi), vector_const_multiply(v.back(), psi));
		xcont.clear();
		for (int i = 0; i < dim; ++i)
			xcont.push_back(sum_vectors(v[0], vector_const_multiply(subtract_vectors(v[i + 1], v[0]), sigma)));
		x.clear();
		x.push_back(xr);
		x.push_back(xe);
		x.push_back(xc);
		x.push_back(xcc);
		for (int i = 0; i < dim; i++)
			x.push_back(xcont[i]);

#pragma omp parallel for
		for (int i = 0; i < x.size(); i++) {
			models_x[i].set_friction_param_value(x[i]);
			models_x[i].valve_simulation();
			residual_x[i] = residual_calc(&models_x[i]);
		}
		nfeval += x.size();

		double fxr = residual_x[0];
		double fxe = residual_x[1];
		double fxc = residual_x[2];
		double fxcc = residual_x[3];
		fxcont.clear();
		for (int i = 4; i < residual_x.size(); i++)
			fxcont.push_back(residual_x[i]);

		if (fxr < residual_v[0]) {
			if (fxe < fxr) {
				v[dim] = xe;
				residual_v[dim] = fxe;
				procedure = "expand";
			}
			else {
				v[dim] = xr;
				residual_v[dim] = fxr;
				procedure = "reflect";
			}
		}
		else {
			if (fxr < residual_v[dim - 1]) {
				v[dim] = xr;
				residual_v[dim] = fxr;
				procedure = "reflect";
			}
			else {
				if (fxr < residual_v[dim]) {
					if (fxc <= fxr) {
						v[dim] = xc;
						residual_v[dim] = fxc;
						procedure = "contract outside";
					}
					else
						procedure = "shrink";
				}
				else {
					if (fxcc < residual_v[dim]) {
						v[dim] = xcc;
						residual_v[dim] = fxcc;
						procedure = "contract inside";
					}
					else
						procedure = "shrink";
				}
				if (procedure.compare("shrink") == 0) {
					for (int i = 1; i < dim + 1; ++i) {
						v[i] = xcont[i - 1];
						residual_v[i] = fxcont[i - 1];
					}
				}
			}
		}
		idx.clear();
		sorted_v.clear();
		idx = sort_indexes(residual_v);
		sort(residual_v.begin(), residual_v.end());
		for (int i = 0; i < idx.size(); i++) {
			sorted_v.push_back(v[idx[i]]);
		}
		v = sorted_v;
		hist_residual.push_back(residual_v[0]);

		iteration += 1;
		std::cout << "Iteration: " << iteration << std::endl;
		std::cout << "Residual: " << residual_v[0] << std::endl;
		std::cout << "Procedure: " << procedure << std::endl;
		for (int i = 0; i < v[0].size(); ++i)
			std::cout << "Best_par[" << i << "]: " << v[0][i] << std::endl;
		std::cout << std::endl;
	}

	delete[] models_x;

	estimator_output results;
	results.parameters.push_back(v[0]);
	results.nFuncEvals = nfeval;
	results.residuals.push_back(residual_v[0]);
	return results;
}

estimator_output Estimator::dif_evolution(std::vector<std::vector<double>> initial_pop) {
	std::cout << "Differential Evolution optimization \n\n";
	// Crossover probability
	double CR = 0.7;
	// Diferential evolution weight
	double DE_weight = 0.8;
	// Maximum number of iterations
	int itermax = 1500;

	// Parameters dimension
	int dim = initial_pop[0].size();

	// Initialize previous population
	std::vector<std::vector<double>> pop_old = initial_pop;

	// Evaluate the objective function for each population member
	std::vector<double> residual_values_old;
	residual_values_old.reserve(pop_old.size());
	for (int i = 0; i < pop_old.size(); ++i) residual_values_old.push_back(0.0);

	ValveModel* models = new ValveModel[pop_old.size()];
	#pragma omp parallel for
	for (int i = 0; i < pop_old.size(); i++) {
		models[i] = ValveModel();
		models[i] = valve;
		models[i].set_friction_param_value(pop_old[i]);
		models[i].valve_simulation();
		residual_values_old[i] = residual_calc(&models[i]);
	}

	// Initialize uniform random number generator
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine gen_uniform(seed);
	std::uniform_real_distribution<double> rand_uniform(0.0, 1.0);

	unsigned seed_shuf = std::chrono::system_clock::now().time_since_epoch().count() + 554887;
	std::default_random_engine e_shuf(seed);


	// Actual population and masks matrices initialization
	std::vector<std::vector<double>> mask, inv_mask, pop, pop_shuf1, pop_shuf2, pop_shuf3, pop_dif;
	std::vector<double> residual_values;
	for (int i = 0; i < pop_old.size(); ++i)
		residual_values.push_back(0.0);
	std::vector<int> idx_pop, rot_nbr;
	idx_pop.reserve(DE_pop);
	for (int i = 0; i < DE_pop; i++) { idx_pop.push_back(i); }
	for (int i = 1; i < 5; i++) rot_nbr.push_back(i); // initialize vector with 1 to 4

	// Initialize best residual and parameters
	std::vector<double> best_ind;
	double best_residual;
	auto idx_best_residual = std::min_element(residual_values_old.begin(), residual_values_old.end()) - residual_values_old.begin();
	best_residual = residual_values_old[idx_best_residual];
	best_ind = pop_old[idx_best_residual];

	std::vector<double> hist_residual;
	hist_residual.push_back(best_residual);

	std::vector<double> member_mask_temp, member_inv_mask_temp;

	// Differential evolution
	int iteration = 1;
	while (iteration < itermax) {
		// Check end criteria - the residual is not decreasing
		if (iteration > int(double(itermax) * 0.2) + 1) {
			if (hist_residual[iteration - int(double(itermax) * 0.2)] - hist_residual[iteration - 1] < mindecrease)
				break;
		}

		// Create vectors of crossover probability
		mask.clear();
		inv_mask.clear();

		for (int i = 0; i < initial_pop.size(); i++) {
			member_mask_temp.clear();
			member_inv_mask_temp.clear();
			for (int j = 0; j < initial_pop[0].size(); j++) {
				if (rand_uniform(gen_uniform) < CR) {
					member_mask_temp.push_back(1);
					member_inv_mask_temp.push_back(0);
				}
				else {
					member_mask_temp.push_back(0);
					member_inv_mask_temp.push_back(1);
				}
			}
			mask.push_back(member_mask_temp);
			inv_mask.push_back(member_inv_mask_temp);
		}

		// Create shuffled populations
		std::shuffle(std::begin(idx_pop), std::end(idx_pop), e_shuf);
		std::shuffle(std::begin(rot_nbr), std::end(rot_nbr), e_shuf);

		pop_shuf1.clear();
		pop_shuf2.clear();
		pop_shuf3.clear();
		for (int i = 0; i < DE_pop; i++) {
			pop_shuf1.push_back(pop_old[idx_pop[(i + rot_nbr[0]) % DE_pop]]);
			pop_shuf2.push_back(pop_old[idx_pop[(i + rot_nbr[0] + rot_nbr[1]) % DE_pop]]);
			pop_shuf3.push_back(pop_old[idx_pop[(i + rot_nbr[0] + rot_nbr[1] + rot_nbr[2]) % DE_pop]]);
		}

		// Differential variation
		pop_dif = sum_matrices(pop_shuf3, matrix_const_multiply(subtract_matrices(pop_shuf1, pop_shuf2), (1 - DE_weight) * rand_uniform(gen_uniform) + DE_weight));
		// Crossover
		pop = sum_matrices(matrix_element_multiply(pop_dif, mask), matrix_element_multiply(pop_old, inv_mask));

		// Ensure the parameters constraints
		for (int i = 0; i < pop.size(); i++) {
			for (int j = 0; j < pop[0].size(); j++) {
				if (pop[i][j] > ub[j])
					pop[i][j] = ub[j] + rand_uniform(gen_uniform) * (pop_shuf3[i][j] - ub[j]);
				if (pop[i][j] < lb[j])
					pop[i][j] = lb[j] + rand_uniform(gen_uniform) * (pop_shuf3[i][j] - lb[j]);
			}
			if (valve.get_model() == friction_model::gms) {
				int offset = 0;
				if (lb.size() >= 15)
					offset = 2;
				while (pop[i][10 + offset] + pop[i][11 + offset] > 1) {
					pop[i][10 + offset] = rand_uniform(gen_uniform);
					pop[i][11 + offset] = rand_uniform(gen_uniform);
				}
			}
		}

		// Evaluate the MSE for every population member
#pragma omp parallel for
		for (int i = 0; i < pop.size(); i++) {
			models[i].set_friction_param_value(pop[i]);
			models[i].valve_simulation();
			residual_values[i] = residual_calc(&models[i]);
		}

		// Create new population
		for (int i = 0; i < DE_pop; i++) {
			if (residual_values[i] < residual_values_old[i]) {
				residual_values_old[i] = residual_values[i];
				pop_old[i] = pop[i];
				if (residual_values[i] < best_residual) {
					best_residual = residual_values[i];
					best_ind = pop[i];
				}
			}
		}
		hist_residual.push_back(best_residual);

		// Print data
		std::cout << "Iteration: " << iteration << std::endl;
		std::cout << "Best residual: " << best_residual << std::endl;
		for (int i = 0; i < best_ind.size(); i++)
			std::cout << "Best_par[" << i << "]: " << best_ind[i] << std::endl;
		std::cout << std::endl;
		iteration += 1;
	}

	delete[] models;

	estimator_output results;
	results.parameters.push_back(best_ind);
	results.nFuncEvals = iteration * DE_pop;
	results.residuals.push_back(best_residual);
	return results;
}


estimator_output Estimator::initial_map() {
	// FUNCTION THAT RETURNS THE MAPPING OF THE GRID SEARCH
	std::cout << "Grid Search \n\n";
	std::vector<std::vector<double>> combinations = find_combinations();
	std::vector<double> mapping_results;

	mapping_results.reserve(combinations.size());
	for (int i = 0; i < combinations.size(); ++i)
		mapping_results.push_back(0.0);

	int size = 16 * 20;
	ValveModel* models = new ValveModel[size];

	int ct{ 0 };
	while (ct < combinations.size()) {
		if (ct + size > combinations.size())
			size = combinations.size() - ct;
		#pragma omp parallel for
		for (int i = 0; i < size; i++) {
			models[i] = ValveModel();
			models[i] = valve;
			models[i].set_friction_param_value(combinations[i + ct]);
			models[i].valve_simulation();
			mapping_results[i + ct] = residual_calc(&models[i]);
			models[i].clear_sim_data();
			std::cout << "Progress: " << double(i + ct) / double(combinations.size()) * 100 << "\t" << "Residuals: " << mapping_results[i + ct] << std::endl;
		}
		ct = ct + size;
	}


	delete[] models;

	std::vector<size_t> idx = sort_indexes(mapping_results);
	std::vector<std::vector<double>> sorted_combinations;
	sort(mapping_results.begin(), mapping_results.end());
	for (int i = 0; i < idx.size(); i++) {
		sorted_combinations.push_back(combinations[idx[i]]);
	}
	estimator_output results;
	results.parameters = sorted_combinations;
	results.residuals = mapping_results;
	results.nFuncEvals = mapping_results.size();

	std::cout << std::endl << std::endl;
	return results;
}


std::vector<std::vector<double>> Estimator::find_combinations() {
	// THIS FUNCTION CREATES THE GRID SEARCH
	std::vector<std::vector<double>> mat_ranges;
	std::vector<double> param_vect;
	std::vector<size_t> size_data;
	size_t matrix_rows{ 1 };

	std::vector<double> lb_aux, ub_aux;
	lb_aux = lb;
	ub_aux = ub;

	// The case where the k and f_init are estimated toghether
	bool estimate_k_finit = false;
	if ((valve.get_model() == friction_model::karnopp && ub.size() == 6) || (valve.get_model() == friction_model::lugre && ub.size() == 8) ||
		(valve.get_model() == friction_model::gms && ub.size() == 15) || (valve.get_model() == friction_model::sgms && ub.size() == 10) ||
		(valve.get_model() == friction_model::gms1 && ub.size() == 9)) {
		estimate_k_finit = true;
		lb_aux.erase(lb_aux.begin(), lb_aux.begin() + 2);
		ub_aux.erase(ub_aux.begin(), ub_aux.begin() + 2);
	}

	if (valve.get_model() == friction_model::gms) {
		lb_aux.erase(lb_aux.begin() + 10, lb_aux.begin() + 12);
		lb_aux.erase(lb_aux.begin() + 7, lb_aux.begin() + 9);
		lb_aux.erase(lb_aux.begin() + 5, lb_aux.begin() + 7);
		ub_aux.erase(ub_aux.begin() + 10, ub_aux.begin() + 12);
		ub_aux.erase(ub_aux.begin() + 7, ub_aux.begin() + 9);
		ub_aux.erase(ub_aux.begin() + 5, ub_aux.begin() + 7);
	}

	if (valve.get_model() == friction_model::sgms) {
		lb_aux.erase(lb_aux.begin() + 5, lb_aux.begin() + 7);
		ub_aux.erase(ub_aux.begin() + 5, ub_aux.begin() + 7);
	}

	if (valve.get_model() == friction_model::gms1) {
		lb_aux.erase(lb_aux.begin() + 6, lb_aux.begin() + 7);
		ub_aux.erase(ub_aux.begin() + 6, ub_aux.begin() + 7);
	}

	double paramAux;
	// create the points of each parameter
	for (size_t i = 0; i < lb_aux.size(); i++) {
		param_vect.clear();
		if (log10(ub_aux[i] / lb_aux[i]) < 2)
			param_vect = pyLinspace(lb_aux[i], ub_aux[i], 3);
		else {
			param_vect.push_back(lb_aux[i]);
			int ct = 1;
			while (param_vect.back() < ub_aux[i]) {
				(param_vect.back() == 0) ? paramAux = 1 : paramAux = 0;
				if (param_vect.back() < 0 && param_vect.back() < -1)
					param_vect.push_back((param_vect.back() + paramAux) / 10);
				else
					param_vect.push_back(std::abs(param_vect.back() + paramAux) * 10);
				ct++;
			}
			param_vect.back() = ub_aux[i];
		}
		size_data.push_back(param_vect.size());
		mat_ranges.push_back(param_vect);
		matrix_rows *= size_data[i];
	}

	std::vector<double> v(lb_aux.size(), 0);
	std::vector<int> v_int(lb_aux.size(), 0);
	std::vector<std::vector<double>> combinations(matrix_rows, v);
	std::vector<std::vector<int>> i_combs(matrix_rows, v_int);

	// find all possible combinations by counting
	int size_div = 1;
	int ct;
	for (size_t j = 0; j < lb_aux.size(); j++) {
		size_div *= mat_ranges[j].size();
		ct = 0;
		for (size_t i = 0; i < matrix_rows; i++) {
			i_combs[i][j] = ct;
			combinations[i][j] = mat_ranges[j][ct];
			if ((i + 1) % (matrix_rows / size_div) == 0)
				ct++;
			if (ct > mat_ranges[j].size() - 1)
				ct = 0;
		}
	}

	if (valve.get_model() != friction_model::kano && valve.get_model() != friction_model::he && valve.get_model() != friction_model::choudhury) {
		for (int i = 0; i < combinations.size(); i++) {
			if (combinations[i][1] < combinations[i][0]) {
				combinations.erase(combinations.begin() + i);
				i--;
			}
		}
	}

	if (valve.get_model() == friction_model::gms) {
		unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
		std::default_random_engine gen_uniform(seed);
		std::uniform_real_distribution<double> rand_uniform(0.0, 1.0);

		seed = std::chrono::system_clock::now().time_since_epoch().count() + 112778;
		std::default_random_engine gen_normal(seed);
		std::normal_distribution<double> randn(0.0, 1.0);
		int comp_idx;
		(estimate_k_finit) ? comp_idx = 2 : comp_idx = 0;
		double kappa2_aux{ ub[4 + comp_idx] * 2 }, kappa3_aux{ ub[4 + comp_idx] * 2 }, nu2_aux{ ub[7 + comp_idx] * 2 }, nu3_aux{ ub[7 + comp_idx] * 2 }, alpha1_aux{ 1.0 }, alpha2_aux{ 1.0 };
		for (int i = 0; i < combinations.size(); i++) {
			while (kappa2_aux > ub[4 + comp_idx] || kappa2_aux < lb[4 + comp_idx])
				kappa2_aux = combinations[i][4] * randn(gen_normal) + combinations[i][4];
			while (kappa3_aux > ub[4 + comp_idx] || kappa3_aux < lb[4 + comp_idx])
				kappa3_aux = combinations[i][4] * randn(gen_normal) + combinations[i][4];
			while (nu2_aux > ub[7 + comp_idx] || nu2_aux < lb[7 + comp_idx])
				nu2_aux = combinations[i][5] * randn(gen_normal) + combinations[i][5];
			while (nu3_aux > ub[7 + comp_idx] || nu3_aux < lb[7 + comp_idx])
				nu3_aux = combinations[i][5] * randn(gen_normal) + combinations[i][5];
			while (alpha1_aux + alpha2_aux > 1.0) {
				alpha1_aux = rand_uniform(gen_uniform);
				alpha2_aux = rand_uniform(gen_uniform);
			}
			combinations[i].insert(combinations[i].begin() + 5, kappa2_aux);
			combinations[i].insert(combinations[i].begin() + 5, kappa3_aux);
			combinations[i].insert(combinations[i].begin() + 8, nu2_aux);
			combinations[i].insert(combinations[i].begin() + 8, nu3_aux);
			combinations[i].insert(combinations[i].begin() + 10, alpha1_aux);
			combinations[i].insert(combinations[i].begin() + 10, alpha2_aux);
			kappa2_aux = ub[4 + comp_idx] * 2;
			kappa3_aux = ub[4 + comp_idx] * 2;
			nu2_aux = ub[7 + comp_idx] * 2;
			nu3_aux = ub[7 + comp_idx] * 2;
			alpha1_aux = 1.0;
			alpha2_aux = 1.0;
		}
	}

	if (valve.get_model() == friction_model::sgms) {
		unsigned seed = std::chrono::system_clock::now().time_since_epoch().count() + 142331;
		std::default_random_engine gen_uniform(seed);
		std::uniform_real_distribution<double> rand_uniform(0.0, 1.0);

		seed = std::chrono::system_clock::now().time_since_epoch().count() + 10154;
		std::default_random_engine gen_normal(seed);
		std::normal_distribution<double> randn(0.0, 1.0);

		int comp_idx;
		(estimate_k_finit) ? comp_idx = 2 : comp_idx = 0;
		double kappa2_aux{ ub[4 + comp_idx] * 2 }, alpha1_aux{ 1.0 };
		for (int i = 0; i < combinations.size(); i++) {
			while (kappa2_aux > ub[4 + comp_idx] || kappa2_aux < lb[4 + comp_idx])
				kappa2_aux = combinations[i][4] * randn(gen_normal) + combinations[i][4];
			alpha1_aux = rand_uniform(gen_uniform);

			combinations[i].insert(combinations[i].begin() + 5, kappa2_aux);
			combinations[i].insert(combinations[i].begin() + 6, alpha1_aux);
			kappa2_aux = ub[4 + comp_idx] * 2;
		}
	}

	if (valve.get_model() == friction_model::gms1) {
		unsigned seed = std::chrono::system_clock::now().time_since_epoch().count() + 142331;
		std::default_random_engine gen_uniform(seed);
		std::uniform_real_distribution<double> rand_uniform(0.0, 1.0);

		int comp_idx;
		double C_aux;
		(estimate_k_finit) ? comp_idx = 2 : comp_idx = 0;
		for (int i = 0; i < combinations.size(); i++) {
			C_aux = rand_uniform(gen_uniform) * (ub[6 + comp_idx] - lb[6 + comp_idx]) + lb[6 + comp_idx];

			combinations[i].insert(combinations[i].begin() + 6, C_aux);
		}
	}

	// The case where the k and f_init are estimated toghether
	if (estimate_k_finit) {
		for (int i = 0; i < combinations.size(); i++) {
			combinations[i].insert(combinations[i].begin(), valve.get_Finit());
			combinations[i].insert(combinations[i].begin(), valve.get_k());
		}
	}

	return combinations;
}