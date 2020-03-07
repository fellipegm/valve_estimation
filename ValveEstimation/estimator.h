#pragma once


#ifndef _ESTIMATOR_H
#define _ESTIMATOR_H

#include <vector>
#include "valve_models.h"


typedef struct estimator_output
{
    std::vector<std::vector<double>> parameters;
    std::vector<double> residuals;
    int nFuncEvals{ 0 };
} estimator_output;


class Estimator {
private:
    // data to match with the model
    std::vector<double> des_data;

    std::vector<double> ub;
    std::vector<double> lb;

	int DE_pop;

    bool pen_search_space = true;

	// Optimization algorithm minimum decrease after 20% of iterations
    double mindecrease { pow(1e-7,2) };
 
    double residual_calc(ValveModel* model);
    estimator_output initial_map();
    std::vector<std::vector<double>> find_combinations();
	estimator_output dif_evolution(std::vector<std::vector<double>>);
	estimator_output simplex(std::vector<double>);

public:
	ValveModel valve;

	Estimator();
    void set_des_data(std::vector<double>);
	void set_pen_search_space(bool pen) { pen_search_space = pen; };
    void set_lb_ub(std::vector<double> lb_input, std::vector<double> ub_input) { lb = lb_input; ub = ub_input; };
    void calc_lbub(double S0, double std_k);

    estimator_output run_estimator();
};


#endif //_ESTIMATOR_H