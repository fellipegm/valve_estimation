

#include "controller.h"
#include "utils.h"
#include <vector>


double Controller::pid(double SP, double PV, int ct){

	if (ct < 5) {
		I_ant = 0;
		e_ant = 0;
		PV_ant = 0;
		D_ant = 0;
		S_pid_ant = 0;
		MV_ant = 0;
	}

	// Compensation in the SP from DCS in case of closed loop estimation
	if (estimation)
		SP = (SP - 50) * 20;

	double e = SP - PV;
	
	// Proportional
	double P = Kp * e;

	// Integral
	double I = I_ant + Ki * 1 / 2 * (e_ant + e) * Ts;

	// Derivative
	double D = (-Kd * (PV - PV_ant) + Tdf * D_ant) / (Tdf + Ts);

	// Controller output
	double S_pid = P + I + D;

	// Output clipping
	double MV;
	if (S_pid > MV_max) {
		MV = MV_max;
		I = I_ant;
	}
	else if (S_pid < MV_min) {
		MV = MV_min;
		I = I_ant;
	}
	else
		MV = S_pid;

	excitation_cl(MV, PV, ct);
	
	// Delayed variables update
	I_ant = I;
	e_ant = e;
	PV_ant = PV;
	D_ant = D;
	S_pid_ant = S_pid;
	MV_ant = MV;

	return MV;
}


void Controller::excitation_cl(double& MV, double PV, int ct) {
	if (ct < 5) {
		flag_start_excitation = false;
		buffer_PV.clear();
		buffer_PV.reserve(1 / Ts);
		for (int i = 0; i < 1 / Ts; ++i)
			buffer_PV.push_back(0.0);
	}
	else {
		for (int i = buffer_PV.size()-1; i > 0 ; --i)
			buffer_PV[i] = buffer_PV[i-1];
		buffer_PV[0] = PV;
	}
	if (estimation && int(ct * Ts) == 100) {
		flag_friction = true;
		flag_calcStop = true;
		time = 0;
	}
	if (flag_friction) {
		if (flag_calcStop) {
			time += Ts;
			MV = MV_ant;
			if (time > 5 * tau_ip + 1) {
				std_PV = std_vec(buffer_PV);
				x_stop = buffer_PV;
				flag_calcStop = false;
				flag_goUp = true;
				time = 0;
			}
		}
		if (flag_goUp) {
			MV = MV_ant + Ts;
			double teste = min_vec(subtract_vect_const(x_stop, PV));
			if (max_vec(subtract_vect_const(x_stop, PV)) < - 2.0 * std_PV) {
				flag_goUp = false;
				flag_wait = true;
			}
		}
		if (flag_wait) {
			MV = MV_ant;
			time += Ts;
			if (time > 5 * tau_ip + 1) {
				flag_wait = false;
				flag_friction = false;
				flag_start_excitation = true;
				time = 0;
				t_exc = ct * Ts;
			}
		}
	}
	if (flag_start_excitation)
		MV = MV + exc_cl[ct];
}


void Controller::set_controller_parameters(std::vector<double> input){
	Kp = input[0];
	Ki = input[1];
	Kd = input[2];
	Tdf = input[3];
	MV_max = input[4];
	MV_min = input[5];
}