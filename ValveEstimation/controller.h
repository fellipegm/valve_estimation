#pragma once
#ifndef _CONTROLLER_H
#define _CONTROLLER_H

#include <vector>

class Controller {
	private:

		// Sampling time
		double Ts{ 1e-3 };

		// PID controller parameters
		double Kp{ 0.1 };
		double Ki{ 0.0 };
		double Kd{ 0.0 };
		double Tdf{ 0.0 };
		double MV_max{ 100.0 };
		double MV_min{ -100.0 };

		double tau_ip;

		// PID internal parameters
		double I_ant{ 0.0 }, e_ant{ 0.0 }, PV_ant{ 0.0 }, D_ant{ 0.0 }, S_pid_ant{ 0.0 }, MV_ant{ 0.0 };

		bool estimation = false;

		std::vector<double> exc_cl;

		// Algorithm to start excitation signal
		bool flag_calcStop{ false }, flag_goUp{ false };
		bool flag_goDown{ false }, flag_wait{ false };
		bool flag_friction{ false }, flag_start_excitation{ false };
		double std_PV{ 0.0 };
		double time{ 0.0 };
		void excitation_cl(double& MV, double PV, int ct);
		std::vector<double> buffer_PV, x_stop;

	public:
		// Time in which the excitation effectively started
		double t_exc{ 0.0 };


		double pid(double SP, double PV, int ct);
		
		void set_controller_parameters(std::vector<double> input);
		void set_excitation(std::vector<double> input) { exc_cl = input; };
		void set_estimation(bool input) { estimation = input; };
		void set_tau_ip(double input) { tau_ip = input; };
};



#endif