#ifndef ROTATING_DISC_H_
#define ROTATING_DISC_H_

#include "../defs.h"

#ifdef ROTATING_DISC

#include "../hydro_grid/hydro_grid.h"
#include <vector>

class RotatingDisc: public HydroGrid {
private:
	double x_in;
	double grid_dim;
	int kick_mode;
	bool ang_mom_cons;
	int int_output_freq;
	int ostep_cnt;
	int verbosity;
	Real time_output_freq;
	static RotatingDisc* root;
	virtual void set_refine_flags();
protected:
	int step_cnt;
public:
	void reset_time();
	void set_output_frequency_by_time(Real);
	void set_output_frequency_by_step(int);
	void set_output_off();
	void set_verbosity(int);
	void step_to_time(Real);
	void set_x_in(double);
	void set_grid_dim(double);
	void set_kick_mode(int);
	void set_ang_mom_cons(bool);
	double get_x_in() const;
	double get_ft_radius() const;
	double get_min_refine_density() const;
	double get_density_floor() const;
	double get_rho_max() const;
	double get_frame_omega() const;
	double get_period() const;
	double get_grid_dim() const;
	bool get_ang_mom_cons() const;
	int get_kick_mode() const;
	void ft_output(double);
	void ft_2D_output(double);
	double convert_theta(double, double);
public:
	RotatingDisc();
	virtual ~RotatingDisc();
	virtual RotatingDisc* new_octnode() const {
		return new RotatingDisc;
	}
	virtual void initialize() {
		RotatingDisc* g = this;
		State U = Vector<Real, STATE_NF>(0.0);
		const Real rho0 = 1.0;
		const Real ei0 = 1.0;
		const Real tau0 = pow(ei0, 1.0 / State::gamma);
		const Real rho1 = root->get_density_floor();
		const Real ei1 = root->get_density_floor();
		const Real tau1 = pow(ei1, 1.0 / State::gamma);

		const double x_in = root->get_x_in();
		const double kappa = KAPPA;
		const double G = BIG_G;
		const double M_c = M_C;
		const double R_outer = R_OUTER;
		const double n = POLYTROPIC_N;
		const double a_0 = A_0;

		const double R_inner = x_in * R_outer;

		const double C = 1.0 / (1.0 + x_in);
		const double j_H = sqrt(2.0 * x_in / (1.0 + x_in));
		double j_here = j_H * sqrt(G * M_c * R_outer); // conversion to
													   // 'real' units

		const double R_0 = R_outer * 2.0 * x_in / (1.0 + x_in);

		const double period = root->get_period();
		const double ft_radius = root->get_ft_radius();
		const double rho_max = root->get_rho_max();

		const double omega_R_0 = root->get_frame_omega();

		const double kick_amplitude = 1e-2;
		const int kick_mode = root->get_kick_mode();

		const double grid_dim = GRID_DIM;

		if (g->get_level() == 0) {
			printf("------------------------\n");
			printf("Physics parameters:\n");
			printf("K = %10.3e\n", KAPPA);
			printf("M_C = %10.3e\n", M_c);
			printf("G = %10.3e\n", G);
			printf("EULER_GAMMA = %10.3e\n", EULER_GAMMA);
			printf("grid_dim = %10.3e\n", grid_dim);
			printf("levels of refinement = %i\n", this->get_max_level_allowed());

			printf("\nProblem parameters:\n");
			printf("x_in = %10.3e\n", x_in);
			printf("R_outer=%10.3e \n", R_outer);
			printf("j_H=%10.3e \n", j_H);
			printf("C=%10.3e \n", C);
			printf("R_0 = %10.3e \n", R_0);
			printf("omega = %10.3e\n", omega_R_0);
			printf("period = %10.3e\n", period);
			//            printf("OUTPUT_TIME_FREQ = %10.3e\n",OUTPUT_TIME_FREQ);
			printf("R_inner=%10.3e \n", R_inner);
			//    printf("density floor = %10.3e\n",DENSITY_FLOOR);
			printf("maximum initial density = %10.3e\n", rho_max);
			printf("------------------------\n");
		}

		for (int k = 0; k < GNX; k++) {
			for (int j = 0; j < GNX; j++) {
				for (int i = 0; i < GNX; i++) {
					const Real x_here = g->xc(i);
					const Real y_here = g->yc(j);
					const Real z_here = abs(g->zc(k));

					const double r = sqrt(pow(x_here, 2.0) + pow(y_here, 2.0)); //cylindrical R
					const double x = r / R_outer;

					// 3d-torus
					if ((R_inner <= r) && (R_outer >= r)) {
						const Real temporary = pow(C + 0.5 * (j_H / x) * (j_H / x), -2.0) - x * x;
						assert( temporary > 0.0);
						// z_max has to be converted to "real"
						// units first.
						const Real z_max = R_outer * sqrt(temporary);
						if (z_here <= z_max) {
							const double z = z_here / R_outer;
							const double H_here = 1.0 / sqrt(x * x + z * z) - 0.5 * (j_H / x) * (j_H / x) - C;
							double rho_here = pow(H_here / ((n + 1) * kappa), n);
							rho_here *= G * M_c / R_outer; // convert
							assert( rho_here > 0.0);
							// Perturbation:
							//                              double g=((double)rand()/(double)RAND_MAX);
							//                              g *= 2.0;
							//                              g -= 1.0;
							// random perturbation
							//                              rho_here *= 1.0+a_0*g;

							double theta_here = atan2(y_here, x_here);
							//mode specific "kick"
							//                              printf("kick_amplitude*cos(kick_mode*theta_here) = %10.3f\n",kick_amplitude*cos(kick_mode*theta_here));

							rho_here *= 1.0 + kick_amplitude * cos(kick_mode * theta_here);
							assert( rho_here > 0.0);

							// Convert from dimensionless units:

							U.set_rho(rho_here);
							//                                  U.set_rho(1.0);
							U.set_et(ei0);
							U.set_tau(tau0);
				//			U.set_sx(-y_here * rho_here * j_here / pow(r, 2.0));
				//			U.set_sy(x_here * rho_here * j_here / pow(r, 2.0));
							U.set_sx(j_here * rho_here);
							U.set_sy(0.0);
						} else {
							U.set_rho(rho1);
							U.set_et(ei1);
							U.set_tau(tau1);
							U.set_sx(0.0);
		//					U.set_sx(0.0);
			//				U.set_sy(0.0);
							U.set_sy(0.0);
						}
					} else {
						U.set_rho(rho1);
						U.set_et(ei1);
						U.set_tau(tau1);
						U.set_sx(0.0);
						U.set_sy(0.0);
						U.set_sy(0.0);
					}
					U.set_sz(0.0);

					// Test blob
					// if ( (abs(x_here-1e-4) < 1e-5) && (abs(y_here-1e-4) < 1e-5) && (abs(z_here) < 1e-5))
					// {
					//     U.set_rho(1.0);
					//     U.set_et(1.0);
					//     U.set_tau(1.0);
					//     U.set_sx(0.0);
					//     U.set_sy(0.0);
					//     U.set_lz(0.0);
					//     U.set_sz(0.0);
					// }
					(*g)(i, j, k) = U;
				}
			}
		}
	}

	int fixed_radius_sweep(double r, std::vector<double>& state, std::vector<double>& theta,
	//                                 double state,
	//                                 double theta,
			int state_index);

	double interpolate(double su, double sl, double yu, double yl, double y) {
		return sl + (su - sl) * (y - yl) / (yu - yl);
	}

};

#endif
#endif /* ROTATING_DISC_H_ */
