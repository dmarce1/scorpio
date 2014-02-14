
#include "../defs.h"

#ifdef ROTATING_DISC
#include "rotating_disc.h"
#include <vector>
#include <fstream>
#include <stdlib.h>

RotatingDisc* RotatingDisc::root = NULL;

RotatingDisc::~RotatingDisc() {
	// TODO Auto-generated destructor stub
}

RotatingDisc::RotatingDisc() {
	if (root == NULL) {
		shadow_off();
		root = this;
	}
	ChildIndex i;
	step_cnt = 0;
	ostep_cnt = 0;
	set_output_off();
	x_in = 0.0;
	grid_dim = 1.5;
	kick_mode = 1;
	ang_mom_cons = true;
}

void RotatingDisc::set_output_frequency_by_time(Real t) {
	int_output_freq = -1;
	time_output_freq = t;
}

void RotatingDisc::reset_time() {
	this->ostep_cnt = 0;
	root->set_time(0.0);
	step_cnt = 0;
}

void RotatingDisc::set_output_frequency_by_step(int a) {
	assert( a >= 1);
	int_output_freq = a;
}

void RotatingDisc::set_x_in(double x) {
	x_in = x;
}

void RotatingDisc::set_grid_dim(double x) {
	grid_dim = x;
}

void RotatingDisc::set_kick_mode(int x) {
	kick_mode = x;
}

void RotatingDisc::set_ang_mom_cons(bool x) {
	ang_mom_cons = x;
}

bool RotatingDisc::get_ang_mom_cons() const {
	return ang_mom_cons;
}

double RotatingDisc::get_ft_radius() const {
	double x_in = get_x_in();
	double R_outer = R_OUTER;
	double R_0 = R_outer * 2.0 * x_in / (1.0 + x_in);
	return R_0;
}

double RotatingDisc::get_frame_omega() const {
	const double x_in = get_x_in();
	const double G = BIG_G;
	const double M_c = M_C;
	const double R_outer = R_OUTER;
	const double j_H = sqrt(2.0 * x_in / (1.0 + x_in));
	const double j_here = j_H * sqrt(G * M_c * R_outer);
	const double omega_R_0 = (G * M_c) * (G * M_c) / (j_here * j_here * j_here);

#ifdef NO_ROTATE
	return 0.0;
#endif
	return omega_R_0;
}

double RotatingDisc::get_period() const {
	const double x_in = get_x_in();
	const double G = BIG_G;
	const double M_c = M_C;
	const double R_outer = R_OUTER;
	const double j_H = sqrt(2.0 * x_in / (1.0 + x_in));
	const double j_here = j_H * sqrt(G * M_c * R_outer);
	const double omega_R_0 = (G * M_c) * (G * M_c) / (j_here * j_here * j_here);
	const double period = 2.0 * M_PI / omega_R_0;

	return period;
}

double RotatingDisc::get_rho_max() const {
	const double x_in = root->get_x_in();

	const double R_outer = R_OUTER;
	const double kappa = KAPPA;
	const double n = POLYTROPIC_N;

	const double R_0 = R_outer * 2.0 * x_in / (1.0 + x_in);
	const double j_H = sqrt(2.0 * x_in / (1.0 + x_in));
	const double x_max = R_0 / R_outer;
	const double C = 1.0 / (1.0 + x_in);
	const double H_max = 1.0 / sqrt(x_max * x_max) - 0.5 * (j_H / x_max) * (j_H / x_max) - C;
	const double rho_max = pow(H_max / ((n + 1) * kappa), n);

	return rho_max;
}

double RotatingDisc::get_density_floor() const {
	const double rho_max = root->get_rho_max();

	return 1e-10 * rho_max;
}

double RotatingDisc::get_min_refine_density() const {
	const double rho_max = root->get_rho_max();

	return 1e-6 * rho_max;
}

double RotatingDisc::get_x_in() const {
	return x_in;
}

double RotatingDisc::get_grid_dim() const {
	return grid_dim;
}

int RotatingDisc::get_kick_mode() const {
	return kick_mode;
}

double RotatingDisc::convert_theta(double theta, double time) {
	double period = get_period();
	double theta_converted = theta - 2.0 * M_PI * (time / period);
	while (theta_converted < 0.0)
		theta_converted += 2.0 * M_PI;
	return theta_converted;

}

void RotatingDisc::ft_output(double time) {
	if (MPI_rank() == 0)
		printf("\nStarting ft_output...");
	double radius = get_ft_radius();
	std::vector<double> density;
	std::vector<double> theta;

	int time_int = nint(get_time() / time_output_freq);

	// Write filename to string.
	char filename[32];
	sprintf(filename, "ft.%i.dat", time_int);

	// Open file for output.
	std::ofstream outfile(filename);

	int err = root->fixed_radius_sweep(radius, density, theta, 0 // State index, 0 = density.
			);

	assert( density.size() == theta.size());

	for (int i = 0; i < density.size(); i++) {
		// Convert theta to intertial frame.
		double theta_here = convert_theta(theta[i], time);
		outfile << theta_here << " " << density[i] << "\n";
	}

	if (MPI_rank() == 0)
		printf("done.\n");

}

void RotatingDisc::set_refine_flags() {
	ChildIndex c;
	if (get_level() < 1) {
		for (int i = 0; i < OCT_NCHILD; i++) {
			set_refine_flag(i, true);
		}
	} else if (get_level() < get_max_level_allowed()) {
		for (int k = BW; k < GNX - BW; k++) {
			for (int j = BW; j < GNX - BW; j++) {
				for (int i = BW; i < GNX - BW; i++) {
					c.set_x(2 * i / GNX);
					c.set_y(2 * j / GNX);
					c.set_z(2 * k / GNX);
					bool test;

					if (true) {
						/* Zach - this is a modifcation I made to keep from refining all z zones */

						double R1 = sqrt(xc(i) * xc(i) + yc(j) * yc(j));
						double R2 = (R_OUTER + root->get_x_in()) / 2.0;
						double radius = sqrt(R1 * R1 + R2 * R2 - 2.0 * R1 * R2 + (zc(k) * zc(k)) *0.7); // 0.7 is fudge factor for vertical height
						test = radius < (R_OUTER - root->get_x_in()) / 2.0;

					} else /*if(false)*/{
						/* This is your original */
						double radius = sqrt(xc(i) * xc(i) + yc(j) * yc(j));
						if (get_level() > get_max_level_allowed() - 2) {
							test = ((radius < root->get_ft_radius()) && (radius > root->get_x_in() * R_OUTER));
						} else {
							test = ((radius < R_OUTER) && (radius > root->get_x_in() * R_OUTER));
						}

//                                printf ("max_refine_level = %i\n"),HydroGrid::max_refine_level;
//                                printf ("max_level = %i\n",max_level());

						/*                                    test = (
						 (radius < R_OUTER) &&
						 (radius > root->get_x_in()*R_OUTER)
						 );
						 */
					}
					set_refine_flag(c, test || get_refine_flag(c));
				}
			}
		}
	}
}

void RotatingDisc::ft_2D_output(double time) {
	if (MPI_rank() == 0)
		printf("\nStarting ft_2D_output...");
	const double x_in = get_x_in();
	const double R_outer = R_OUTER;

	int time_int = nint(get_time() / time_output_freq);

	int n_radial = 32;

	double r_out = 1.0 * R_outer;
	double r_in = x_in * R_outer;

	double dr = (r_out - r_in) / (n_radial - 1.0);

	// Write filename to string.
	char filename[32];
	sprintf(filename, "ft2D.%i.dat", time_int);

	// Open file for output.
	std::ofstream outfile(filename);

	// Loop over radii.
	for (int j = 0; j < n_radial; j++) {
		double radius = dr * j + r_in;

		std::vector<double> density;
		std::vector<double> theta;

		int err = root->fixed_radius_sweep(radius, density, theta, 0 // State index, 0 = density.
				);

		assert( density.size() == theta.size());

		for (int i = 0; i < density.size(); i++) {
			// Convert theta to intertial frame.
			double theta_here = convert_theta(theta[i], time);
			outfile << radius << " " << theta_here << " " << density[i] << "\n";
		}

	}
	if (MPI_rank() == 0)

		printf("done.\n");

}

void RotatingDisc::step_to_time(Real tmax) {
	root->set_time(0.0);
	bool last_step = false;
	bool do_output;
	Real tleft, dt, last_dt;
	Real next_output_time;
	dt = 0.0;
	Vector<State, 4> state_sum;
	int iters = 0;
	iters = 0;
	/*******/
	setup_grid_structure();
	dt = next_dt(&do_output, &last_step, &ostep_cnt);
	/*********/
	dt *= min(1.0, MAXINITDT);
	int last_max_level;
	int last_max_level2;
	if (MPI_rank() == 0)

		printf("opening file...\n");
	FILE* fp_ts = fopen("ts.dat", "w");
	if (MPI_rank() == 0)
		printf("done opening file\n");

#ifdef OUTPUT_NONFLAT
	root->output_nonflat("X", 0.0);
#endif
#ifndef OUTPUT_NONFLAT
	output("X", 0.0, GNX, BW);
#endif

//	ft_output(0.0);
//	ft_2D_output(0.0);

	do {
		dt = max_dt_driver();
		double dt_cfl = dt;
		// Output immediately if max_dt is too small.
		if (dt < 1e-15) {
#ifdef OUTPUT_NONFLAT
			root->output_nonflat("laststep", 0.0);
#endif
#ifndef OUTPUT_NONFLAT
			root->output("laststep", 0.0, GNX, BW);
#endif
			exit(1);
		}
		if (get_time() == 0.0) {
			dt = MAXINITDT * dt;
		} else if (step_cnt != 0) {
			dt = min(last_dt * MAXDTINC, dt);
		}
		last_dt = dt;
		do_output = false;
		tleft = tmax - get_time();
		if (dt + get_time() >= next_output_time) {
			next_output_time = Real(ostep_cnt + 1) * time_output_freq;
//                dt = next_output_time - get_time();
			ostep_cnt++;
			do_output = true;
		}
		/*		if (int_output_freq < 0) {
		 next_output_time = Real(ostep_cnt + 1) * time_output_freq;
		 if (dt + get_time() >= next_output_time) {
		 dt = next_output_time - get_time();
		 ostep_cnt++;
		 do_output = true;
		 } else {
		 dt = min(dt,(next_output_time - get_time()) / int(
		 (next_output_time - get_time()) / dt + 1.0));
		 }
		 } else if (int_output_freq != 0) {
		 do_output = bool(step_cnt % int_output_freq == 0);
		 }
		 */
		//		printf( "%e\n", dt );
		if (tleft <= dt) {
			dt = tleft;
			last_step = true;
			if (int_output_freq < 0 && next_output_time == tmax) {
				do_output = true;
			}
		}

		this->step(dt);
		step_cnt++;
		if (MPI_rank() == 0)

			printf("step=%i t=%e dt=%e lmax=%i ngrids=%i \n", step_cnt, HydroGrid::get_time(), dt, OctNode::get_max_level(), get_node_cnt());
		if (do_output) {

#ifdef OUTPUT_NONFLAT
			root->output_nonflat("X", nint(get_time() / time_output_freq));
#endif
#ifndef OUTPUT_NONFLAT
			root->output("X", nint(get_time() / time_output_freq), GNX, BW);
#endif
			// Call function to write out neccesary data for FT.
			//ft_output(nint(get_time()/time_output_freq));
			//		ft_output(get_time());
			//		ft_2D_output(get_time());

			//	printf("\n --- output %i ---", nint(get_time() / time_output_freq));
		}
	} while (!last_step);
	fclose(fp_ts);
	fp_ts = NULL;
}

void RotatingDisc::set_output_off() {
	int_output_freq = 0;
}

int RotatingDisc::fixed_radius_sweep(double r, std::vector<double>& state, std::vector<double>& theta,
//                                 double state,
//                                 double theta,
		int state_index) {
	int k, j, i;
	double dx = get_dx();
	double z_here = 1e-7;

	// Make sure we are near the equitorial plane.
	if (!((zc(BW) - 0.5 * dx < z_here) && (zc(GNX - BW - 1) + 0.5 * dx >= z_here)))
		return 0;
	// Now determine our z index, k.
	k = -1;
	for (i = BW; i < GNX - BW; i++)
		if ((z_here >= zc(i) - 0.5 * dx) && (z_here < zc(i) + 0.5 * dx))
			k = i;

//    printf("k = %i\n",k);

	// Loop over all x's in this grid node object.
	for (i = BW; i < GNX - BW; i++) {
		double x_here = xc(i);
//        printf("x_here = %f r = %f\n", x_here, r);

//        double y_here_squared = r*r - x_here*x_here;
//        double y_here;

		// Loop over y values.
		// There may be more than one point for each x.
		// Or there may be none.
		for (int j = BW; j < GNX - BW; j++) {
			double y_squared_here = r * r - x_here * x_here;

			if (y_squared_here > 0.0) {
				double y_here = sqrt(y_squared_here);
				if (((y_here >= yc(j) - 0.5 * dx) && (y_here < yc(j) + 0.5 * dx)) || ((-y_here >= yc(j) - 0.5 * dx) && (-y_here < yc(j) + 0.5 * dx))) {
					if (!zone_is_refined(i, j, k)) {
						if (yc(j) < 0.0)
							y_here = -y_here;

						int j_L = -1;
						int j_U = -1;
						if (y_here > yc(j))
							j_L = j;
						if (y_here <= yc(j))
							j_L = j - 1;
						j_U = j_L + 1;
						assert(j_L > -1);
						assert(j_U > -1);
//                printf("i = %i j = %i k = %i\n",i,j_U,k);
						double state_upper = (*this)(i, j_U, k)[state_index];
						double state_lower = (*this)(i, j_L, k)[state_index];
						double y_upper = yc(j_U);
						double y_lower = yc(j_L);

						double theta_here = atan2(y_here, x_here);
						theta_here += M_PI;

						double state_here = interpolate(state_upper, state_lower, y_upper, y_lower, y_here);
//                        state.push_back(yc(j));
//                        theta.push_back(x_here);
						double radius_here = sqrt(x_here * x_here + y_here * y_here);
						state.push_back(state_here);
						theta.push_back(theta_here);
					}
				}

			}

		}
	}

	// Now to the reverse process
	for (i = BW; i < GNX - BW; i++) {
		double y_here = yc(i);

		// Loop over x values.
		// There may be more than one point for each y.
		// Or there may be none.
		for (int j = BW; j < GNX - BW; j++) {
			double x_squared_here = r * r - y_here * y_here;

			if (x_squared_here > 0.0) {
				double x_here = sqrt(x_squared_here);
				if (((x_here >= xc(j) - 0.5 * dx) && (x_here < xc(j) + 0.5 * dx)) || ((-x_here >= xc(j) - 0.5 * dx) && (-x_here < xc(j) + 0.5 * dx))) {
					if (!zone_is_refined(i, j, k)) {
						if (xc(j) < 0.0)
							x_here = -x_here;

						int j_L = -1;
						int j_U = -1;
						if (x_here > xc(j))
							j_L = j;
						if (x_here <= xc(j))
							j_L = j - 1;
						j_U = j_L + 1;
						assert(j_L > -1);
						assert(j_U > -1);
//                printf("i = %i j = %i k = %i\n",i,j_U,k);
						double state_upper = (*this)(j_U, i, k)[state_index];
						double state_lower = (*this)(j_L, i, k)[state_index];
						double x_upper = xc(j_U);
						double x_lower = xc(j_L);

						double theta_here = atan2(y_here, x_here);
						theta_here += M_PI;

						double state_here = interpolate(state_upper, state_lower, x_upper, x_lower, x_here);
//                        state.push_back(xc(j));
//                        theta.push_back(x_here);
						double radius_here = sqrt(x_here * x_here + y_here * y_here);
						state.push_back(state_here);
						theta.push_back(theta_here);
					}
				}

			}

		}
	}

	// This block of code recurses to all the children.
	for (int i = 0; i < OCT_NCHILD; i++) {
		if (this->get_child(i) != NULL) {
			int err = dynamic_cast<RotatingDisc*>(get_child(i))->fixed_radius_sweep(r, state, theta, state_index);
		}

	}

	return 0;

}

#endif
