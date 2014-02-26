#include "binary_star.h"
#include "../FMM/FMM.h"
#include <list>

#ifdef HYDRO_GRAV_GRID
#ifdef BINARY_STAR

#define CHECKPT_FREQ 250

//static Real MA = 1.04;
//static Real MD = 0.20;
static Real MA = 1.04;
static Real MD = 0.20;
Real BinaryStar::refine_floor = 1.0e-4;
bool BinaryStar::scf_code = false;
Real BinaryStar::q_ratio = MD / MA;
Real BinaryStar::polyK;
Real BinaryStar::M1 = 5.0e-6, BinaryStar::M2 = 2.5e-6;
Real BinaryStar::Kd, BinaryStar::Ka;
Real BinaryStar::Ax, BinaryStar::Bx, BinaryStar::Cx, BinaryStar::Aphi, BinaryStar::Bphi, BinaryStar::Cphi, BinaryStar::rho_max_a, BinaryStar::rho_max_d;
_3Vec BinaryStar::a0, BinaryStar::d0, BinaryStar::com_vel_correction = 0.0;
Real BinaryStar::lz_t0 = 0.0;
Real BinaryStar::code_to_cm, BinaryStar::code_to_s, BinaryStar::code_to_K, BinaryStar::code_to_g;

void find_eigenvectors(double q[3][3], double e[3][3], double lambda[3]) {
    double b0[3], b1[3], A, bdif;
    int iter = 0;
    for (int l = 0; l < 3; l++) {
        b0[0] = b0[1] = b0[2] = 0.0;
        b0[l] = 1.0;
        do {
            iter++;
            for (int i = 0; i < 3; i++) {
                b1[i] = 0.0;
            }
            for (int i = 0; i < 3; i++) {
                for (int m = 0; m < 3; m++) {
                    b1[i] += q[i][m] * b0[m];
                }
            }
            A = 0.0;
            for (int i = 0; i < 3; i++) {
                A += b1[i] * b1[i];
            }
            A = sqrt(A);
            bdif = 0.0;
            for (int i = 0; i < 3; i++) {
                b1[i] = b1[i] / A;
                bdif += pow(b0[i] - b1[i], 2);
            }
            for (int i = 0; i < 3; i++) {
                b0[i] = b1[i];
            }

        } while (fabs(bdif) > 1.0e-10);
        for (int m = 0; m < 3; m++) {
            e[l][m] = b0[m];
        }
        for (int i = 0; i < 3; i++) {
            A += b0[i] * q[l][i];
        }
        lambda[l] = sqrt(A) / sqrt(e[l][0] * e[l][0] + e[l][1] * e[l][1] + e[l][2] * e[l][2]);
    }
}

void BinaryStar::compute_axis(Real dt, Real* theta, Real* theta_dot) {

    BinaryStar* g;
    Real Q[3][3] = { { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 } }, e[3][3], lambda[3];
    Real Q0[3][3] = { { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 } }, e0[3][3], lambda0[3];
    Real outbuf[6], inbuf[6];
    int cnt;
    double tstart = MPI_Wtime();
    for (int l = 0; l < get_local_node_cnt(); l++) {
        g = dynamic_cast<BinaryStar*>(get_local_node(l));
        for (int k = BW; k < GNX - BW; k++) {
            for (int j = BW; j < GNX - BW; j++) {
                for (int i = BW; i < GNX - BW; i++) {
                    if (g->zone_is_refined(i, j, k)) {
                        continue;
                    }
                    _3Vec x = g->X(i, j, k);
                    Real r2 = x.dot(x);
                    Real dm = g->U(i, j, k).rho() * pow(g->get_dx(), 3);
                    Real dm0 = g->U0(i, j, k).rho() * pow(g->get_dx(), 3);
                    for (int n = 0; n < 3; n++) {
                        for (int m = n; m < 3; m++) {
                            Q[n][m] += 3.0 * x[n] * x[m] * dm;
                            Q0[n][m] += 3.0 * x[n] * x[m] * dm0;
                        }
                    }
                }
            }
        }
    }
    cnt = 0;
    for (int n = 0; n < 3; n++) {
        for (int m = n; m < 3; m++) {
            outbuf[cnt] = Q[n][m];
            cnt++;
        }
    }
    MPI_Allreduce(outbuf, inbuf, 6, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
    cnt = 0;
    for (int n = 0; n < 3; n++) {
        for (int m = n; m < 3; m++) {
            Q[n][m] = Q[m][n] = inbuf[cnt];
            cnt++;
        }
    }
    find_eigenvectors(Q, e, lambda);
    cnt = 0;
    for (int n = 0; n < 3; n++) {
        for (int m = n; m < 3; m++) {
            outbuf[cnt] = Q0[n][m];
            cnt++;
        }
    }
    MPI_Allreduce(outbuf, inbuf, 6, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
    cnt = 0;
    for (int n = 0; n < 3; n++) {
        for (int m = n; m < 3; m++) {
            Q0[n][m] = Q0[m][n] = inbuf[cnt];
            cnt++;
        }
    }
    find_eigenvectors(Q0, e0, lambda0);
    int ei;
    if (lambda[0] > lambda[1] && lambda[0] > lambda[2]) {
        ei = 0;
    } else if (lambda[1] > lambda[0] && lambda[1] > lambda[2]) {
        ei = 1;
    } else {
        ei = 2;
    }
    *theta = asin(e[ei][1]);
    *theta_dot = (asin(e[ei][1]) - asin(e0[ei][1])) / dt;
    Real w = 10.0 * State::omega;
    State::omega_dot = w * (w * *theta + 2.0 * *theta_dot);
    //  State::omega_dot = 0.0;
    Real omega0 = State::omega;
    Real omega1 = omega0 + dt * State::omega_dot;
    if (MPI_rank() == 0) {
        FILE* fp = fopen("eigen.dat", "at");
        fprintf(fp, "%.6e %.12e %.12e %.12e %.12e %.12e   \n", dynamic_cast<BinaryStar*>(get_root())->get_time(), State::omega, State::omega_dot, *theta, *theta_dot, dtheta);
        fclose(fp);
    }
    for (int l = 0; l < get_local_node_cnt(); l++) {
        g = dynamic_cast<BinaryStar*>(get_local_node(l));
        for (int k = 0; k < GNX; k++) {
            for (int j = 0; j < GNX; j++) {
                for (int i = 0; i < GNX; i++) {
                    _3Vec x = g->X(i, j, k);
                    Real r2 = x[0] * x[0] + x[1] * x[1];
                    State u = g->U(i, j, k);
                    g->U(i, j, k)[State::et_index] -= ((omega1 - omega0) * u.lz() - 0.5 * u[State::d_index] * r2 * ((omega1 * omega1) - (omega0 * omega0)));
                }
            }
        }
    }
    State::omega = omega1;

}

void BinaryStar::adjust_frame(Real dt) {
    _3Vec com, com_dot;
    Real w;
    Real theta, theta_dot;
    com_dot = compute_Mdot(&com);
    w = 10.0 * State::omega;
    Real x, y, z, xdot, ydot, zdot, r, a, tmp, ax, ay, az;
    x = com[0];
    y = com[1];
    z = com[2];
    xdot = com_dot[0];
    ydot = com_dot[1];
    zdot = com_dot[2];
    r = sqrt(x * x + y * y + z * z);
    tmp = (x * xdot + y * ydot + z * zdot) / r;
    ax = -w * (w * r + 2.0 * tmp) * (x / r);
    ay = -w * (w * r + 2.0 * tmp) * (y / r);
    az = -w * (w * r + 2.0 * tmp) * (z / r);

    for (int l = 0; l < get_local_node_cnt(); l++) {
        BinaryStar* g = dynamic_cast<BinaryStar*>(get_local_node(l));
        Real dv = pow(g->get_dx(), 3);
        for (int k = BW; k < GNX - BW; k++) {
            for (int j = BW; j < GNX - BW; j++) {
                for (int i = BW; i < GNX - BW; i++) {
                    State& u = g->U(i, j, k);
                    _3Vec x = g->X(i, j, k);
                    Real R2 = x[0] * x[0] + x[1] * x[1];
                    u.set_sx(u.sx() + u.rho() * ax * dt);
                    u.set_sy(u.sy() + u.rho() * ay * dt);
                    u.set_sz(u.sz() + u.rho() * az * dt);
                    u.set_lz(u.lz() + u.rho() * (x[0] * ay - x[1] * ax) * dt);
                }
            }
        }
    }
    if (MPI_rank() == 0) {
        FILE* fp = fopen("adjust.dat", "at");
        fprintf(fp, "%e %e %e %e %e %e %e %e %e %e \n", dynamic_cast<BinaryStar*>(get_root())->get_time(), com[0], com[1], com[2], com_dot[0], com_dot[1],
                com_dot[2], ax, ay, az);
        fclose(fp);
    }

}

Real BinaryStar::compute_I() {
    BinaryStar* g;
    Real I, tmp;
    I = 0.0;
    for (int l = 0; l < get_local_node_cnt(); l++) {
        g = dynamic_cast<BinaryStar*>(get_local_node(l));
        g->euler_force_coeff.allocate();
        Real dv = pow(g->get_dx(), 3);
        for (int k = BW; k < GNX - BW; k++) {
            for (int j = BW; j < GNX - BW; j++) {
                for (int i = BW; i < GNX - BW; i++) {
                    State U = g->U(i, j, k);
                    _3Vec x = g->X(i, j, k);
                    Real R2 = x[0] * x[0] + x[1] * x[1];
                    g->euler_force_coeff(i, j, k) = -(U.lz() - R2 * State::omega * U.rho());
                    if (!g->zone_is_refined(i, j, k)) {
                        I += R2 * State::omega * U.rho() * dv;
                    }
                }
            }
        }
    }
    tmp = I;
    MPI_Allreduce(&tmp, &I, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
    return I;
}

_3Vec BinaryStar::compute_Mdot(_3Vec *com) {
    BinaryStar* g;
    _3Vec Mdot;
    Real M;
    Real a[4], b[4];
    Mdot = 0.0;
    *com = 0.0;
    M = 0.0;
    for (int l = 0; l < get_local_node_cnt(); l++) {
        g = dynamic_cast<BinaryStar*>(get_local_node(l));
        Real dv = pow(g->get_dx(), 3);
        for (int k = BW; k < GNX - BW; k++) {
            for (int j = BW; j < GNX - BW; j++) {
                for (int i = BW; i < GNX - BW; i++) {
                    if (g->zone_is_refined(i, j, k)) {
                        continue;
                    }
                    _3Vec x = g->X(i, j, k);
                    Mdot += x * (g->U(i, j, k).rho() - g->U0(i, j, k).rho()) * dv;
                    M += g->U(i, j, k).rho() * dv;
                    *com += x * g->U(i, j, k).rho() * dv;
                }
            }
        }
    }
    for (int i = 0; i < 3; i++) {
        a[i] = Mdot[i];
    }
    a[3] = M;
    MPI_Allreduce(a, b, 4, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
    for (int i = 0; i < 3; i++) {
        Mdot[i] = b[i];
    }
    M = b[3];
    for (int i = 0; i < 3; i++) {
        a[i] = (*com)[i];
    }
    MPI_Allreduce(a, b, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
    for (int i = 0; i < 3; i++) {
        (*com)[i] = b[i];
    }
    Mdot /= M;
    *com /= M;
    // printf( "%e %e %e\n", (*com)[0], (*com)[1], (*com)[2]);
    // Mdot -= com * State::omega * 10.0;
    return Mdot;
}

void BinaryStar::apply_Mdot(_3Vec mdot) {
    BinaryStar* g;
    Real fx, fy, fz;
    for (int l = 0; l < get_local_node_cnt(); l++) {
        g = dynamic_cast<BinaryStar*>(get_local_node(l));
        for (int k = BW; k < GNX - BW; k++) {
            for (int j = BW; j < GNX - BW; j++) {
                for (int i = BW; i < GNX - BW; i++) {
                    _3Vec x = g->X(i, j, k);
                    fx = g->U(i, j, k).rho() * mdot[0];
                    fy = g->U(i, j, k).rho() * mdot[1];
                    fz = g->U(i, j, k).rho() * mdot[2];
                    g->U(i, j, k)[State::sx_index] -= fx;
                    g->U(i, j, k)[State::sy_index] -= fy;
                    g->U(i, j, k)[State::sz_index] -= fz;
                    g->U(i, j, k)[State::lz_index] -= x[0] * fy - x[1] * fx;
                }
            }
        }
        g->euler_force_coeff.deallocate();
    }
    if (MPI_rank() == 0) {
        FILE* fp = fopen("mdot.dat", "at");
        fprintf(fp, "%e %e %e %e\n", dynamic_cast<BinaryStar*>(get_root())->get_time(), mdot[0], mdot[1], mdot[2]);
        fclose(fp);
    }
}
void BinaryStar::apply_omega_dot(Real odot, Real dt, Real beta) {
    BinaryStar* g;
    Real d;
    for (int l = 0; l < get_local_node_cnt(); l++) {
        g = dynamic_cast<BinaryStar*>(get_local_node(l));
        for (int k = BW; k < GNX - BW; k++) {
            for (int j = BW; j < GNX - BW; j++) {
                for (int i = BW; i < GNX - BW; i++) {
                    d = g->euler_force_coeff(i, j, k) * odot;
                    g->U(i, j, k)[State::et_index] += beta * dt * d;
                }
            }
        }
        g->euler_force_coeff.deallocate();
    }
}

void BinaryStar::assign_fracs(Real Hfrac, Real min_phi, Real max_phi) {
    BinaryStar* g;
    Real d, a, minc, maxc, midc, dm, tmp;
    minc = min_phi;
    maxc = max_phi;
    while (log(fabs(minc / maxc)) > 1.0e-6) {
        a = 0.0;
        d = 0.0;
        midc = 0.5 * (maxc + minc);
        for (int l = 0; l < get_local_node_cnt(); l++) {
            g = dynamic_cast<BinaryStar*>(get_local_node(l));
            Real dv = pow(g->get_dx(), 3);
            for (int k = BW; k < GNX - BW; k++) {
                for (int j = BW; j < GNX - BW; j++) {
                    for (int i = BW; i < GNX - BW; i++) {
                        if (!g->zone_is_refined(i, j, k)) {
                            dm = (*g)(i, j, k).rho() * dv;
                            if (g->X(i, j, k)[0] < 0.0 || (*g)(i, j, k).phi_eff() > midc) {
                                (*g)(i, j, k).set_frac(0, 0.0);
                                (*g)(i, j, k).set_frac(1, (*g)(i, j, k).rho());
                                d += dm;
                            } else {
                                (*g)(i, j, k).set_frac(0, (*g)(i, j, k).rho());
                                (*g)(i, j, k).set_frac(1, 0.0);
                                a += dm;
                            }
                        }
                    }
                }
            }
        }
        tmp = a;
        MPI_Allreduce(&tmp, &a, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
        tmp = d;
        MPI_Allreduce(&tmp, &d, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
        if (d / (a + d) > Hfrac) {
            minc = midc;
        } else {
            maxc = midc;
        }
        if (MPI_rank() == 0) {
            printf("%e %e %e %e\n", minc, midc, maxc, d / (a + d));
        }
    }
}

Real BinaryStar::virial_error() {
    BinaryStar *g;
    Real tmp, w, t, ek_tot, p_tot, w_tot1, w_tot2, w_tot3;
    ek_tot = p_tot = w_tot1 = w_tot2 = w_tot3 = 0.0;
    for (int l = 0; l < get_local_node_cnt(); l++) {
        g = dynamic_cast<BinaryStar*>(get_local_node(l));
        Real dv = pow(g->get_dx(), 3);
        for (int k = BW; k < GNX - BW; k++) {
            for (int j = BW; j < GNX - BW; j++) {
                for (int i = BW; i < GNX - BW; i++) {
                    if (!g->zone_is_refined(i, j, k)) {
                        _3Vec X = g->X(i, j, k);
                        ek_tot += -(*g)(i, j, k).rot_pot(X) * dv;
                        p_tot += (*g)(i, j, k).pd() * dv;
                        //					w_tot1 += g->get_phi(i, j, k) * (*g)(i, j, k).rho() * 0.5 * dv;
                        w_tot2 += g->gx(i, j, k) * g->HydroGrid::xc(i) * (*g)(i, j, k).rho() * dv;
                        w_tot2 += g->gy(i, j, k) * g->HydroGrid::yc(j) * (*g)(i, j, k).rho() * dv;
                        w_tot2 += g->gz(i, j, k) * g->HydroGrid::zc(k) * (*g)(i, j, k).rho() * dv;
                        //				w_tot3 -= g->gz(i, j, k) * g->gz(i, j, k) * 0.5 * dv / (4.0 * M_PI);
                        //				w_tot3 -= g->gy(i, j, k) * g->gy(i, j, k) * 0.5 * dv / (4.0 * M_PI);
                        //				w_tot3 -= g->gx(i, j, k) * g->gx(i, j, k) * 0.5 * dv / (4.0 * M_PI);
                    }
                }
            }
        }
    }
    tmp = p_tot;
    MPI_Allreduce(&tmp, &p_tot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
//	tmp = w_tot1;
//	MPI_Allreduce(&tmp, &w_tot1, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
    tmp = w_tot2;
    MPI_Allreduce(&tmp, &w_tot2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
//	tmp = w_tot3;
//	MPI_Allreduce(&tmp, &w_tot3, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
    tmp = ek_tot;
    MPI_Allreduce(&tmp, &ek_tot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
    /*	printf("%e %e %e %e %e %e %e %e\n", 2.0 * ek_tot, 3.0 * p_tot, w_tot1, w_tot2, w_tot3, (2.0 * ek_tot + 3.0 * p_tot + w_tot1) / w_tot1,
     (2.0 * ek_tot + 3.0 * p_tot + w_tot2) / w_tot2, (2.0 * ek_tot + 3.0 * p_tot + w_tot3) / w_tot3);*/
    return (2.0 * ek_tot + 3.0 * p_tot + w_tot2) / w_tot2;
}

void BinaryStar::find_rho_max(Real* rho1, Real* rho2) {
    BinaryStar *g;
    Real tmp;
    *rho1 = 0.0;
    *rho2 = 0.0;
    for (int l = 0; l < get_local_node_cnt(); l++) {
        g = dynamic_cast<BinaryStar*>(get_local_node(l));
        for (int k = BW; k < GNX - BW; k++) {
            for (int j = BW; j < GNX - BW; j++) {
                for (int i = BW; i < GNX - BW; i++) {
                    if (!g->zone_is_refined(i, j, k)) {
                        *rho1 = max(*rho1, (*g)(i, j, k).frac(0));
                        *rho2 = max(*rho2, (*g)(i, j, k).frac(1));
                    }
                }
            }
        }
    }
    tmp = *rho1;
    MPI_Allreduce(&tmp, rho1, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );
    tmp = *rho2;
    MPI_Allreduce(&tmp, rho2, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );
}

#ifdef ZTWD
Real BinaryStar::find_K(int frac, Real phi0, Real center, Real l1_x, Real* span) {
    BinaryStar* g;
    const Real n = 1.5;
    Real K, M, tmp, dv, xmin, xmax, A, df, M0, f;
    _3Vec x0, dx, x;
    bool test;
    x0 = 0.0;
    x0[0] = center;
    if (frac == 0) {
        M0 = M1;
        A = Ka;
    } else {
        M0 = M2;
        A = Kd;
    }
    do {
        xmin = 1.0e+99;
        xmax = -1.0e-99;
        Real dphi;
        M = 0.0;
        df = 0.0;
        for (int l = 0; l < get_local_node_cnt(); l++) {
            g = dynamic_cast<BinaryStar*>(get_local_node(l));
            dv = pow(g->get_dx(), 3);
            for (int k = BW; k < GNX - BW; k++) {
                for (int j = BW; j < GNX - BW; j++) {
                    for (int i = BW; i < GNX - BW; i++) {
                        if (!g->zone_is_refined(i, j, k)) {
                            x = g->X(i, j, k);
                            dx = x - x0;
                            test = ((frac == 0 && x[0] > l1_x) || (frac == 1 && x[0] < l1_x));
                            if (test && (g->geff(i, j, k)).dot(dx) < 0.0) {
                                dphi = phi0 - (*g)(i, j, k).phi_eff();
                                if (dphi > 0.0) {
                                    const Real y0 = dphi / (8.0 * A);
                                    const Real rho0 = pow(y0 * y0 + 2.0 * y0, 3.0 / 2.0);
                                    df += -3.0 * y0 * (y0 + 1.0) * pow(rho0, 1.0 / 3.0) * dv / A;
                                    M += rho0 * dv;
                                    xmax = max(xmax, x[0]);
                                    xmin = min(xmin, x[0]);
                                }
                            }
                        }
                    }
                }
            }
        }
        tmp = df;
        MPI_Allreduce(&tmp, &df, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
        tmp = M;
        MPI_Allreduce(&tmp, &M, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
        tmp = xmax;
        MPI_Allreduce(&tmp, &xmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );
        tmp = xmin;
        MPI_Allreduce(&tmp, &xmin, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD );
        *span = xmax - xmin;
        f = M - M0;
        A -= f / df;
    } while (fabs(f / df / A) > 1.0e-6);
    return A;
}

#else

Real BinaryStar::find_K(int frac, Real phi0, Real center, Real l1_x, Real* span) {
    BinaryStar* g;
    const Real n = 1.5;
    Real K, M, tmp, dv, xmin, xmax;
    _3Vec x0, dx, x;
    bool test;
    x0 = 0.0;
    x0[0] = center;
    if (frac == 0) {
        M = M1;
    } else {
        M = M2;
    }
    xmin = 1.0e+99;
    xmax = -1.0e-99;
    Real dphi;
    K = 0.0;
    for (int l = 0; l < get_local_node_cnt(); l++) {
        g = dynamic_cast<BinaryStar*>(get_local_node(l));
        dv = pow(g->get_dx(), 3);
        for (int k = BW; k < GNX - BW; k++) {
            for (int j = BW; j < GNX - BW; j++) {
                for (int i = BW; i < GNX - BW; i++) {
                    if (!g->zone_is_refined(i, j, k)) {
                        x = g->X(i, j, k);
                        dx = x - x0;
                        test = ((frac == 0 && x[0] > l1_x) || (frac == 1 && x[0] < l1_x));
                        if (test && (g->geff(i, j, k)).dot(dx) < 0.0) {
                            dphi = phi0 - (*g)(i, j, k).phi_eff();
                            if (dphi > 0.0) {
                                K += pow(dphi / ((n + 1.0)), n) * dv;
                                xmax = max(xmax, x[0]);
                                xmin = min(xmin, x[0]);
                            }
                        }
                    }
                }
            }
        }
    }
    tmp = K;
    MPI_Allreduce(&tmp, &K, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
    tmp = xmax;
    MPI_Allreduce(&tmp, &xmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );
    tmp = xmin;
    MPI_Allreduce(&tmp, &xmin, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD );
    K = pow(K / M, 1.0 / n);
    *span = xmax - xmin;
    return K;
}
#endif

void BinaryStar::find_mass(int frac, Real* mass, _3Vec* com_ptr) {
    BinaryStar* g;
    Real M, tmp[3], com[3] = { 0.0, 0.0, 0.0 }, dm, dv, tmp0;
    M = 0.0;
    for (int l = 0; l < get_local_node_cnt(); l++) {
        g = dynamic_cast<BinaryStar*>(get_local_node(l));
        dv = pow(g->get_dx(), 3);
        for (int k = BW; k < GNX - BW; k++) {
            for (int j = BW; j < GNX - BW; j++) {
                for (int i = BW; i < GNX - BW; i++) {
                    if (!g->zone_is_refined(i, j, k)) {
                        dm = dv * (*g)(i, j, k).frac(frac);
                        M += dm;
                        com[0] += g->HydroGrid::X(i, j, k)[0] * dm;
                        com[1] += g->HydroGrid::X(i, j, k)[1] * dm;
                        com[2] += g->HydroGrid::X(i, j, k)[2] * dm;
                    }
                }
            }
        }
    }
    tmp[0] = com[0];
    tmp[1] = com[1];
    tmp[2] = com[2];
    MPI_Allreduce(tmp, com, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
    tmp0 = M;
    MPI_Allreduce(&tmp0, &M, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
    (*com_ptr)[0] = com[0];
    (*com_ptr)[1] = com[1];
    (*com_ptr)[2] = com[2];
    *com_ptr /= M;
    *mass = M;
}

int scf_iter;

Real hd(double rho) {
    const Real x = pow(rho / PhysicalConstants::B, 1.0 / 3.0);
    Real h;
    if (x < 0.01) {
        h = 4.0 * PhysicalConstants::A / PhysicalConstants::B * x * x;
    } else {
        h = 8.0 * PhysicalConstants::A / PhysicalConstants::B * (sqrt(x * x + 1.0) - 1.0);
    }
    return max(h, 0.0);
}

void BinaryStar::scf_run(int argc, char* argv[]) {
    PhysicalConstants::set_cgs();
    setup_grid_structure();

    Real& A = PhysicalConstants::A;
    Real& B = PhysicalConstants::B;
    const Real Bover8A = 0.125 * B / A;
    BinaryStar* g0;
    Real xc_d, xc_a;
    _3Vec XC_d, XC_a;
    Real rho0_d, tmp, h0_a, h0_d, h;
    Real xm_d, xp_d, phip_d, phi0_a, omega, C_a, C_d, new_rho, o2;
    const Real w = 0.9;
    Real m_a, m_d, l1_phi, l1_x;
    _3Vec com_a, com_d, com_x, last_com_x;
    Real last_w0, w0;
#ifdef USE_FMM
    const int o = 0;
#else
    const int o = BW - 1;
#endif

#ifdef USE_FMM
    FMM_solve();
#else
    solve_poisson();
#endif
    find_mass(0, &m_a, &com_a);
    find_mass(1, &m_d, &com_d);
    find_l(com_a[0], com_d[0], &l1_phi, &l1_x, 1);
    com_x = (com_a * m_a + com_d * m_d) / (m_a + m_d);

    xc_d = bparam.x2;
    xc_a = bparam.x1;
    XC_d = 0.0;
    XC_a = 0.0;
    XC_d[0] = xc_d;
    XC_a[0] = xc_a;
    xp_d = bparam.x2 + bparam.r2;
    w0 = 0.0;
    get_root()->output("S", 0, GNX, BW);
    int iter = -1;
    _3Vec O = 0.0;
    for (int iter = 0;; iter++) {
        find_l(com_a[0], com_d[0], &l1_phi, &l1_x, 1);
        phip_d = get_phi_at(xp_d, 0.0, 0.0);
        phi0_a = get_phi_at(xc_a, 0.0, 0.0);
        rho0_d = 0.0;
        for (int l = 0; l < get_local_node_cnt(); l++) {
            g0 = dynamic_cast<BinaryStar*>(get_local_node(l));
            for (int k = BW; k < GNX - BW; k++) {
                for (int j = BW; j < GNX - BW; j++) {
                    for (int i = BW; i < GNX - BW; i++) {
                        if (g0->zone_is_refined(i, j, k)) {
                            if (fabs(g0->HydroGrid::zc(k)) < g0->get_dx() / 2.0 && fabs(g0->HydroGrid::yc(j)) < g0->get_dx() / 2.0) {
                                if (fabs(g0->HydroGrid::xc(i) - xc_d) < g0->get_dx() / 2.0) {
                                    rho0_d = (*g0)(i, j, k).rho();
                                }
                            }
                        }
                    }
                }
            }
        }
        tmp = rho0_d;
        MPI_Allreduce(&tmp, &rho0_d, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
        Real phi_min_a, phi_min_d, xa, xd;
        find_phimins(&phi_min_a, &xa, &phi_min_d, &xd);
        C_d = l1_phi * (DONOR_FILL) + phi_min_d * (1.0 - DONOR_FILL);
        xm_d = l1_x;
        o2 = 2.0 * (phip_d - C_d) / (xp_d * xp_d - l1_x * l1_x);
        omega = sqrt(o2);
        State::set_omega(omega);
        h0_d = C_d - get_phi_at(xc_d, 0.0, 0.0) + 0.5 * o2 * xc_d * xc_d;
        //   printf("h0_d %e Bover8A %e get_phi_at(xc_d, 0.0, 0.0) %e xc_d %e\n", h0_d, Bover8A, get_phi_at(xc_d, 0.0, 0.0), xc_d);
        B = bparam.rho2 / pow(pow(Bover8A * h0_d + 1.0, 2) - 1.0, 1.5);
        A = B / Bover8A / 8.0;
        C_a = phi0_a - 0.5 * omega * omega * xc_a * xc_a + hd(bparam.rho1);
        XC_d[0] = xd;
        XC_a[0] = xa;
        Real verr;
        Real min_dx = HydroGrid::h0 / Real(1 << get_max_level_allowed());
        if (iter == 0) {
            verr = dynamic_cast<BinaryStar*>(get_root())->virial_error();
        }
        if (MPI_rank() == 0) {
            printf("%i %i: l1_phi=%e C_d=%e l1_phi-C_d=%e q=%e mtot=%e omega=%e  C_a=%e B=%e A=%e com_x=%e  virial_error=%e\n", get_root()->get_node_cnt(),
                    iter, l1_phi, C_d, l1_phi - C_d, m_d / m_a, m_d + m_a, omega, C_a, B, A, com_x[0], verr);
        }
        if (fabs(verr) < 2.0e-2) {
            find_mass(0, &m_a, &com_a);
            find_mass(1, &m_d, &com_d);
            com_x = (com_a * m_a + com_d * m_d) / (m_a + m_d);
            O += com_x / 2.0;
            set_origin(O);
        }
        for (int l = 0; l < get_local_node_cnt(); l++) {
            g0 = dynamic_cast<BinaryStar*>(get_local_node(l));
            for (int k = BW; k < GNX - BW; k++) {
                for (int j = BW; j < GNX - BW; j++) {
                    for (int i = BW; i < GNX - BW; i++) {
                        Real R2 = g0->HydroGrid::xc(i) * g0->HydroGrid::xc(i) + g0->HydroGrid::yc(j) * g0->HydroGrid::yc(j);
                        Real phi_eff = g0->get_phi(i - o, j - o, k - o) - 0.5 * o2 * R2;
                        if ((g0->HydroGrid::xc(i) >= xm_d)
                                && (g0->geff(i, j, k).dot(g0->X(i, j, k) - XC_d) < 0.0 || (g0->X(i, j, k) - XC_d).mag() < 2.0 * min_dx) && (phi_eff < C_d)) {
                            h = max(C_d - phi_eff, 0.0);
                            new_rho = B * pow(pow(h * Bover8A + 1.0, 2) - 1.0, 1.5);
                        } else if ((g0->HydroGrid::xc(i) <= xm_d)
                                && (g0->geff(i, j, k).dot(g0->X(i, j, k) - XC_a) < 0.0 || (g0->X(i, j, k) - XC_a).mag() < 2.0 * min_dx) && (phi_eff < C_a)) {
                            h = max(C_a - phi_eff, 0.0);
                            new_rho = B * pow(pow(h * Bover8A + 1.0, 2) - 1.0, 1.5);
                        } else {
                            new_rho = 0.0;
                        }
                        new_rho = max(new_rho, State::rho_floor);
                        (*g0)(i, j, k)[State::d_index] *= 1.0 - w;
                        (*g0)(i, j, k)[State::d_index] += w * new_rho;
                        (*g0)(i, j, k)[State::d_index] = (*g0)(i, j, k).rho();
                        if (g0->HydroGrid::xc(i) > xm_d) {
                            (*g0)(i, j, k).set_frac(1, max((*g0)(i, j, k).rho(), State::rho_floor * 0.5));
                            (*g0)(i, j, k).set_frac(0, State::rho_floor * 0.5);
                        } else {
                            (*g0)(i, j, k).set_frac(1, State::rho_floor * 0.5);
                            (*g0)(i, j, k).set_frac(0, max((*g0)(i, j, k).rho(), State::rho_floor * 0.5));
                        }
                    }
                }
            }
        }
        check_for_refine();
#ifdef USE_FMM
        FMM_solve();
#else
        solve_poisson();
#endif
        max_dt_driver();
        if (iter % 5 == 0) {
            get_root()->output("S", iter / 5 + 1, GNX, BW);
        }
        find_mass(0, &m_a, &com_a);
        find_mass(1, &m_d, &com_d);
        last_com_x = com_x;
        com_x = (com_a * m_a + com_d * m_d) / (m_a + m_d);
        verr = dynamic_cast<BinaryStar*>(get_root())->virial_error();
        //	printf("%e %e %e \n", verr, com_x, 1.0e-6 * h0 / Real(1 << get_max_level_allowed()));
        if (iter >= 15) {
            if (fabs(verr) < 1.0e-4 && fabs(com_x[0]) < 5.0e-4 * min_dx) {
                break;
            }
            if (iter > 30) {
                break;
            }
        }
    };

    find_mass(0, &m_a, &com_a);
    find_mass(1, &m_d, &com_d);
    last_com_x = com_x;
    com_x = (com_a * m_a + com_d * m_d) / (m_a + m_d);
    O += com_x;
    set_origin(O);
    find_mass(0, &m_a, &com_a);
    find_mass(1, &m_d, &com_d);
    last_com_x = com_x;
    com_x = (com_a * m_a + com_d * m_d) / (m_a + m_d);
    O += com_x;
    set_origin(O);
    Real g = pow(A, -1.5) * B * B;
    Real s = sqrt(B);
    Real cm = B / sqrt(A);
//	printf( "%e %e\n", A,B);
    PhysicalConstants::set_cgs();
//	printf( "%e %e\n", A,B);
    g *= pow(A, 1.5) / B / B;
    s *= 1.0 / sqrt(B);
    cm *= 1.0 / B * sqrt(A);
//	printf( "%e %e %e\n", g, s, cm);
    for (int l = 0; l < get_local_node_cnt(); l++) {
        g0 = dynamic_cast<BinaryStar*>(get_local_node(l));
        for (int k = BW; k < GNX - BW; k++) {
            for (int j = BW; j < GNX - BW; j++) {
                for (int i = BW; i < GNX - BW; i++) {
                    (*g0)(i, j, k)[State::d_index] *= g / (cm * cm * cm);
                }
            }
        }
    }
    omega *= 1.0 / s;
    State::set_omega(omega);
    dynamic_cast<HydroGrid*>(get_root())->HydroGrid::mult_dx(cm);
#ifndef USE_FMM
    dynamic_cast<MultiGrid*>(get_root())->MultiGrid::mult_dx(cm);
#endif
    Real dxmin = HydroGrid::h0 / Real(1 << get_max_level_allowed());
    for (int l = 0; l < get_local_node_cnt(); l++) {
        g0 = dynamic_cast<BinaryStar*>(get_local_node(l));
        for (int k = BW; k < GNX - BW; k++) {
            for (int j = BW; j < GNX - BW; j++) {
                for (int i = BW; i < GNX - BW; i++) {
                    if (g0->HydroGrid::xc(i) < xm_d) {
                        (*g0)(i, j, k)[State::frac_index + 0] = (*g0)(i, j, k).rho() - State::rho_floor / 2.0;
                        (*g0)(i, j, k)[State::frac_index + 1] = State::rho_floor / 2.0;
                    } else {
                        (*g0)(i, j, k)[State::frac_index + 1] = (*g0)(i, j, k).rho() - State::rho_floor / 2.0;
                        (*g0)(i, j, k)[State::frac_index + 0] = State::rho_floor / 2.0;
                    }
                    Real ei, ek, sy, sx, lz;
                    Real x = g0->HydroGrid::xc(i);
                    Real y = g0->HydroGrid::yc(j);
                    Real R2 = x * x + y * y;
                    sx = -(*g0)(i, j, k).rho() * omega * g0->HydroGrid::yc(j);
                    sy = +(*g0)(i, j, k).rho() * omega * g0->HydroGrid::xc(i);
                    lz = (*g0)(i, j, k).rho() * R2 * omega;
                    if ((*g0)(i, j, k).rho() < 10.0 * State::rho_floor) {
#ifdef USE_FMM
                        ei = -g0->get_phi(i, j, k) / 100.0;
#else
                        ei = -g0->get_phi(i - BW + 1, j - BW + 1, k - BW + 1) / 100.0;
#endif
                    } else {
                        ei = State::ei_floor;
                    }
                    (*g0)(i, j, k).set_lz(lz);
                    (*g0)(i, j, k).set_sx(sx);
                    (*g0)(i, j, k).set_sy(sy);
                    (*g0)(i, j, k).set_sz(0.0);
                    ek = (*g0)(i, j, k).ek(g0->X(i, j, k));
                    (*g0)(i, j, k).set_et((*g0)(i, j, k).ed() + ei + ek);
                    (*g0)(i, j, k).set_tau(pow(ei, 1.0 / State::gamma));
                }
            }
        }
    }
    State::omega0 = State::omega;
    find_mass(0, &m_a, &com_a);
    find_mass(1, &m_d, &com_d);
    com_x = (com_a * m_a + com_d * m_d) / (m_a + m_d);
#ifdef USE_FMM
    FMM_solve();
    FMM_from_children();
#else
    solve_poisson();
#endif
    max_dt_driver();
    Real verr = dynamic_cast<BinaryStar*>(get_root())->virial_error();
    if (MPI_rank() == 0) {
        FILE* fp = fopen("scf.dat", "wt");
        fprintf(fp, "\n SCF Parameters\n");
        fprintf(fp, "Virial Error   = %e\n", verr);
        fprintf(fp, "Omega          = %0.5f Hertz \n", omega);
        fprintf(fp, "Period         = %0.3f m\n", 2.0 * M_PI / omega / 60.0);
        fprintf(fp, "Accretor Mass  = %0.3f Msol\n", m_a / 1.98e+33);
        fprintf(fp, "Donor Mass     = %0.3f Msol\n", m_d / 1.98e+33);
        fprintf(fp, "Mass Ratio     = %0.3f\n", m_d / m_a);
        fprintf(fp, "Separation     = %0.3f M_R\n", (com_d[0] - com_a[0]) / 7e+10);
        fprintf(fp, "Center of Mass = %e  %e  %e cm\n", com_x[0], com_x[1], com_x[2]);
        fclose(fp);
        system("cat scf.dat\n");
    }
    get_root()->output("S", iter / 5 + 2, GNX, BW);
    if (MPI_rank() == 0) {
        system("mkdir checkpoint\n");
    }
    write_to_file(0, 0, "init");
    if (MPI_rank() == MPI_size() - 1) {
        system("mkdir output\n");
        system("mkdir scf\n");
        system("cp ./checkpoint/init.*.bin ./scf\n");
        system("mv S.*.silo ./scf\n");
        system("mv scf.dat ./scf\n");
    }
}

void BinaryStar::read_from_file(const char* str, int* i1, int* i2) {
    if (MPI_rank() == 0) {
        printf("Reading checkpoint file\n");
    }
    MPI_File fh;
    PhysicalConstants::set_cgs();
    char* fname;
    Real omega;
    FILE* fp;
    _3Vec O;
    asprintf(&fname, "./checkpoint/%s.%i.bin", str, MPI_rank());

    MPI_File_open(MPI_COMM_SELF, fname, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
    free(fname);
    MPI_File_read(fh, &refine_floor, sizeof(Real), MPI_BYTE, MPI_STATUS_IGNORE );
    MPI_File_read(fh, &State::ei_floor, sizeof(Real), MPI_BYTE, MPI_STATUS_IGNORE );
    MPI_File_read(fh, &State::rho_floor, sizeof(Real), MPI_BYTE, MPI_STATUS_IGNORE );
    MPI_File_read(fh, &com_vel_correction, sizeof(_3Vec), MPI_BYTE, MPI_STATUS_IGNORE );
    MPI_File_read(fh, &O, sizeof(_3Vec), MPI_BYTE, MPI_STATUS_IGNORE );
    MPI_File_read(fh, &State::com_correction, sizeof(_3Vec), MPI_BYTE, MPI_STATUS_IGNORE );
    MPI_File_read(fh, &State::omega, sizeof(Real), MPI_BYTE, MPI_STATUS_IGNORE );
    MPI_File_read(fh, &State::omega0, sizeof(Real), MPI_BYTE, MPI_STATUS_IGNORE );
    MPI_File_read(fh, &State::omega_dot, sizeof(Real), MPI_BYTE, MPI_STATUS_IGNORE );
    MPI_File_read(fh, &dtheta, sizeof(Real), MPI_BYTE, MPI_STATUS_IGNORE );
    MPI_File_read(fh, &last_dt, sizeof(Real), MPI_BYTE, MPI_STATUS_IGNORE );
    MPI_File_read(fh, i1, sizeof(int), MPI_BYTE, MPI_STATUS_IGNORE );
    MPI_File_read(fh, i2, sizeof(int), MPI_BYTE, MPI_STATUS_IGNORE );
    set_origin(O);
    get_root()->read_checkpoint(&fh);
    MPI_File_close(&fh);

#ifndef USE_FMM
    int count = OctNode::get_local_node_cnt();
    for (int i = 0; i < count; i++) {
        dynamic_cast<MultiGrid*>(OctNode::get_local_node(i))->phi_calc_amr_bounds();
    }
#endif
#ifdef USE_FMM
    FMM_solve();
#else
    solve_poisson();
#endif
}

void BinaryStar::write_to_file(int i1, int i2, const char* idname) {
    char* fname;
    Real omega = State::get_omega();
    FILE* fp;
    MPI_File fh;
    _3Vec O = get_origin();
    char dummy;
    asprintf(&fname, "./checkpoint/%s.%i.bin", idname, MPI_rank());
    MPI_File_open(MPI_COMM_SELF, fname, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
    free(fname);
    MPI_File_write(fh, &refine_floor, sizeof(Real), MPI_BYTE, MPI_STATUS_IGNORE );
    MPI_File_write(fh, &State::ei_floor, sizeof(Real), MPI_BYTE, MPI_STATUS_IGNORE );
    MPI_File_write(fh, &State::rho_floor, sizeof(Real), MPI_BYTE, MPI_STATUS_IGNORE );
    MPI_File_write(fh, &com_vel_correction, sizeof(_3Vec), MPI_BYTE, MPI_STATUS_IGNORE );
    MPI_File_write(fh, &O, sizeof(_3Vec), MPI_BYTE, MPI_STATUS_IGNORE );
    MPI_File_write(fh, &State::com_correction, sizeof(_3Vec), MPI_BYTE, MPI_STATUS_IGNORE );
    MPI_File_write(fh, &State::omega, sizeof(Real), MPI_BYTE, MPI_STATUS_IGNORE );
    MPI_File_write(fh, &State::omega0, sizeof(Real), MPI_BYTE, MPI_STATUS_IGNORE );
    MPI_File_write(fh, &State::omega_dot, sizeof(Real), MPI_BYTE, MPI_STATUS_IGNORE );
    MPI_File_write(fh, &dtheta, sizeof(Real), MPI_BYTE, MPI_STATUS_IGNORE );
    MPI_File_write(fh, &last_dt, sizeof(Real), MPI_BYTE, MPI_STATUS_IGNORE );
    MPI_File_write(fh, &i1, sizeof(int), MPI_BYTE, MPI_STATUS_IGNORE );
    MPI_File_write(fh, &i2, sizeof(int), MPI_BYTE, MPI_STATUS_IGNORE );
    get_root()->write_checkpoint(&fh);
    MPI_File_close(&fh);
}

void BinaryStar::analyze() {
    BinaryStar* g;
    Real tmp1, tmp2, dv;
#ifndef USE_FMM
    State U;
    _3Vec X;
    tmp1 = tmp2 = 0.0;
    for (int m = 0; m < get_local_node_cnt(); m++) {
        g = dynamic_cast<BinaryStar*>(get_local_node(m));
        dv = pow(g->get_dx(), 3);
        for (int k = BW; k < GNX - BW; k++) {
            for (int j = BW; j < GNX - BW; j++) {
                for (int i = BW; i < GNX - BW; i++) {
                    if (!g->zone_is_refined(i, j, k)) {
                        U = (*g)(i, j, k);
                        X = g->X(i, j, k);
                        tmp1 += U.lz() * dv;
                        tmp2 += U.rho() * (X[0] * X[0] + X[1] * X[1]) * dv;
                    }
                }
            }
        }
    }
    tmp1 /= tmp2;
    State::set_omega(tmp1);
//  printf("Omega = %e\n", tmp1);
    for (int m = 0; m < get_local_node_cnt(); m++) {
        g = dynamic_cast<BinaryStar*>(get_local_node(m));
        dv = pow(g->get_dx(), 3);
        for (int k = BW; k < GNX - BW; k++) {
            for (int j = BW; j < GNX - BW; j++) {
                for (int i = BW; i < GNX - BW; i++) {
                    if (!g->zone_is_refined(i, j, k)) {
                        X = g->X(i, j, k);
                        (*g)(i, j, k)[State::pot_index] = g->phi(i + 1 - BW, j + 1 - BW, k + 1 - BW) - 0.5 * (X[0] * X[0] + X[1] * X[1]) * tmp1 * tmp1;
                    }
                }
            }
        }
    }
    printf("OMega = %e\n", State::omega);
    boundary_driver();
    for (int m = 0; m < get_local_node_cnt(); m++) {
        g = dynamic_cast<BinaryStar*>(get_local_node(m));
        dv = pow(g->get_dx(), 3);
        for (int k = BW - 1; k < GNX - BW + 1; k++) {
            for (int j = BW - 1; j < GNX - BW + 1; j++) {
                for (int i = BW - 1; i < GNX - BW + 1; i++) {
                    if (!g->zone_is_refined(i, j, k)) {
                        g->phi(i + 1 - BW, j + 1 - BW, k + 1 - BW) = (*g)(i, j, k)[State::pot_index];
                    }
                }
            }
        }
    }

    Real rho_max_a = 0.0;
    _3Vec rho_max_a_loc;
    Real rho_max_d = 0.0;
    _3Vec rho_max_d_loc;
    int i0, j0, k0;
    for (int m = 0; m < get_local_node_cnt(); m++) {
        g = dynamic_cast<BinaryStar*>(get_local_node(m));
        dv = pow(g->get_dx(), 3);
        for (int k = BW; k < GNX - BW; k++) {
            if (g->HydroGrid::zc(k) < 0.0 || g->HydroGrid::zc(k) > g->get_dx()) {
                continue;
            }
            for (int j = BW; j < GNX - BW; j++) {
                if (g->HydroGrid::yc(j) < 0.0 || g->HydroGrid::yc(j) > g->get_dx()) {
                    continue;
                }
                for (int i = BW; i < GNX - BW; i++) {
                    if (!g->zone_is_refined(i, j, k)) {
                        U = (*g)(i, j, k);
                        X = g->X(i, j, k);
                        if (U.frac(0) > U.frac(1)) {
                            if (U.rho() > rho_max_a) {
                                rho_max_a_loc = X;
                                rho_max_a = U.rho();
                            }
                        } else {
                            if (U.rho() > rho_max_d) {
                                rho_max_d_loc = X;
                                rho_max_d = U.rho();
                            }
                        }
                    }
                }
            }
        }
    }
    printf("rho_acc = %e at %e %e %e\n", rho_max_a, rho_max_a_loc[0], rho_max_a_loc[1], rho_max_a_loc[2]);
    printf("rho_don = %e at %e %e %e\n", rho_max_d, rho_max_d_loc[0], rho_max_d_loc[1], rho_max_d_loc[2]);
    Real l1phi = -1.0e+99;
    _3Vec l1;
    for (int m = 0; m < get_local_node_cnt(); m++) {
        g = dynamic_cast<BinaryStar*>(get_local_node(m));
        dv = pow(g->get_dx(), 3);
        for (int k = BW; k < GNX - BW; k++) {
            if (g->HydroGrid::zc(k) < 0.0 || g->HydroGrid::zc(k) > g->get_dx()) {
                continue;
            }
            k0 = k + 1 - BW;
            for (int j = BW; j < GNX - BW; j++) {
                if (g->HydroGrid::yc(j) < 0.0 || g->HydroGrid::yc(j) > g->get_dx()) {
                    continue;
                }
                j0 = j + 1 - BW;
                for (int i = BW; i < GNX - BW; i++) {
                    if (g->HydroGrid::xc(i) < rho_max_a_loc[0] || g->HydroGrid::xc(i) > rho_max_d_loc[0]) {
                        continue;
                    }
                    i0 = i + 1 - BW;
                    if (!g->zone_is_refined(i, j, k)) {
                        U = (*g)(i, j, k);
                        X = g->X(i, j, k);
                        Real phi0 = g->phi(i0, j0, k0);
                        if (phi0 > l1phi) {
                            l1phi = phi0;
                            l1 = X;
                        }
                    }
                }
            }
        }
    }
    printf("l1 is %e at %e %e %e\n", l1phi, l1[0], l1[1], l1[2]);
    Real Mc, M1, M2;
    _3Vec a, b;
    Mc = M1 = M2 = 0.0;
    _3Vec com_a, com_d;
    com_a = 0.0;
    com_d = 0.0;
    Real J1, J2;
    for (int m = 0; m < get_local_node_cnt(); m++) {
        Real dm;
        g = dynamic_cast<BinaryStar*>(get_local_node(m));
        dv = pow(g->get_dx(), 3);
        for (int k = BW; k < GNX - BW; k++) {
            for (int j = BW; j < GNX - BW; j++) {
                for (int i = BW; i < GNX - BW; i++) {
                    if (!g->zone_is_refined(i, j, k)) {
                        dm = dv * (*g)(i, j, k).rho();
                        if (g->phi(i + 1 - BW, j + 1 - BW, k + 1 - BW) > l1phi) {
                            Mc += dm;
                        } else {
                            a = l1 - rho_max_a_loc;
                            b = g->X(i, j, k) - l1;
                            if (a.dot(b) < 0.0) {
                                M1 += dm;
                                com_a += g->X(i, j, k) * dm;
                            } else {
                                M2 += dm;
                                com_d += g->X(i, j, k) * dm;
                            }
                        }
                    }
                }
            }
        }
    }
    com_a /= M1;
    com_d /= M2;
    a = com_a - com_d;
    Real sep = sqrt(a.dot(a));
    J1 = J2 = 0.0;
    Real Jorb = 0.0;
    for (int m = 0; m < get_local_node_cnt(); m++) {
        Real dm;
        g = dynamic_cast<BinaryStar*>(get_local_node(m));
        dv = pow(g->get_dx(), 3);
        for (int k = BW; k < GNX - BW; k++) {
            for (int j = BW; j < GNX - BW; j++) {
                for (int i = BW; i < GNX - BW; i++) {
                    if (!g->zone_is_refined(i, j, k)) {
                        dm = dv * (*g)(i, j, k).rho();
                        if (g->phi(i + 1 - BW, j + 1 - BW, k + 1 - BW) > l1phi) {
                        } else {
                            a = l1 - rho_max_a_loc;
                            b = g->X(i, j, k) - l1;
                            State U = g->U(i, j, k);
                            _3Vec x = g->X(i, j, k);
                            Jorb += (x[0] * U.sy() - x[1] * U.sx()) * dv;
                            if (a.dot(b) < 0.0) {
                                x -= rho_max_a_loc;
                                J1 += (x[0] * U.sy() - x[1] * U.sx()) * dv;
                            } else {
                                x -= rho_max_d_loc;
                                J2 += (x[0] * U.sy() - x[1] * U.sx()) * dv;
                            }
                        }
                    }
                }
            }
        }
    }
    Jorb -= J1 + J2;
    FILE* fp = fopen("silo.dat", "at");
    fprintf(fp, "%e %e %e %e %e %e %e %e %e %e %e %e\n", dynamic_cast<const BinaryStar*>(get_root())->get_time(), State::get_omega(), l1[3], l1[0], l1[1], M1,
            M2, Mc, sep, J1, J2, Jorb);
    fclose(fp);

    /*Real dphi_dx_p, dphi_dx_m, dphi_dy_p, dphi_dy_m, d2phi_dx2, d2phi_dy2;
     Vector<Real, 4> minima[256];
     Vector<Real, 4> saddles[256];
     int n_min, n_sad;
     n_min = n_sad = 0;
     for (int m = 0; m < get_local_node_cnt(); m++) {
     g = dynamic_cast<BinaryStar*>(get_local_node(m));
     dv = pow(g->get_dx(), 3);
     for (int k = BW; k < GNX - BW; k++) {
     if (g->HydroGrid::zc(k) < 0.0 || g->HydroGrid::zc(k) > g->get_dx()) {
     continue;
     }
     k0 = k + 1 - BW;
     for (int j = BW; j < GNX - BW; j++) {
     if (g->HydroGrid::yc(j) < 0.0 || g->HydroGrid::yc(j) > g->get_dx()) {
     continue;
     }
     j0 = j + 1 - BW;
     for (int i = BW; i < GNX - BW; i++) {
     i0 = i + 1 - BW;
     if (!g->zone_is_refined(i, j, k)) {
     U = (*g)(i, j, k);
     X = g->X(i, j, k);
     dphi_dx_p = g->phi(i0 + 1, j0, k0) - g->phi(i0, j0, k0);
     dphi_dx_m = g->phi(i0, j0, k0) - g->phi(i0 - 1, j0, k0);
     dphi_dy_p = g->phi(i0, j0 + 1, k0) - g->phi(i0, j0, k0);
     dphi_dy_m = g->phi(i0, j0, k0) - g->phi(i0, j0 - 1, k0);
     if (dphi_dx_p * dphi_dx_m < 0.0) {
     Vector<Real, 4> v;
     d2phi_dx2 = dphi_dx_p - dphi_dx_m;
     d2phi_dy2 = dphi_dy_p - dphi_dy_m;
     v[0] = X[0] - 0.5 * g->get_dx() * (dphi_dx_p + dphi_dx_m) / d2phi_dx2;
     v[1] = X[1] - 0.5 * g->get_dx() * (dphi_dy_p + dphi_dy_m) / d2phi_dy2;
     v[2] = X[2];
     v[3] = g->phi(i0, j0, k0);
     if (d2phi_dx2 < 0.0) {
     saddles[n_sad] = v;
     printf("sad %e %e %e %e\n", saddles[n_sad][0], saddles[n_sad][1], saddles[n_sad][2], saddles[n_sad][3]);
     n_sad++;
     } else if (d2phi_dy2 > 0.0) {
     minima[n_min] = v;
     printf("min %e %e %e %e\n", minima[n_min][0], minima[n_min][1], minima[n_min][2], minima[n_min][3]);
     n_min++;
     }
     }
     }
     }
     }
     }
     }

     if (n_sad == 3 && n_min == 2) {
     Vector<Real, 4> l1;
     _3Vec a, b;
     for (int i = 0; i < n_sad; i++) {
     for (int j = 0; j < 3; j++) {
     a[j] = saddles[i][j] - minima[0][j];
     b[j] = saddles[i][j] - minima[1][j];
     }
     if (a.dot(b) < 0.0) {
     l1 = saddles[i];
     }
     printf("%e %e %e %e\n", saddles[i][0], saddles[i][1], saddles[i][2], saddles[i][3]);
     }
     Real M1, M2, Mc, dm;
     M1 = M2 = Mc = 0.0;
     _3Vec M1x;
     int index;
     if (minima[0][3] < minima[1][3]) {
     index = 0;
     } else {
     index = 1;
     }
     for (int j = 0; j < 3; j++) {
     M1x[j] = minima[index][j];
     }
     for (int m = 0; m < get_local_node_cnt(); m++) {
     g = dynamic_cast<BinaryStar*>(get_local_node(m));
     dv = pow(g->get_dx(), 3);
     for (int k = BW; k < GNX - BW; k++) {
     for (int j = BW; j < GNX - BW; j++) {
     for (int i = BW; i < GNX - BW; i++) {
     if (!g->zone_is_refined(i, j, k)) {
     dm = dv * (*g)(i, j, k).rho();
     if (g->phi(i + 1 - BW, j + 1 - BW, k + 1 - BW) > l1[3]) {
     Mc += dm;
     } else {
     for (int l = 0; l < 3; l++) {
     a[l] = l1[l] - M1x[l];
     b[l] = g->X(i, j, k)[l] - l1[l];
     }
     if (a.dot(b) < 0.0) {
     M1 += dm;
     } else {
     M2 += dm;
     }
     }
     }
     }
     }
     }
     }

     FILE* fp = fopen("silo.dat", "at");
     fprintf(fp, "%e %e %e %e %e %e %e %e\n", dynamic_cast<const BinaryStar*>(get_root())->get_time(), State::get_omega(), l1[3], l1[0], l1[1], M1, M2, Mc);
     fclose(fp);
     } else {
     printf("Bad Roch %i %i\n", n_sad, n_min);
     }*/
#endif
}

#define RK2

void BinaryStar::step(Real dt) {
    Real start_time;
#ifdef RK2
    Real beta[2] = { 1.0, 1.0 / 2.0 };
    const int nstep = 2;
#else
    Real beta[3] = {1.0, 0.25, 2.0 / 3.0};
    const int nstep = 2;
#endif
    HydroGrid::set_dt(dt);
    store();
    store_pot();
    Real dtheta0 = dtheta;
    Real omega0 = State::omega;
    for (int i = 0; i < nstep; i++) {
        HydroGrid::set_beta(beta[i]);
        start_time = MPI_Wtime();
        substep_driver();
#ifdef USE_FMM
        FMM_solve_dot();
        update();
        dtheta = (dtheta + dt * (State::omega - State::omega0)) * beta[i] + dtheta0 * (1.0 - beta[i]);
        //  State::omega = (State::omega + dt * State::omega_dot) * beta[i] + omega0 * (1.0 - beta[i]);
        FMM_solve();
        account_pot();
#else
        solve_poisson();
#endif
    }
    Real theta, dethea;
    adjust_frame(dt);
    compute_axis(dt, &theta, &dtheta);
#ifdef USE_FMM
    FMM_from_children();
#endif
    pot_to_hydro_grid();
    set_time(get_time() + dt);

}

void BinaryStar::run(int argc, char* argv[]) {
    int mingrids, maxgrids;
    char* silo_in_file;
    mingrids = 0;
    maxgrids = 1000000000;
    if (argc == 4) {
        bparam.m1 = atof(argv[2]);
        bparam.m2 = atof(argv[3]);
    } else if (argc == 5) {
        mingrids = atoi(argv[3]);
        maxgrids = atoi(argv[4]);
    } else if (argc == 2) {
        read_silo(argv[1]);
        HydroGrid::redistribute_grids();
        analyze();
        return;
    } else if (argc != 3) {
        printf("Missing arguments\n");
        return;
    }
//	State::turn_off_total_energy();

    /*	if (scf_code) {
     scf_run(argc, argv);
     return;
     }*/
#ifndef USE_FMM
    set_poisson_tolerance(1.0e-8);
#endif
    shadow_off();

    Real dt;
    bool do_output, last_step;
    int step_cnt = 0;
    int ostep_cnt = 0;
    int nnodes;

    if (argc == 4) {
//		scf_code = true;
        scf_run(argc, argv);
        return;
        scf_code = false;
        PhysicalConstants::set_cgs();
//						setup_grid_structure();
    } else {
        read_from_file(argv[2], &step_cnt, &ostep_cnt);
    }

    PhysicalConstants::set_cgs();

#ifndef USE_FMM
    set_poisson_tolerance(1.0e-10);
#endif
    if (MPI_rank() == 0) {
        printf("Beginning evolution with period = %e, omega = %e, refine_floor = %e\n", 2.0 * M_PI / State::get_omega(), State::get_omega(), refine_floor);
    }
#ifndef USE_FMM
    hydro_time = poisson_boundary_time = poisson_interior_time = 0.0;
#endif
    State sum;
    if (get_time() == 0.0) {
#ifdef USE_FMM
        FMM_solve();
#else
        solve_poisson();
#endif
        get_root()->output("./output/X", 0.0, GNX, BW, dtheta);
        integrate_conserved_variables(&sum);
        lz_t0 = sum[State::lz_index];
    }
    int iter = 0;
    Real time_start = MPI_Wtime();
    do {
        iter++;
        if (iter == 10) {
            //		break;
        }
        integrate_conserved_variables(&sum);
        if (MPI_rank() == 0) {
            FILE* fp = fopen("sums.dat", "at");
            fprintf(fp, "%e ", get_time());
            for (int l = 0; l < STATE_NF; l++) {
                fprintf(fp, "%.10e ", sum[l]);
            }
            fprintf(fp, "\n");
            fclose(fp);
        }
        Real ofreq = (OUTPUT_TIME_FREQ * (2.0 * M_PI) / State::omega0);
//		Real ofreq = (OUTPUT_TIME_FREQ * (2.0 * M_PI) / State::get_omega());

        if (step_cnt % CHECKPT_FREQ == 0) {
            if (step_cnt % (CHECKPT_FREQ * 2) == 0) {
                write_to_file(step_cnt, ostep_cnt, "hello");
            } else {
                write_to_file(step_cnt, ostep_cnt, "goodbye");
            }
        }
        dt = next_dt(&do_output, &last_step, &ostep_cnt, ofreq);
        step(dt);
        //compute_omega_dot(dt);
        pot_to_hydro_grid();
        Real theta, theta_dot;
        //     compute_axis(dt, &theta, &theta_dot);
#ifndef USE_FMM
        _3Vec com = get_center_of_mass();
#else
        _3Vec com = system_com();
#endif
        Real com_omega = State::get_omega() * 200.0;
        com_vel_correction -= (com_vel_correction * 2.0 + com * com_omega) * dt * com_omega;
        if (MPI_rank() == 0) {
            FILE* fp = fopen("com.dat", "at");
            fprintf(fp, "%.12e %.12e %.12e %.12e\n", get_time(), com[0], com[1], com[2]);
            fclose(fp);
        }

        diagnostics(dt);

        if (get_root()->get_max_level() < get_root()->get_max_level_allowed()) {
            printf("FATAL ERROR DETECTED - levels\n");
            abort();
        }
        if (dt < 0.0) {
            printf("Negative DT\n");
            abort();
        }

        if (step_cnt % (GNX - 2 * BW)== 0){
            check_for_refine();
            if (get_root()->get_node_cnt() > maxgrids) {
                BinaryStar::refine_floor *= pow(10.0, 0.1);
                check_for_refine();
            } else if (get_root()->get_node_cnt() < mingrids) {
                BinaryStar::refine_floor /= pow(10.0, 0.1);
                check_for_refine();
            }
        }
        step_cnt++;

        if (MPI_rank() == 0) {
            nnodes = OctNode::get_node_cnt();
            FILE* fp = fopen("dt.dat", "at");
            fprintf(fp, "%e %e %i\n", HydroGrid::get_time(), dt, nnodes);
            fclose(fp);
            printf("step=%i t=%e dt=%e lmax=%i ngrids=%i step_time=%e s refine_floor=%e", step_cnt, HydroGrid::get_time(), dt, OctNode::get_max_level(), nnodes,
                    MPI_Wtime() - time_start, BinaryStar::refine_floor);
        }
        time_start = MPI_Wtime();
        if (do_output) {
            if (MPI_rank() == 0) {
                printf("*");
            }
            get_root()->output("./output/X", nint(HydroGrid::get_time() / ofreq), GNX, BW, dtheta);
        }
        if (MPI_rank() == 0) {
            printf("\n");
        }
    } while (!last_step);
    /*if (MPI_rank() == 0) {
     char* str;
     if (asprintf(&str, "time.%i.txt", get_max_level_allowed()) == 0) {
     printf("Unable to create filename\n");
     } else {
     FILE* fp = fopen(str, "at");
     fprintf(fp, "%i %e %e %e %e\n", MPI_size(), hydro_time + poisson_boundary_time + poisson_interior_time, hydro_time, poisson_boundary_time,
     poisson_interior_time);
     fclose(fp);
     free(str);
     }
     }*/
}

void BinaryStar::integrate_conserved_variables(Vector<Real, STATE_NF>* sum_ptr) {
    int i;
    BinaryStar* g;
    Vector<Real, STATE_NF> my_sum = 0.0;
    Real dv;
    for (int m = 0; m < get_local_node_cnt(); m++) {
        g = dynamic_cast<BinaryStar*>(get_local_node(m));
        dv = pow(g->get_dx(), 3);
        for (int k = BW; k < GNX - BW; k++) {
            for (int j = BW; j < GNX - BW; j++) {
                for (int i = BW; i < GNX - BW; i++) {
                    if (!g->zone_is_refined(i, j, k)) {
                        for (int l = 0; l < STATE_NF; l++) {
                            if (l == State::et_index) {
                                my_sum[l] += (*g)(i, j, k).conserved_energy(g->HydroGrid::X(i, j, k)) * dv;
                            } else {
                                my_sum[l] += (*g)(i, j, k)[l] * dv;
                            }
                        }
                    }
                }
            }
        }
    }
    MPI_Allreduce(&my_sum, sum_ptr, STATE_NF, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
}

void BinaryStar::diagnostics(Real dt) {
    int i;
    BinaryStar* b;
    Real dv;
    State u, dudt;
    _3Vec x;
    Real tmp1[6], tmp2[6];
    Real sz = 0.0, sx = 0.0, sy = 0.0;
    Real lz = 0.0;
    Real da;
    Real etall = 0.0;
    Real mtot = 0.0;
    for (int m = 0; m < get_local_node_cnt(); m++) {
        b = dynamic_cast<BinaryStar*>(get_local_node(m));
        dv = pow(b->get_dx(), 3);
        da = pow(b->get_dx(), 2);
        for (int k = BW; k < GNX - BW; k++) {
            for (int j = BW; j < GNX - BW; j++) {
                for (int i = BW; i < GNX - BW; i++) {
                    if (!b->zone_is_refined(i, j, k)) {
                        u = (*b)(i, j, k);
                        dudt = b->get_dudt(i, j, k);
                        x = b->X(i, j, k);
                        etall += u.conserved_energy(x) * dv;
                        mtot += u[State::d_index] * dv;
                        lz += u.lz() * dv;
                        sx += u.sx() * dv;
                        sy += u.sy() * dv;
                        sz += u.sz() * dv;
                        const int l = State::lz_index;
                    }
                }
            }
        }
    }
    tmp1[0] = sx;
    tmp1[1] = sy;
    tmp1[2] = sz;
    tmp1[3] = mtot;
    tmp1[4] = lz;
    tmp1[5] = etall;
    MPI_Reduce(tmp1, tmp2, 6, MPI_DOUBLE, MPI_SUM, get_root()->proc(), MPI_COMM_WORLD );
    if (MPI_rank() == get_root()->proc()) {
        sx = tmp2[0];
        sy = tmp2[1];
        sz = tmp2[2];
        mtot = tmp2[3];
        lz = tmp2[4];
        etall = tmp2[5];
        etall += State::omega * lz;
        etall += dynamic_cast<const BinaryStar*>(get_root())->get_flow_off()[State::et_index];
        //  etall += dynamic_cast<const BinaryStar*>(get_root())->get_flow_off()[State::pot_index];
        // etall += State::omega * dynamic_cast<const BinaryStar*>(get_root())->get_flow_off()[State::lz_index];
        lz += dynamic_cast<const BinaryStar*>(get_root())->get_flow_off()[State::lz_index];
        sx += dynamic_cast<const BinaryStar*>(get_root())->get_flow_off()[State::sx_index];
        sy += dynamic_cast<const BinaryStar*>(get_root())->get_flow_off()[State::sy_index];
        sz += dynamic_cast<const BinaryStar*>(get_root())->get_flow_off()[State::sz_index];
        mtot += dynamic_cast<const BinaryStar*>(get_root())->get_flow_off()[State::d_index];
        FILE* fp = fopen("diag.dat", "at");
        fprintf(fp, "%.9e %.15e %.15e %.15e %.15e %.15e %.15e \n", get_time(), lz, etall, mtot, sx, sy, sz);
        fclose(fp);
    }
}

void BinaryStar::find_phimins(Real* phi_min_a, Real* a_x, Real* phi_min_d, Real* d_x) {
    BinaryStar* g;
    Vector<int, 3> loc;
    Real pma = +1.0e+99;
    Real pmd = +1.0e+99;
    int i, j, k;
    int plane_index;
    Real xa, xd;
    Real r[2], r0[2];
    for (int l = 0; l < get_local_node_cnt(); l++) {
        g = dynamic_cast<BinaryStar*>(get_local_node(l));
        loc = g->get_location();
        plane_index = (1 << g->get_level()) / 2;
        if (loc[1] == plane_index && loc[2] == plane_index) {
            j = k = BW;
            for (i = BW; i < GNX - BW; i++) {
                if (!g->zone_is_refined(i, j, k) && (*g)(i, j, k).rho() > 10.0 * State::rho_floor) {
                    if ((*g)(i, j, k).frac(0) < (*g)(i, j, k).frac(1)) {
                        if ((*g)(i, j, k).phi_eff() < pmd) {
                            pmd = (*g)(i, j, k).phi_eff();
                            xd = g->HydroGrid::xc(i);
                        }
                    } else {
                        if ((*g)(i, j, k).phi_eff() < pma) {
                            pma = (*g)(i, j, k).phi_eff();
                            xa = g->HydroGrid::xc(i);
                        }
                    }
                }
            }
        }
    }
    r[0] = pma;
    r[1] = xa;
    MPI_Allreduce(r, r0, 1, MPI_2DOUBLE_PRECISION, MPI_MINLOC, MPI_COMM_WORLD );
    *phi_min_a = r0[0];
    *a_x = r0[1];
    r[0] = pmd;
    r[1] = xd;
    MPI_Allreduce(r, r0, 1, MPI_2DOUBLE_PRECISION, MPI_MINLOC, MPI_COMM_WORLD );
    *phi_min_d = r0[0];
    *d_x = r0[1];
//	printf("%e %e \n", xd, pmd);
//	printf("phi_min_a = %e\n", phi_min_a);
}

Real BinaryStar::find_acc_com() {
    BinaryStar* g;
    Real com = 0.0;
    Real tmp, dv, mtot = 0.0, dm;
    int i, j, k;
    for (int l = 0; l < get_local_node_cnt(); l++) {
        g = dynamic_cast<BinaryStar*>(get_local_node(l));
        dv = pow(g->get_dx(), 3);
        for (k = BW; k < GNX - BW; k++) {
            for (j = BW; j < GNX - BW; j++) {
                for (i = BW; i < GNX - BW; i++) {
                    dm = (*g)(i, j, k).frac(0) * dv;
                    com += g->HydroGrid::xc(i) * dm;
                    mtot += dm;
                }
            }
        }
    }
    tmp = com;
    MPI_Allreduce(&tmp, &com, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
    tmp = mtot;
    MPI_Allreduce(&tmp, &mtot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
    return com / mtot;
}

void BinaryStar::next_omega(Real Kd) {
    BinaryStar* g;
    int i, j, k;
    Real frot, ftot, x, hp, hm, phip, phim, tmp;
    Real n = 1.5, d, dp, dm;
#ifdef USE_FMM
    const int o = 0;
#else
    const int o = BW - 1;
#endif
    frot = ftot = 0.0;
    for (int l = 0; l < get_local_node_cnt(); l++) {
        g = dynamic_cast<BinaryStar*>(get_local_node(l));
        Real dv = g->get_dx() * g->get_dx() * g->get_dx();
        for (int k = BW; k < GNX - BW; k++) {
            for (int j = BW; j < GNX - BW; j++) {
                for (int i = BW; i < GNX - BW; i++) {
                    if (!g->zone_is_refined(i, j, k) && g->HydroGrid::xc(i) < 0.0) {
                        d = (*g)(i, j, k).rho();
                        dp = (*g)(i + 1, j, k).rho();
                        dm = (*g)(i - 1, j, k).rho();
                        hp = Kd * (n + 1.0) * pow(max(dp, 0.0), 1.0 / n);
                        hm = Kd * (n + 1.0) * pow(max(dm, 0.0), 1.0 / n);
                        phip = g->get_phi(i - o + 1, j - o, k - o);
                        phim = g->get_phi(i - o - 1, j - o, k - o);
                        tmp = d * g->HydroGrid::xc(i) * dv;
                        ftot += d * ((hp - hm) + (phip - phim)) / (2.0 * g->get_dx()) * dv;
                        frot += tmp;
                    }
                }
            }
        }
    }
    tmp = ftot;
    MPI_Allreduce(&tmp, &ftot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
    tmp = frot;
    MPI_Allreduce(&tmp, &frot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
//printf("%e %e\n", frot, ftot);
    Real oldO = State::get_omega();
    Real newO = sqrt(ftot / frot);
    State::set_omega(oldO * pow(newO / oldO, 0.05));
}

void BinaryStar::find_l(Real m1_x, Real m2_x, Real* l1_phi, Real* l1_x, int lnum) {
#ifdef USE_FMM
    const int o = 0;
#else
    const int o = BW - 1;
#endif
    BinaryStar* g;
    Vector<int, 3> loc;
    Real phi_max = -1.0e+99;
    int i, j, k;
    int plane_index;
    Real x, x_max;
    Real r[2], r0[2];
    Real O2 = pow(State::get_omega(), 2);
    bool test;
    for (int l = 0; l < get_local_node_cnt(); l++) {
        g = dynamic_cast<BinaryStar*>(get_local_node(l));
        loc = g->get_location();
        plane_index = (1 << g->get_level()) / 2;
        if (loc[1] == plane_index && loc[2] == plane_index) {
            j = k = BW;
            for (i = BW; i < GNX - BW; i++) {
                if (!g->zone_is_refined(i, j, k)) {
                    x = g->HydroGrid::xc(i);
                    if (lnum == 1) {
                        test = x < m2_x && x > m1_x;
                    } else {
                        test = x > m2_x;
                    }
                    if (test) {
                        Real phi_eff = g->get_phi(i - o, j - o, k - o) - 0.5 * (x * x + 0.25 * g->get_dx() * g->get_dx()) * O2;
                        if (phi_eff > phi_max) {
                            phi_max = phi_eff;
                            x_max = x;
                        }
                    }
                }
            }
        }
    }

    r[0] = phi_max;
    r[1] = x_max;
    MPI_Allreduce(r, r0, 1, MPI_2DOUBLE_PRECISION, MPI_MAXLOC, MPI_COMM_WORLD );
    *l1_phi = r0[0];
    *l1_x = r0[1];
}
#endif
#endif
