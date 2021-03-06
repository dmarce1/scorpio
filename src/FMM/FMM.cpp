/*
 * hydro_FMM_grid.cpp
 *
 *  Created on: Sep 20, 2013
 *      Author: dmarce1
 */

#include "FMM.h"

#ifdef USE_FMM

#include "../tag.h"
#include <mpi.h>
#include <omp.h>

#define XLYL 0
#define XLYU 1
#define XLZL 2
#define XLZU 3
#define XUYL 4
#define XUYU 5
#define XUZL 6
#define XUZU 7
#define YLZL 8
#define YLZU 9
#define YUZL 10
#define YUZU 11
#define XLYLZL 0
#define XLYLZU 1
#define XLYUZL 2
#define XLYUZU 3
#define XUYLZL 4
#define XUYLZU 5
#define XUYUZL 6
#define XUYUZU 7

static int face_id[26];
static int face_opp_id[26];

FMM::ifunc_t FMM::cs[FSTAGE + 1] = { &FMM::moments_recv, &FMM::moments_recv_wait, &FMM::moments_send, &FMM::moments_send_wait, &FMM::moments_communicate_all,
        &FMM::moments_communicate_wait_all, &FMM::compute_interactions, &FMM::expansion_recv, &FMM::expansion_recv_wait, &FMM::expansion_send,
        &FMM::expansion_send_wait, &FMM::null };

FMM::ifunc_t FMM::cs_dot[FSTAGE + 1] = { &FMM::moments_recv_dot, &FMM::moments_recv_wait, &FMM::moments_send_dot, &FMM::moments_send_wait_dot,
        &FMM::moments_communicate_all_dot, &FMM::moments_communicate_wait_all, &FMM::compute_interactions_dot, &FMM::expansion_recv_dot,
        &FMM::expansion_recv_wait_dot, &FMM::expansion_send_dot, &FMM::expansion_send_wait, &FMM::null };
FMM::ifunc_t FMM::cs_children[5] = { &FMM::_4force_recv, &FMM::_4force_recv_wait, &FMM::_4force_send, &FMM::_4force_send_wait, &FMM::null };

FMM::ifunc_t FMM::fmmhydro[40] = { &Hydro::flux_bnd_comm, &Hydro::physical_boundary, &Hydro::flux_bnd_recv_wait, &Hydro::amr_bnd_send,
        &Hydro::recv_from_children_dt, &Hydro::flux_compute, &Hydro::recv_from_children_dt_wait, &Hydro::amr_bnd_send_wait, &Hydro::flux_cf_adjust_recv,
        &Hydro::flux_cf_adjust_recv_wait, &Hydro::flux_cf_adjust_send, &Hydro::compute_dudt, &FMM::moments_recv_dot, &FMM::moments_recv_wait,
        &FMM::moments_send_dot, &FMM::moments_communicate_all_dot, &FMM::moments_communicate_wait_all, &FMM::expansion_recv_dot, &FMM::compute_interactions_dot,
        &FMM::expansion_recv_wait_dot, &FMM::expansion_send_dot, &FMM::compute_gravity_update, &FMM::moments_recv, &FMM::moments_recv_wait, &FMM::moments_send,
        &FMM::moments_communicate_all, &FMM::moments_communicate_wait_all, &FMM::expansion_recv, &FMM::compute_interactions, &FMM::expansion_recv_wait,
        &FMM::expansion_send, &FMM::copy_pot_to_hydro_grid, &FMM::apply_pot_change, &Hydro::inject_from_children_recv, &Hydro::inject_from_children_recv_wait,
        &Hydro::flux_cf_adjust_send_wait, &Hydro::inject_from_children_send, &FMM::moments_send_wait_dot, &FMM::moments_send_wait,
        &Hydro::inject_from_children_send_wait, };

MPI_Datatype FMM::MPI_multipole_t;
MPI_Datatype FMM::MPI_send_bnd_t[26];
MPI_Datatype FMM::MPI_recv_bnd_t[26];
MPI_Datatype FMM::MPI_comm_child_poles_t[8];
MPI_Datatype FMM::MPI_comm_taylor_t[8];
MPI_Datatype FMM::MPI_taylor_t;
MPI_Datatype FMM::MPI_multipole_dot_t;
MPI_Datatype FMM::MPI_send_bnd_dot_t[26];
MPI_Datatype FMM::MPI_recv_bnd_dot_t[26];
MPI_Datatype FMM::MPI_comm_child_poles_dot_t[8];
MPI_Datatype FMM::MPI_comm_taylor_dot_t[8];
MPI_Datatype FMM::MPI_taylor_dot_t;
MPI_Datatype FMM::MPI_comm_child3_t[8];
MPI_Datatype FMM::MPI_4force_t;
Real FMM::d0_array[INX][INX][INX];
Real FMM::d1_array[2 * INX + 1][2 * INX + 1][2 * INX + 1][3];
bool FMM::solve_on = true;

void FMM::moments_recv(int) {

    int tag;
    FMM* child;
    const Real dv = pow(get_dx(), 3);
    for (int k = 0; k < FNX; k++) {
        for (int j = 0; j < FNX; j++) {
            for (int i = 0; i < FNX; i++) {
                poles(i, j, k).M = 0.0;
                poles(i, j, k).X = 0.0;
                poles(i, j, k).is_leaf = false;
            }
        }
    }
    for (ChildIndex ci = 0; ci < 8; ci++) {
        child = dynamic_cast<FMM*>(get_child(ci));
        if (child != NULL) {
            tag = tag_gen(TAG_FMM_MOMENT, get_id(), ci);
            MPI_Irecv(poles.ptr(), 1, MPI_comm_child_poles_t[ci], get_child(ci)->proc(), tag, MPI_COMM_WORLD, recv_request + ci);
        } else {
            recv_request[ci] = MPI_REQUEST_NULL;
            const int xlb = FBW + ci.get_x() * (INX / 2);
            const int ylb = FBW + ci.get_y() * (INX / 2);
            const int zlb = FBW + ci.get_z() * (INX / 2);
            const int xub = xlb + (INX / 2) - 1;
            const int yub = ylb + (INX / 2) - 1;
            const int zub = zlb + (INX / 2) - 1;
//#pragma omp parallel for collapse(2)
            for (int k = zlb; k <= zub; k++) {
                for (int j = ylb; j <= yub; j++) {
                    for (int i = xlb; i <= xub; i++) {
                        poles(i, j, k).X[0] = xc(i + BW - FBW);
                        poles(i, j, k).X[1] = yc(j + BW - FBW);
                        poles(i, j, k).X[2] = zc(k + BW - FBW);
                        poles(i, j, k).M() = ((*this)(i + BW - FBW, j + BW - FBW, k + BW - FBW).rho() - State::rho_floor) * dv;
                        poles(i, j, k).is_leaf = true;
                    }
                }
            }
        }
    }
    inc_instruction_pointer();
}

void FMM::moments_recv_dot(int) {

    int tag;
    FMM* child;
    const Real dv = pow(get_dx(), 3);
    for (int k = 0; k < FNX; k++) {
        for (int j = 0; j < FNX; j++) {
            for (int i = 0; i < FNX; i++) {
                poles_dot(i, j, k).M_dot = 0.0;
            }
        }
    }
    for (ChildIndex ci = 0; ci < 8; ci++) {
        child = dynamic_cast<FMM*>(get_child(ci));
        if (child != NULL) {
            tag = tag_gen(TAG_FMM_MOMENT, get_id(), ci);
            MPI_Irecv(poles_dot.ptr(), 1, MPI_comm_child_poles_dot_t[ci], get_child(ci)->proc(), tag, MPI_COMM_WORLD, recv_request + ci);
        } else {
            recv_request[ci] = MPI_REQUEST_NULL;
            const int xlb = FBW + ci.get_x() * (INX / 2);
            const int ylb = FBW + ci.get_y() * (INX / 2);
            const int zlb = FBW + ci.get_z() * (INX / 2);
            const int xub = xlb + (INX / 2) - 1;
            const int yub = ylb + (INX / 2) - 1;
            const int zub = zlb + (INX / 2) - 1;
//#pragma omp parallel for collapse(2)
            for (int k = zlb; k <= zub; k++) {
                for (int j = ylb; j <= yub; j++) {
                    for (int i = xlb; i <= xub; i++) {
                        Real rnf = 0.0;
                        for (int fi = 0; fi < NFRAC; fi++) {
                            rnf += get_dudt(i + BW - FBW, j + BW - FBW, k + BW - FBW)[State::frac_index + fi];
                        }
                        poles_dot(i, j, k).M_dot() = rnf * dv;
                    }
                }
            }
        }
    }
    inc_instruction_pointer();
}

void FMM::moments_recv_wait(int) {
    int flag;
    MPI_Testall(8, recv_request, &flag, MPI_STATUS_IGNORE);
    if (flag) {
        inc_instruction_pointer();
    }
}

void FMM::moments_send(int) {
    if (get_level() != 0) {
        int cnt;
        _3Vec Y;
        multipole_t tmp;
        Real mass;
        cnt = 0;
        moment_buffer = new multipole_t[INX * INX * INX / 8];
        for (int k = FBW; k < FNX - FBW; k += 2) {
            for (int j = FBW; j < FNX - FBW; j += 2) {
                for (int i = FBW; i < FNX - FBW; i += 2) {
                    tmp.M = 0.0;
                    tmp.X = 0.0;
                    tmp.is_leaf = false;
                    Real mass = 0.0;
//#pragma omp parallel for collapse(2) reduction(+:mass)
                    for (int a = i; a < i + 2; a++) {
                        for (int b = j; b < j + 2; b++) {
                            for (int c = k; c < k + 2; c++) {
                                multipole_t& child = poles(a, b, c);
                                mass += child.M();
//#pragma omp critical
                                tmp.X += child.X * child.M();
                            }
                        }
                    }
                    if (mass != 0.0) {
                        tmp.X /= mass;
                    } else {
                        tmp.X = 0.0;
                        for (int a = i; a < i + 2; a++) {
                            for (int b = j; b < j + 2; b++) {
                                for (int c = k; c < k + 2; c++) {
                                    tmp.X += poles(a, b, c).X;
                                }
                            }
                        }
                        tmp.X *= 0.125;
                    }
//#pragma omp parallel for collapse(2) private(Y)
                    for (int a = i; a < i + 2; a++) {
                        for (int b = j; b < j + 2; b++) {
                            for (int c = k; c < k + 2; c++) {
                                multipole_t& child = poles(a, b, c);
                                Y = child.X - tmp.X;
//#pragma omp critical
                                tmp.M += child.M >> Y;
                            }
                        }
                    }
                    Xp(i / 2, j / 2, k / 2) = tmp.X;
                    moment_buffer[cnt] = tmp;
                    cnt++;
                }
            }
        }
        assert(cnt==INX * INX * INX / 8);
        int tag = tag_gen(TAG_FMM_MOMENT, get_parent()->get_id(), my_child_index());
        MPI_Isend(moment_buffer, cnt * sizeof(multipole_t), MPI_BYTE, get_parent()->proc(), tag, MPI_COMM_WORLD, send_request);
    }
    inc_instruction_pointer();
}

void FMM::moments_send_dot(int) {
    if (get_level() != 0) {
        int cnt;
        _3Vec Y;
        multipole_dot_t tmp;
        Real mass;
        cnt = 0;
        moment_dot_buffer = new multipole_dot_t[INX * INX * INX / 8];
        for (int k = FBW; k < FNX - FBW; k += 2) {
            for (int j = FBW; j < FNX - FBW; j += 2) {
                for (int i = FBW; i < FNX - FBW; i += 2) {
                    tmp.M_dot = 0.0;
                    Real mass = 0.0;
//#pragma omp parallel for collapse(2) private(Y)
                    for (int a = i; a < i + 2; a++) {
                        for (int b = j; b < j + 2; b++) {
                            for (int c = k; c < k + 2; c++) {
                                multipole_t& child = poles(a, b, c);
                                Y = child.X - Xp(i / 2, j / 2, k / 2);
//#pragma omp critical
                                tmp.M_dot += poles_dot(a, b, c).M_dot >> Y;
                            }
                        }
                    }
                    moment_dot_buffer[cnt] = tmp;
                    cnt++;
                }
            }
        }
        assert(cnt==INX * INX * INX / 8);
        int tag = tag_gen(TAG_FMM_MOMENT, get_parent()->get_id(), my_child_index());
        MPI_Isend(moment_dot_buffer, cnt * sizeof(multipole_dot_t), MPI_BYTE, get_parent()->proc(), tag, MPI_COMM_WORLD, send_request);
    }
    inc_instruction_pointer();
}

void FMM::moments_send_wait(int) {
    int flag;
    if (get_level() != 0) {
        MPI_Test(send_request, &flag, MPI_STATUS_IGNORE);
        if (flag) {
            delete[] moment_buffer;
            inc_instruction_pointer();
        }
    } else {
        inc_instruction_pointer();
    }
}
void FMM::moments_send_wait_dot(int) {
    int flag;
    if (get_level() != 0) {
        MPI_Test(send_request, &flag, MPI_STATUS_IGNORE);
        if (flag) {
            delete[] moment_dot_buffer;
            inc_instruction_pointer();
        }
    } else {
        inc_instruction_pointer();
    }
}

void FMM::moments_communicate_all(int) {
    int dir, tag_send, tag_recv;
    OctNode* neighbor;
    find_neighbors();
    for (dir = 0; dir < 26; dir++) {
        neighbor = neighbors[dir];
        if (neighbor == NULL) {
            recv_request[dir] = MPI_REQUEST_NULL;
        } else {
            assert(dir==face_id[dir]);
            tag_send = tag_gen(TAG_FMM_BOUND, neighbor->get_id(), face_opp_id[dir]);
            tag_recv = tag_gen(TAG_FMM_BOUND, get_id(), face_id[dir]);
            MPI_Request req;
            MPI_Isend(poles.ptr(), 1, MPI_send_bnd_t[dir], neighbor->proc(), tag_send, MPI_COMM_WORLD, &req);
            MPI_Irecv(poles.ptr(), 1, MPI_recv_bnd_t[dir], neighbor->proc(), tag_recv, MPI_COMM_WORLD, recv_request + dir);
            MPI_Request_free(&req);
        }
    }
    inc_instruction_pointer();
}

void FMM::moments_communicate_all_dot(int) {
    int dir, tag_send, tag_recv;
    OctNode* neighbor;
    find_neighbors();
    for (dir = 0; dir < 26; dir++) {
        neighbor = neighbors[dir];
        if (neighbor == NULL) {
            recv_request[dir] = MPI_REQUEST_NULL;
        } else {
            assert(dir==face_id[dir]);
            tag_send = tag_gen(TAG_FMM_BOUND, neighbor->get_id(), face_opp_id[dir]);
            tag_recv = tag_gen(TAG_FMM_BOUND, get_id(), face_id[dir]);
            MPI_Request req;
            MPI_Isend(poles_dot.ptr(), 1, MPI_send_bnd_dot_t[dir], neighbor->proc(), tag_send, MPI_COMM_WORLD, &req);
            MPI_Irecv(poles_dot.ptr(), 1, MPI_recv_bnd_dot_t[dir], neighbor->proc(), tag_recv, MPI_COMM_WORLD, recv_request + dir);
            MPI_Request_free(&req);
        }
    }
    inc_instruction_pointer();
}

void FMM::moments_communicate_wait_all(int) {
    int flag_recv;
    bool ready;
    MPI_Testall(26, recv_request, &flag_recv, MPI_STATUS_IGNORE);
    ready = flag_recv;
    if (ready) {
        inc_instruction_pointer();
    }
}

void FMM::compute_interactions(int) {
    int xlb, xub, ylb, yub, zlb, zub;
    const Real dxinv = 1.0 / get_dx();
    multipole_t* n1;
    taylor_t *l, *l1;
//#pragma omp for collapse(2)
    for (int k0 = FBW; k0 < FNX - FBW; k0++) {
        for (int j0 = FBW; j0 < FNX - FBW; j0++) {
            for (int i0 = FBW; i0 < FNX - FBW; i0++) {
                l = L.ptr(i0, j0, k0);
                l->phi = 0.0;
                l->X = poles(i0, j0, k0).X;
            }
        }
    }
    for (int k0 = FBW; k0 < FNX - FBW; k0++) {
        for (int j0 = FBW; j0 < FNX - FBW; j0++) {
            for (int i0 = FBW; i0 < FNX - FBW; i0++) {
                if (get_level() != 0) {
                    xlb = 2 * ((i0 / 2) - FORDER);
                    ylb = 2 * ((j0 / 2) - FORDER);
                    zlb = 2 * ((k0 / 2) - FORDER);
                    xub = 2 * ((i0 / 2) + FORDER) + 1;
                    yub = 2 * ((j0 / 2) + FORDER) + 1;
                    zub = 2 * ((k0 / 2) + FORDER) + 1;
                } else {
                    xlb = ylb = zlb = FBW;
                    xub = yub = zub = FNX - FBW - 1;
                }
                n1 = poles.ptr(i0, j0, k0);
                l1 = L.ptr(i0, j0, k0);
                if (n1->is_leaf) {
                    const Real factor = (3.0 + 3.0 / sqrt(2.0) + 1.0 / sqrt(3.0)) / 3.0;
                    l1->phi() -= factor * n1->M() / get_dx();
                }
                for (int k = zlb; k <= zub; k++) {
                    const int k1 = k - k0 + INX;
                    for (int j = ylb; j <= yub; j++) {
                        const int j1 = j - j0 + INX;
                        l1 = L.ptr(i0, j0, k0);
                        for (int i = xlb; i <= xub; i++) {
                            multipole_t* n2 = poles.ptr(i, j, k);
                            taylor_t* l2 = L.ptr(i, j, k);
                            bool interaction_pair;
                            if (n1->is_leaf || n2->is_leaf) {
                                interaction_pair = !((i == i0) && (j == j0) && (k == k0));
                            } else {
                                interaction_pair = !((abs(i0 - i) <= FORDER) && (abs(j - j0) <= FORDER) && (abs(k - k0) <= FORDER));
                            }
                            if (interaction_pair) {
                                if (!n1->is_leaf || !n2->is_leaf) {
                                    const _3Vec Y = n1->X - n2->X;
                                    expansion_t D;
                                    D.compute_D(Y);
                                    l1->phi += n2->M * D;
                                } else {
                                    const int i1 = i - i0 + INX;
                                    const Real* d1 = d1_array[k1][j1][i1];
                                    Real m2 = n2->M() * dxinv;
                                    Real m2dxinv = m2 * dxinv;
                                    l1->phi() += d0_array[abs(k - k0)][abs(j - j0)][abs(i - i0)] * m2;
                                    for (int i = 0; i < 3; i++) {
                                        l1->phi(i) += d1[i] * m2dxinv;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    for (int k0 = FBW; k0 < FNX - FBW; k0++) {
        for (int j0 = FBW; j0 < FNX - FBW; j0++) {
            for (int i0 = FBW; i0 < FNX - FBW; i0++) {
                l = L.ptr(i0, j0, k0);
                n1 = poles.ptr(i0, j0, k0);
                l->phi *= PhysicalConstants::G;
                l->g_lz = 0.0;
                const _3Vec R = n1->X;
                for (int a = 0; a < 2; a++) {
                    const int b = 1 - a;
                    const Real eps = Real(2 * b - 1);
                    l->g_lz() -= eps * R[a] * l->phi(b);
                    for (int m = 0; m < 3; m++) {
                        l->g_lz(m) -= eps * (+delta[a][m] * l->phi(b) + R[a] * l->phi(b, m));
                        for (int n = 0; n <= m; n++) {
                            l->g_lz(m, n) -= eps * (delta[a][n] * l->phi(b, m) + delta[a][m] * l->phi(b, n) + R[a] * l->phi(b, n, m));
                        }
                    }
                }
            }
        }
    }
    inc_instruction_pointer();
}

void FMM::compute_interactions_dot(int) {
    int xlb, xub, ylb, yub, zlb, zub;
    const Real dxinv = 1.0 / get_dx();
    multipole_dot_t* n1;
    taylor_dot_t *l, *l1;
    for (int k0 = FBW; k0 < FNX - FBW; k0++) {
        for (int j0 = FBW; j0 < FNX - FBW; j0++) {
            for (int i0 = FBW; i0 < FNX - FBW; i0++) {
                l = L_dot.ptr(i0, j0, k0);
                l->phi_dot = 0.0;
            }
        }
    }
    for (int k0 = FBW; k0 < FNX - FBW; k0++) {
        for (int j0 = FBW; j0 < FNX - FBW; j0++) {
            for (int i0 = FBW; i0 < FNX - FBW; i0++) {
                if (get_level() != 0) {
                    xlb = 2 * ((i0 / 2) - FORDER);
                    ylb = 2 * ((j0 / 2) - FORDER);
                    zlb = 2 * ((k0 / 2) - FORDER);
                    xub = 2 * ((i0 / 2) + FORDER) + 1;
                    yub = 2 * ((j0 / 2) + FORDER) + 1;
                    zub = 2 * ((k0 / 2) + FORDER) + 1;
                } else {
                    xlb = ylb = zlb = FBW;
                    xub = yub = zub = FNX - FBW - 1;
                }
                n1 = poles_dot.ptr(i0, j0, k0);
                l1 = L_dot.ptr(i0, j0, k0);
                if (poles(i0, j0, k0).is_leaf) {
                    const Real factor = (3.0 + 3.0 / sqrt(2.0) + 1.0 / sqrt(3.0)) / 3.0;
                    l1->phi_dot() -= factor * n1->M_dot() / get_dx();
                }
                for (int k = zlb; k <= zub; k++) {
                    const int k1 = k - k0 + INX;
                    for (int j = ylb; j <= yub; j++) {
                        const int j1 = j - j0 + INX;
                        l1 = L_dot.ptr(i0, j0, k0);
                        for (int i = xlb; i <= xub; i++) {
                            multipole_dot_t* n2 = poles_dot.ptr(i, j, k);
                            taylor_dot_t* l2 = L_dot.ptr(i, j, k);
                            bool interaction_pair;
                            if (poles(i, j, k).is_leaf || poles(i0, j0, k0).is_leaf) {
                                interaction_pair = !((i == i0) && (j == j0) && (k == k0));
                            } else {
                                interaction_pair = !((abs(i0 - i) <= FORDER) && (abs(j - j0) <= FORDER) && (abs(k - k0) <= FORDER));
                            }
                            if (interaction_pair) {
                                if (!poles(i, j, k).is_leaf || !poles(i0, j0, k0).is_leaf) {
                                    const _3Vec Y = poles(i0, j0, k0).X - poles(i, j, k).X;
                                    expansion_t D;
                                    D.compute_D(Y);
                                    l1->phi_dot += n2->M_dot * D;
                                } else {
                                    const int i1 = i - i0 + INX;
                                    const Real d0 = d0_array[abs(k - k0)][abs(j - j0)][abs(i - i0)];
                                    Real dm2 = n2->M_dot() * dxinv;
                                    l1->phi_dot() += d0 * dm2;
                                    Real dm1 = n1->M_dot() * dxinv;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    for (int k0 = FBW; k0 < FNX - FBW; k0++) {
        for (int j0 = FBW; j0 < FNX - FBW; j0++) {
            for (int i0 = FBW; i0 < FNX - FBW; i0++) {
                l = L_dot.ptr(i0, j0, k0);
                l->phi_dot *= PhysicalConstants::G;
            }
        }
    }
    inc_instruction_pointer();
}

void FMM::expansion_recv(int) {
    if (get_level() != 0) {
        const ChildIndex ci = my_child_index();
        const FMM* p = dynamic_cast<const FMM*>(get_parent());
        int i, j, k, sz, tag;
        Real a;
        tag = tag_gen(TAG_FMM_EXPANSION, get_id(), ci);
        taylor_buffer = new taylor_t[INX * INX * INX / 8];
        MPI_Irecv(taylor_buffer, INX * INX * INX * sizeof(taylor_t) / 8, MPI_BYTE, p->proc(), tag, MPI_COMM_WORLD, recv_request + 0);
    }
    inc_instruction_pointer();

}

void FMM::expansion_recv_wait(int) {
    bool rc;
    int cnt;
    if (get_level() != 0) {
        int flag;
        MPI_Test(recv_request, &flag, MPI_STATUS_IGNORE);
        rc = flag;
        if (rc) {
            cnt = 0;
            for (int k0 = FBW; k0 < FNX - FBW; k0 += 2) {
                for (int j0 = FBW; j0 < FNX - FBW; j0 += 2) {
                    for (int i0 = FBW; i0 < FNX - FBW; i0 += 2) {
                        taylor_t pL;
                        pL.phi = taylor_buffer[cnt].phi;
                        pL.g_lz = taylor_buffer[cnt].g_lz;
                        pL.X = taylor_buffer[cnt].X;
                        cnt++;
//#pragma omp parallel for collapse(2)
                        for (int k = k0; k < k0 + 2; k++) {
                            for (int j = j0; j < j0 + 2; j++) {
                                for (int i = i0; i < i0 + 2; i++) {
                                    _3Vec dX = L(i, j, k).X - pL.X;
                                    L(i, j, k).phi += pL.phi << dX;
                                    L(i, j, k).g_lz += pL.g_lz << dX;
                                }
                            }
                        }
                    }
                }
            }
            delete[] taylor_buffer;
        }
    } else {
        rc = true;
    }
    if (rc) {
        inc_instruction_pointer();
    }

}

void FMM::expansion_send(int) {
    FMM* child;
    int tag;
    for (int ci = 0; ci < 8; ci++) {
        child = dynamic_cast<FMM*>(get_child(ci));
        if (child == NULL) {
            send_request[ci] = MPI_REQUEST_NULL;
//#pragma omp parallel for  collapse(2)
            for (int k = FBW; k < FNX - FBW; k++) {
                for (int j = FBW; j < FNX - FBW; j++) {
                    for (int i = FBW; i < FNX - FBW; i++) {
                        _4force(i, j, k).phi = L(i, j, k).phi();
                        for (int a = 0; a < 3; a++) {
                            _4force(i, j, k).g[a] = -L(i, j, k).phi(a);
                        }
                    }
                }
            }
        } else {
            tag = tag_gen(TAG_FMM_EXPANSION, child->get_id(), ci);
            MPI_Request req;
            MPI_Isend(L.ptr(), 1, MPI_comm_taylor_t[ci], child->proc(), tag, MPI_COMM_WORLD, &req);
            MPI_Request_free(&req);
        }
    }
    inc_instruction_pointer();

}

void FMM::expansion_recv_dot(int) {
    if (get_level() != 0) {
        const ChildIndex ci = my_child_index();
        const FMM* p = dynamic_cast<const FMM*>(get_parent());
        int i, j, k, sz, tag;
        Real a;
        tag = tag_gen(TAG_FMM_EXPANSION, get_id(), ci);
        taylor_dot_buffer = new taylor_dot_t[INX * INX * INX / 8];
        MPI_Irecv(taylor_dot_buffer, INX * INX * INX * sizeof(taylor_dot_t) / 8, MPI_BYTE, p->proc(), tag, MPI_COMM_WORLD, recv_request + 0);
    }
    inc_instruction_pointer();

}

void FMM::expansion_recv_wait_dot(int) {
    bool rc;
    int cnt;
    if (get_level() != 0) {
        int flag;
        MPI_Test(recv_request, &flag, MPI_STATUS_IGNORE);
        rc = flag;
        if (rc) {
            cnt = 0;
            for (int k0 = FBW; k0 < FNX - FBW; k0 += 2) {
                for (int j0 = FBW; j0 < FNX - FBW; j0 += 2) {
                    for (int i0 = FBW; i0 < FNX - FBW; i0 += 2) {
                        taylor_dot_t pL;
                        pL.phi_dot = taylor_dot_buffer[cnt].phi_dot;
                        cnt++;
//#pragma omp parallel for collapse(2)
                        for (int k = k0; k < k0 + 2; k++) {
                            for (int j = j0; j < j0 + 2; j++) {
                                for (int i = i0; i < i0 + 2; i++) {
                                    _3Vec dX = L(i, j, k).X - Xp(i0 / 2, j0 / 2, k0 / 2);
                                    L_dot(i, j, k).phi_dot += pL.phi_dot << dX;
                                }
                            }
                        }
                    }
                }
            }
            delete[] taylor_dot_buffer;
        }
    } else {
        rc = true;
    }
    if (rc) {
        inc_instruction_pointer();
    }

}

void FMM::expansion_send_dot(int) {
    FMM* child;
    int tag;
    for (int ci = 0; ci < 8; ci++) {
        child = dynamic_cast<FMM*>(get_child(ci));
        if (child != NULL) {
            tag = tag_gen(TAG_FMM_EXPANSION, child->get_id(), ci);
            MPI_Request req;
            MPI_Isend(L_dot.ptr(), 1, MPI_comm_taylor_dot_t[ci], child->proc(), tag, MPI_COMM_WORLD, &req);
            MPI_Request_free(&req);
        }
    }
    inc_instruction_pointer();

}

void FMM::expansion_send_wait(int) {
    inc_instruction_pointer();
}

void FMM::FMM_solve() {
//	printf( "Solve\n");
//	printf("!\n");
    MPI_datatypes_init();
    FMM** list = new FMM*[get_local_node_cnt()];
    for (int i = 0; i < get_local_node_cnt(); i++) {
        list[i] = dynamic_cast<FMM*>(get_local_node(i));
        list[i]->ip = 0;
    }
    int last_ip;
    bool done;
    const int nprocs = get_local_node_cnt();

    do {
        done = true;
        for (int i = 0; i < nprocs; i++) {
            if (list[i]->ip < FSTAGE) {
                //			do {
                //	printf("%i %i %i\n", list[i]->get_id(), list[i]->get_level(), list[i]->ip[0]);
                last_ip = list[i]->ip;
                (list[i]->*cs[list[i]->ip])(0);
                //		} while (last_ip != list[i]->ip[0] && list[i]->ip[0] <= FSTAGE);
                if (list[i]->ip <= FSTAGE) {
                    done = false;
                }
            }
        }
    } while (!done);

    delete[] list;
    pot_to_hydro_grid();
}

Real FMM::substep_driver() {
    // MPI_rank() ? 0 : printf("Starting substep driver\n");
    DFO = Vector<Real, STATE_NF>(0.0);
    MPI_datatypes_init();
    FMM** list = new FMM*[get_local_node_cnt()];
    for (int i = 0; i < get_local_node_cnt(); i++) {
        list[i] = dynamic_cast<FMM*>(get_local_node(i));
        list[i]->ip = 0;
    }
    int last_ip;
    bool done;
    const int nprocs = get_local_node_cnt();

    if (_beta == 1.0) {
        do {
            done = true;
            for (int i = 0; i < nprocs; i++) {
                if (list[i]->ip != 21) {
                    (list[i]->*fmmhydro[list[i]->ip])(0);
                    if (list[i]->ip != 21) {
                        done = false;
                    }
                }
            }
        } while (!done);
        MPI_Bcast(&root_dt, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        Hydro::set_dt(root_dt);

    }
    do {
        done = true;
        for (int i = 0; i < nprocs; i++) {
            if (list[i]->ip < 40) {
                (list[i]->*fmmhydro[list[i]->ip])(0);
                if (list[i]->ip < 40) {
                    done = false;
                }
            }
        }
    } while (!done);
    delete[] list;
    Real tmp[STATE_NF];
    for (int i = 0; i < STATE_NF; i++) {
        tmp[i] = DFO[i];
    }
    Real tmp2[STATE_NF];
    MPI_Allreduce(tmp, tmp2, STATE_NF, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD);
    for (int i = 0; i < STATE_NF; i++) {
        DFO[i] = tmp2[i];
    }
    DFO[State::et_index] += State::omega * DFO[State::lz_index] + DFO[State::pot_index];
    FO = (FO + DFO * _dt) * _beta + FO0 * (1.0 - _beta);
//    MPI_rank() ? 0 : printf("Ending substep driver\n");
    return root_dt;
}

void FMM::FMM_from_children() {
//	printf( "Solve\n");
//	printf("!\n");
    MPI_datatypes_init();
    FMM** list = new FMM*[get_local_node_cnt()];
    for (int i = 0; i < get_local_node_cnt(); i++) {
        list[i] = dynamic_cast<FMM*>(get_local_node(i));
        list[i]->ip = 0;
    }
    int last_ip;
    bool done;
    const int nprocs = get_local_node_cnt();

    do {
        done = true;
        for (int i = 0; i < nprocs; i++) {
            if (list[i]->ip < 4) {
                last_ip = list[i]->ip;
                (list[i]->*cs_children[list[i]->ip])(0);
                if (list[i]->ip <= 4) {
                    done = false;
                }
            }
        }
    } while (!done);

    delete[] list;
    pot_to_hydro_grid();
}

void FMM::FMM_solve_dot() {
//	printf( "Solve_dot\n");
//	printf("!\n");
    MPI_datatypes_init();
    FMM** list = new FMM*[get_local_node_cnt()];
    for (int i = 0; i < get_local_node_cnt(); i++) {
        list[i] = dynamic_cast<FMM*>(get_local_node(i));
        list[i]->ip = 0;
    }
    int last_ip;
    bool done;
    const int nprocs = get_local_node_cnt();

    do {
        done = true;
        for (int i = 0; i < nprocs; i++) {
            if (list[i]->ip < FSTAGE) {
                //			do {
                //	printf("%i %i %i\n", list[i]->get_id(), list[i]->get_level(), list[i]->ip[0]);
                last_ip = list[i]->ip;
                (list[i]->*cs_dot[list[i]->ip])(0);
                //		} while (last_ip != list[i]->ip[0] && list[i]->ip[0] <= FSTAGE);
                if (list[i]->ip <= FSTAGE) {
                    done = false;
                }
            }
        }
    } while (!done);

    delete[] list;
}

FMM::FMM() {
    MPI_datatypes_init();

}

void FMM::initialize() {
}

FMM* FMM::new_octnode() const {
    return new FMM;
}

void FMM::allocate_arrays() {
    Hydro::allocate_arrays();
    Xp.allocate();
    poles.allocate();
    poles_dot.allocate();
    L.allocate();
    L_dot.allocate();
    _4force.allocate();
    old_pot.allocate();
    dpot.allocate();
}

void FMM::deallocate_arrays() {
    Hydro::deallocate_arrays();
    Xp.deallocate();
    poles.deallocate();
    L.deallocate();
    poles_dot.deallocate();
    L_dot.deallocate();
    _4force.deallocate();
    old_pot.deallocate();
    dpot.deallocate();
}

FMM::~FMM() {
// TODO Auto-generated destructor stub
}
void FMM::_4force_recv(int) {
    int tag;
    FMM* child;
    const Real dv = pow(get_dx(), 3);
    for (ChildIndex ci = 0; ci < 8; ci++) {
        child = dynamic_cast<FMM*>(get_child(ci));

        if (child != NULL) {
            tag = tag_gen(TAG_FMM_4FORCE, get_id(), ci);
            MPI_Irecv(_4force.ptr(), 1, MPI_comm_child3_t[ci], get_child(ci)->proc(), tag, MPI_COMM_WORLD, recv_request + ci);
        } else {
            recv_request[ci] = MPI_REQUEST_NULL;
        }
    }
    inc_instruction_pointer();
}

void FMM::_4force_recv_wait(int) {
    int flag;
    MPI_Testall(8, recv_request, &flag, MPI_STATUS_IGNORE);
    if (flag) {
        inc_instruction_pointer();
    }
}

void FMM::_4force_send(int) {
    if (get_level() != 0) {
        int cnt;
        _3Vec Y;
        _4force_t tmp;
        cnt = 0;
        _4force_buffer = new _4force_t[INX * INX * INX / 8];
        for (int k = FBW; k < FNX - FBW; k += 2) {
            for (int j = FBW; j < FNX - FBW; j += 2) {
                for (int i = FBW; i < FNX - FBW; i += 2) {
                    tmp.phi = 0.0;
                    tmp.g = 0.0;
                    for (int a = i; a < i + 2; a++) {
                        for (int b = j; b < j + 2; b++) {
                            for (int c = k; c < k + 2; c++) {
                                tmp.phi += _4force(a, b, c).phi;
                                tmp.g += _4force(a, b, c).g;
                            }
                        }
                    }
                    tmp.phi *= 0.125;
                    tmp.g *= 0.125;
                    _4force_buffer[cnt] = tmp;
                    cnt++;
                }
            }
        }
        assert(cnt==INX * INX * INX / 8);
        int tag = tag_gen(TAG_FMM_4FORCE, get_parent()->get_id(), my_child_index());
        MPI_Isend(_4force_buffer, cnt * sizeof(_4force_t), MPI_BYTE, get_parent()->proc(), tag, MPI_COMM_WORLD, send_request);
    }
    inc_instruction_pointer();
}
void FMM::_4force_send_wait(int) {
    int flag;
    if (get_level() != 0) {
        MPI_Test(send_request, &flag, MPI_STATUS_IGNORE);
        if (flag) {
            delete[] _4force_buffer;
            inc_instruction_pointer();
        }
    } else {
        inc_instruction_pointer();
    }
}

Real FMM::g_energy(int ii, int jj, int kk) const {

    const int i = ii + FBW - BW;
    const int j = jj + FBW - BW;
    const int k = kk + FBW - BW;
    const Real dvinv = 1.0 / pow(get_dx(), 3);
//	return _4force(i, j, k).phi_dot * poles(i, j, k).M()*dvinv;
    Real term1 = L(i, j, k).phi() * poles_dot(i, j, k).M_dot();
    Real term2 = L_dot(i, j, k).phi_dot() * poles(i, j, k).M();
    return (term2 - term1) * 0.5 * dvinv;
}

Real FMM::set_phi(int i, int j, int k, Real phi) {
    return _4force(i + FBW - BW, j + FBW - BW, k + FBW - BW).phi = phi;
}

Real FMM::get_phi(int i, int j, int k) const {
    return _4force(i + FBW - BW, j + FBW - BW, k + FBW - BW).phi;
}

Real FMM::gx(int i, int j, int k) const {
    return _4force(i + FBW - BW, j + FBW - BW, k + FBW - BW).g[0];
}

Real FMM::gy(int i, int j, int k) const {
    return _4force(i + FBW - BW, j + FBW - BW, k + FBW - BW).g[1];
}

Real FMM::gz(int i, int j, int k) const {
    return _4force(i + FBW - BW, j + FBW - BW, k + FBW - BW).g[2];
}

void FMM::null(int) {
    inc_instruction_pointer();
}

#define RK2

void FMM::step(Real dt) {
    Real start_time;
#ifdef RK2
    int nsteps = 2;
    Real beta[2] = { 1.0, 1.0 / 2.0 };
#else
    int nsteps = 3;
    Real beta[3] = {1.0, 0.25, 2.0 / 3.0};
#endif
    Hydro::set_dt(dt);
    store();
    store_pot();
    for (int i = 0; i < nsteps; i++) {

        Hydro::set_beta(beta[i]);
        start_time = MPI_Wtime();
        substep_driver();
        FMM_solve_dot();
        update();
        FMM_solve();
        account_pot();
    }
    FMM_from_children();
    pot_to_hydro_grid();
    set_time(get_time() + dt);
}

void FMM::compute_update(int) {
    Hydro::inc_instruction_pointer();
}

bool FMM::check_for_refine() {
    bool rc;
    to_conserved_energy();
    inject_from_children();
    rc = OctNode::check_for_refine();
    if (rc) {
        pot_to_hydro_grid();
        Hydro::redistribute_grids();
        FMM_solve();
        FMM_from_children();
    }
    from_conserved_energy();

    if (dynamic_cast<Hydro*>(get_root())->check_for_expand()) {
        Hydro::redistribute_grids();
//        OctNode::check_for_refine();
        inject_from_children();
        FMM_solve();
    }

    return rc;

}

void FMM::store_pot() {
    FMM* g0;
    for (int n = 0; n < get_local_node_cnt(); n++) {
        g0 = dynamic_cast<FMM*>(get_local_node(n));
//#pragma omp parallel for collapse(2)
        for (int k = BW; k < GNX - BW; k++) {
            for (int j = BW; j < GNX - BW; j++) {
                for (int i = BW; i < GNX - BW; i++) {
                    g0->dpot(i, j, k) = 0.0;
                    Real pot = g0->get_phi(i, j, k);
                    pot *= (*g0)(i, j, k).rho();
                    pot += (*g0)(i, j, k).rot_pot(g0->X(i, j, k));
                    g0->old_pot(i, j, k) = pot - 0.5 * g0->get_phi(i, j, k) * (*g0)(i, j, k).rho();
                }
            }
        }
    }
}

void FMM::apply_pot_change(int) {
//#pragma omp parallel for collapse(2)

    for (int k = BW; k < GNX - BW; k++) {
        for (int j = BW; j < GNX - BW; j++) {
            for (int i = BW; i < GNX - BW; i++) {
                Real pot = get_phi(i, j, k);
                pot *= U(i, j, k).rho();
                pot += U(i, j, k).rot_pot(X(i, j, k));
                U(i, j, k).set_pot(pot);
                Real new_pot = U(i, j, k).pot() - 0.5 * get_phi(i, j, k) * U(i, j, k).rho();
                dpot(i, j, k) = (new_pot - old_pot(i, j, k));
                U(i, j, k)[State::et_index] -= (new_pot - old_pot(i, j, k));
            }
        }
    }
    inc_instruction_pointer();
}

void FMM::account_pot() {
    FMM* g0;
    for (int n = 0; n < get_local_node_cnt(); n++) {
        g0 = dynamic_cast<FMM*>(get_local_node(n));
//#pragma omp parallel for collapse(2)
        for (int k = BW; k < GNX - BW; k++) {
            for (int j = BW; j < GNX - BW; j++) {
                for (int i = BW; i < GNX - BW; i++) {
                    Real pot = g0->get_phi(i, j, k);
                    pot *= (*g0)(i, j, k).rho();
                    pot += (*g0)(i, j, k).rot_pot(g0->X(i, j, k));
                    (*g0)(i, j, k).set_pot(pot);
                    Real new_pot = (*g0)(i, j, k).pot() - 0.5 * g0->get_phi(i, j, k) * (*g0)(i, j, k).rho();
                    g0->dpot(i, j, k) = (new_pot - g0->old_pot(i, j, k));
                    (*g0)(i, j, k)[State::et_index] -= (new_pot - g0->old_pot(i, j, k));
                }
            }
        }
    }
}

void FMM::update() {
    FMM* g0;
    for (int n = 0; n < get_local_node_cnt(); n++) {
        g0 = dynamic_cast<FMM*>(get_local_node(n));
//#pragma omp parallel for collapse(2)
        for (int k = BW; k < GNX - BW; k++) {
            for (int j = BW; j < GNX - BW; j++) {
                for (int i = BW; i < GNX - BW; i++) {
                    _3Vec x = g0->X(i, j, k);
                    Real d = (*g0)(i, j, k).rho() - State::rho_floor;
                    Real R = sqrt(x[0] * x[0] + x[1] * x[1]);
                    g0->D(i, j, k)[State::lz_index] += d * g0->g_lz(i, j, k);
                    g0->D(i, j, k)[State::sx_index] += d * g0->gx(i, j, k);
                    g0->D(i, j, k)[State::sy_index] += d * g0->gy(i, j, k);
                    g0->D(i, j, k)[State::sz_index] += d * g0->gz(i, j, k);
                    g0->D(i, j, k)[State::et_index] += g0->g_energy(i, j, k) + g0->D(i, j, k)[State::pot_index];
                    g0->D(i, j, k)[State::pot_index] = 0.0;
                }
            }
        }
//#pragma omp parallel for collapse(2)
        for (int k = BW; k < GNX - BW; k++) {
            for (int j = BW; j < GNX - BW; j++) {
                for (int i = BW; i < GNX - BW; i++) {
                    g0->U(i, j, k)[State::et_index] += g0->dpot(i, j, k);
                    g0->U(i, j, k) = (g0->U(i, j, k) + g0->D(i, j, k) * _dt) * _beta + g0->U0(i, j, k) * (1.0 - _beta);
                    g0->U(i, j, k).floor(g0->X(i, j, k), g0->get_dx());
                }
            }
        }
    }
}

void FMM::compute_gravity_update(int) {
//#pragma omp parallel for collapse(2)
    for (int k = BW; k < GNX - BW; k++) {
        for (int j = BW; j < GNX - BW; j++) {
            for (int i = BW; i < GNX - BW; i++) {
                _3Vec x = X(i, j, k);
                Real d = U(i, j, k).rho() - State::rho_floor;
                Real R = sqrt(x[0] * x[0] + x[1] * x[1]);
                D(i, j, k)[State::lz_index] += d * g_lz(i, j, k);
                D(i, j, k)[State::sx_index] += d * gx(i, j, k);
                D(i, j, k)[State::sy_index] += d * gy(i, j, k);
                D(i, j, k)[State::sz_index] += d * gz(i, j, k);
                D(i, j, k)[State::et_index] += g_energy(i, j, k) + D(i, j, k)[State::pot_index];
                D(i, j, k)[State::pot_index] = 0.0;
            }
        }
    }

//#pragma omp parallel for collapse(2)
    for (int k = BW; k < GNX - BW; k++) {
        for (int j = BW; j < GNX - BW; j++) {
            for (int i = BW; i < GNX - BW; i++) {
                U(i, j, k)[State::et_index] += dpot(i, j, k);
                U(i, j, k) = (U(i, j, k) + D(i, j, k) * _dt) * _beta + U0(i, j, k) * (1.0 - _beta);
                U(i, j, k).floor(X(i, j, k), get_dx());
            }
        }
    }
    inc_instruction_pointer();
}

_3Vec FMM::system_com() {
    FMM* g0;
    Real send[4], recv[4];
    _3Vec rc;
    Real mtot, mx, my, mz;
    mtot = mx = my = mz = 0.0;
    for (int n = 0; n < get_local_node_cnt(); n++) {
        g0 = dynamic_cast<FMM*>(get_local_node(n));
        if (g0->get_level() == 0) {
//#pragma omp parallel for collapse(2) reduction(+:mx,my,mz,mtot)
            for (int k = BW; k < GNX - BW; k++) {
                for (int j = BW; j < GNX - BW; j++) {
                    for (int i = BW; i < GNX - BW; i++) {
                        Real m = g0->poles(i, j, k).M();
                        mx += g0->poles(i, j, k).X[0] * m;
                        my += g0->poles(i, j, k).X[1] * m;
                        mz += g0->poles(i, j, k).X[2] * m;
                        mtot += m;
                    }
                }
            }
        }
    }
    send[0] = mtot;
    send[1] = mx;
    send[2] = my;
    send[3] = mz;
    MPI_Allreduce(send, recv, 4, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD);
    mtot = recv[0];
    mx = recv[1];
    my = recv[2];
    mz = recv[3];
    rc[0] = mx / mtot;
    rc[1] = my / mtot;
    rc[2] = mz / mtot;
    return rc;
}

void FMM::to_conserved_energy() {
    Hydro* g0;
    for (int n = 0; n < get_local_node_cnt(); n++) {
        g0 = dynamic_cast<Hydro*>(get_local_node(n));
//#pragma omp parallel for collapse(2)
        for (int k = BW - 1; k < GNX - BW + 1; k++) {
            for (int j = BW - 1; j < GNX - BW + 1; j++) {
                for (int i = BW - 1; i < GNX - BW + 1; i++) {
                    (*g0)(i, j, k).to_con(g0->X(i, j, k));
                }
            }
        }
    }
}

void FMM::from_conserved_energy() {
    Hydro* g0;
    for (int n = 0; n < get_local_node_cnt(); n++) {
        g0 = dynamic_cast<Hydro*>(get_local_node(n));
//#pragma omp parallel for collapse(2)
        for (int k = BW - 1; k < GNX - BW + 1; k++) {
            for (int j = BW - 1; j < GNX - BW + 1; j++) {
                for (int i = BW - 1; i < GNX - BW + 1; i++) {
                    (*g0)(i, j, k).from_con(g0->X(i, j, k));
                }
            }
        }
    }
}

Vector<Real, 6> FMM::momentum_sum() {
    FMM* g0;
    Real sum[3];
    Real lz, lz2, sx2, sy2, R, x, y;
    lz = 0.0;
    lz2 = 0.0;
    sx2 = sy2 = 0.0;
    sum[0] = sum[1] = sum[2] = 0.0;
    Real norm = 0.0;
    Real energy = 0.0;
    Real norme = 0.0;
    Real norml = 0.0, this_sr, this_lz;
    for (int n = 0; n < get_local_node_cnt(); n++) {
        g0 = dynamic_cast<FMM*>(get_local_node(n));
        for (int k = FBW; k < FNX - FBW; k++) {
            for (int j = FBW; j < FNX - FBW; j++) {
                for (int i = FBW; i < FNX - FBW; i++) {
                    if (g0->poles(i, j, k).is_leaf) {
                        Real de = g0->g_energy(i + BW - FBW, j + BW - FBW, k + BW - FBW) * g0->get_dx() * g0->get_dx() * g0->get_dx();
                        energy += de;
                        norme += fabs(de);
                        norml += fabs(g0->g_lz(i + BW - FBW, j + BW - FBW, k + BW - FBW) * g0->poles(i, j, k).M());
                        lz2 += g0->g_lz(i + BW - FBW, j + BW - FBW, k + BW - FBW) * g0->poles(i, j, k).M();
                        lz += (g0->_4force(i, j, k).g[1] * g0->poles(i, j, k).X[0] - g0->_4force(i, j, k).g[0] * g0->poles(i, j, k).X[1])
                                * g0->poles(i, j, k).M();
                        sum[0] += g0->_4force(i, j, k).g[0] * g0->poles(i, j, k).M();
                        sum[1] += g0->_4force(i, j, k).g[1] * g0->poles(i, j, k).M();
                        sum[2] += g0->_4force(i, j, k).g[2] * g0->poles(i, j, k).M();
                        norm += sqrt((g0->_4force(i, j, k).g).dot(g0->_4force(i, j, k).g)) * g0->poles(i, j, k).M();
                    }
                }
            }
        }
    }
    Real sum0[3];
    Real norm0 = norm;
    Real lz0 = lz;
    sum0[0] = sum[0];
    sum0[1] = sum[1];
    sum0[2] = sum[2];
    MPI_Allreduce(sum0, sum, 3, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&norm0, &norm, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&lz0, &lz, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD);
    lz0 = energy;
    MPI_Allreduce(&lz0, &energy, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD);
    lz0 = norme;
    MPI_Allreduce(&lz0, &norme, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD);
    lz0 = sx2;
    MPI_Allreduce(&lz0, &sx2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD);
    lz0 = sy2;
    MPI_Allreduce(&lz0, &sy2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD);
    lz0 = lz2;
    MPI_Allreduce(&lz0, &lz2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD);
    norm0 = norml;
    MPI_Allreduce(&norm0, &norml, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD);
    if (MPI_rank() == 0 && energy != 0.0) {
        printf("Momentum Sum = %e %e %e %e %e %e %e %e\n", sum[0] / norm, sum[1] / norm, sum[2] / norm, lz / norml, lz2 / norml, energy, norme, energy / norme);
    }
    Vector<Real, 6> ret_vec;
    ret_vec[0] = lz / norml;
    ret_vec[1] = lz2 / norml;
    ret_vec[2] = sum[0] / norm;
    ret_vec[3] = sx2 / norm;
    ret_vec[4] = sum[1] / norm;
    ret_vec[5] = sy2 / norm;
    return ret_vec;
}

Vector<Real, 4> FMM::com_sum() {
    FMM* g0;
    Real sum[4];
    sum[0] = sum[1] = sum[2] = sum[3] = 0.0;
    Real norm = 0.0;
    const Real floor = State::rho_floor;
    Real d;
    for (int n = 0; n < get_local_node_cnt(); n++) {
        g0 = dynamic_cast<FMM*>(get_local_node(n));
        const Real dv = pow(g0->get_dx(), 3);
        for (int k = BW; k < GNX - BW; k++) {
            for (int j = BW; j < GNX - BW; j++) {
                for (int i = BW; i < GNX - BW; i++) {
                    if (!g0->zone_is_refined(i, j, k)) {
                        d = (*g0)(i, j, k).rho() - floor;
                        sum[0] += d * dv * g0->xc(i);
                        sum[1] += d * dv * g0->yc(j);
                        sum[2] += d * dv * g0->zc(k);
                        sum[3] += d * dv;
                    }
                }
            }
        }
    }
    Real sum0[4];
    sum0[0] = sum[0];
    sum0[1] = sum[1];
    sum0[2] = sum[2];
    sum0[3] = sum[3];
    MPI_Allreduce(sum0, sum, 4, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD);
    sum[0] /= sum[3];
    sum[1] /= sum[3];
    sum[2] /= sum[3];
    if (MPI_rank() == 0) {
        FILE* fp = fopen("com.dat", "at");
        fprintf(fp, "%e %.14e %.14e %.14e  %.14e\n", get_time(), sum[0], sum[1], sum[2], sum[3]);
        fclose(fp);
    }
    Vector<Real, 4> v;
    v[0] = sum[0];
    v[1] = sum[1];
    v[2] = sum[2];
    v[3] = sum[3];
    return v;
}

bool FMM::is_leaf(int i, int j, int k) const {
    const FMM* ptr = this;
    for (int l = 0; l < 3; l++) {
        if (i < FBW) {
            if (ptr->get_sibling(XL) != NULL) {
                ptr = dynamic_cast<const FMM*>(ptr->get_sibling(XL));
                i += INX;
            }
        } else if (i >= FNX - FBW) {
            if (ptr->get_sibling(XU) != NULL) {
                ptr = dynamic_cast<const FMM*>(ptr->get_sibling(XU));
                i -= INX;
            }
        }
        if (j < FBW) {
            if (ptr->get_sibling(YL) != NULL) {
                ptr = dynamic_cast<const FMM*>(ptr->get_sibling(YL));
                j += INX;
            }
        } else if (j >= FNX - FBW) {
            if (ptr->get_sibling(YU) != NULL) {
                ptr = dynamic_cast<const FMM*>(ptr->get_sibling(YU));
                j -= INX;
            }
        }
        if (k < FBW) {
            if (ptr->get_sibling(ZL) != NULL) {
                ptr = dynamic_cast<const FMM*>(ptr->get_sibling(ZL));
                k += INX;
            }
        } else if (k >= FNX - FBW) {
            if (ptr->get_sibling(ZU) != NULL) {
                ptr = dynamic_cast<const FMM*>(ptr->get_sibling(ZU));
                k -= INX;
            }
        }
    }
    return !ptr->zone_is_refined(i + BW - FBW, j + BW - FBW, k + BW - FBW);
}

void FMM::find_neighbors() {
    OctNode* corners[8];
    OctNode* edges[12];
//#pragma omp parallel for
    for (int i = 0; i < 12; i++) {
        edges[i] = NULL;
    }
//#pragma omp parallel for
    for (int i = 0; i < 8; i++) {
        corners[i] = NULL;
    }
    if (get_sibling(XL) != NULL) {
        edges[0] = get_sibling(XL)->get_sibling(YL);
        edges[1] = get_sibling(XL)->get_sibling(YU);
        edges[2] = get_sibling(XL)->get_sibling(ZL);
        edges[3] = get_sibling(XL)->get_sibling(ZU);
    }
    if (get_sibling(XU) != NULL) {
        edges[4] = get_sibling(XU)->get_sibling(YL);
        edges[5] = get_sibling(XU)->get_sibling(YU);
        edges[6] = get_sibling(XU)->get_sibling(ZL);
        edges[7] = get_sibling(XU)->get_sibling(ZU);
    }
    if (get_sibling(YL) != NULL) {
        edges[0] = get_sibling(YL)->get_sibling(XL);
        edges[4] = get_sibling(YL)->get_sibling(XU);
        edges[8] = get_sibling(YL)->get_sibling(ZL);
        edges[9] = get_sibling(YL)->get_sibling(ZU);
    }
    if (get_sibling(YU) != NULL) {
        edges[1] = get_sibling(YU)->get_sibling(XL);
        edges[5] = get_sibling(YU)->get_sibling(XU);
        edges[10] = get_sibling(YU)->get_sibling(ZL);
        edges[11] = get_sibling(YU)->get_sibling(ZU);
    }
    if (get_sibling(ZL) != NULL) {
        edges[2] = get_sibling(ZL)->get_sibling(XL);
        edges[6] = get_sibling(ZL)->get_sibling(XU);
        edges[8] = get_sibling(ZL)->get_sibling(YL);
        edges[10] = get_sibling(ZL)->get_sibling(YU);
    }
    if (get_sibling(ZU) != NULL) {
        edges[3] = get_sibling(ZU)->get_sibling(XL);
        edges[7] = get_sibling(ZU)->get_sibling(XU);
        edges[9] = get_sibling(ZU)->get_sibling(YL);
        edges[11] = get_sibling(ZU)->get_sibling(YU);
    }
//#pragma omp parallel for
    for (int i = 0; i < 8; i++) {
        int l = (i & 1) + ZL;
        int k = ((i & 2) >> 1) + YL;
        int j = ((i & 4) >> 2) + XL;
        corners[i] = NULL;
        if (get_sibling(j) != NULL && corners[i] == NULL) {
            if (get_sibling(j)->get_sibling(k) != NULL && corners[i] == NULL) {
                corners[i] = get_sibling(j)->get_sibling(k)->get_sibling(l);
            }
            if (get_sibling(j)->get_sibling(l) != NULL && corners[i] == NULL) {
                corners[i] = get_sibling(j)->get_sibling(l)->get_sibling(k);
            }
        }
        if (get_sibling(k) != NULL && corners[i] == NULL) {
            if (get_sibling(k)->get_sibling(j) != NULL && corners[i] == NULL) {
                corners[i] = get_sibling(k)->get_sibling(j)->get_sibling(l);
            }
            if (get_sibling(k)->get_sibling(l) != NULL && corners[i] == NULL) {
                corners[i] = get_sibling(k)->get_sibling(l)->get_sibling(j);
            }
        }
        if (get_sibling(l) != NULL && corners[i] == NULL) {
            if (get_sibling(l)->get_sibling(j) != NULL && corners[i] == NULL) {
                corners[i] = get_sibling(l)->get_sibling(j)->get_sibling(k);
            }
            if (get_sibling(l)->get_sibling(k) != NULL && corners[i] == NULL) {
                corners[i] = get_sibling(l)->get_sibling(k)->get_sibling(j);
            }
        }
    }
//#pragma omp parallel for
    for (int i = 0; i < 6; i++) {
        if (get_sibling(i) != NULL) {
            neighbors[i] = dynamic_cast<FMM*>(get_sibling(i));
        } else {
            neighbors[i] = NULL;
        }
    }
//#pragma omp parallel for
    for (int i = 0; i < 12; i++) {
        if (edges[i] != NULL) {
            neighbors[i + 6] = dynamic_cast<FMM*>(edges[i]);
        } else {
            neighbors[i + 6] = NULL;
        }
    }
//#pragma omp parallel for
    for (int i = 0; i < 8; i++) {
        if (corners[i] != NULL) {
            neighbors[i + 18] = dynamic_cast<FMM*>(corners[i]);
        } else {
            neighbors[i + 18] = NULL;
        }
    }
}

Real FMM::g_lz(int i, int j, int k) const {
    return L(i + FBW - BW, j + FBW - BW, k + FBW - BW).g_lz();
}

void FMM::copy_pot_to_hydro_grid(int) {
    FMM* p;
    Hydro* g;
    Real pot;
    for (int k = BW - 1; k < GNX - BW + 1; k++) {
        for (int j = BW - 1; j < GNX - BW + 1; j++) {
            for (int i = BW - 1; i < GNX - BW + 1; i++) {
                pot = get_phi(i, j, k);
                pot *= U(i, j, k).rho();
                pot += U(i, j, k).rot_pot(X(i, j, k));
                U(i, j, k).set_pot(pot);
            }
        }
    }
    inc_instruction_pointer();
}

void FMM::pot_to_hydro_grid() {
    FMM* p;
    Hydro* g;
    Real pot;
    for (int n = 0; n < get_local_node_cnt(); n++) {
        p = dynamic_cast<FMM*>(get_local_node(n));
        g = dynamic_cast<Hydro*>(get_local_node(n));
//#pragma omp parallel for
        for (int k = BW - 1; k < GNX - BW + 1; k++) {
            for (int j = BW - 1; j < GNX - BW + 1; j++) {
                for (int i = BW - 1; i < GNX - BW + 1; i++) {
                    pot = p->get_phi(i, j, k);
                    pot *= (*g)(i, j, k).rho();
                    pot += (*g)(i, j, k).rot_pot(g->X(i, j, k));
                    (*g)(i, j, k).set_pot(pot);
                }
            }
        }
    }
}

void FMM::MPI_datatypes_init() {
    static bool initialized = false;
    if (!initialized) {
        for (int k = 0; k < INX; k++) {
            for (int j = 0; j < INX; j++) {
                for (int i = 0; i < INX; i++) {
                    if (!(i == 0 && j == 0 && k == 0)) {
                        d0_array[k][j][i] = -1.0 / sqrt(i * i + j * j + k * k);
                    }
                }
            }
        }
        for (int k = -INX; k < INX; k++) {
            for (int j = -INX; j < INX; j++) {
                for (int i = -INX; i < INX; i++) {
                    if (!(i == 0 && j == 0 && k == 0)) {
                        d1_array[k + INX][j + INX][i + INX][0] = -i / pow(i * i + j * j + k * k, 1.5);
                        d1_array[k + INX][j + INX][i + INX][1] = -j / pow(i * i + j * j + k * k, 1.5);
                        d1_array[k + INX][j + INX][i + INX][2] = -k / pow(i * i + j * j + k * k, 1.5);
                    }
                }
            }
        }
        const int lbi0 = FBW;
        const int ubi0 = 2 * FBW - 1;
        const int lbi1 = FNX - 2 * FBW;
        const int ubi1 = FNX - FBW - 1;
        const int lbe0 = 0;
        const int ube0 = FBW - 1;
        const int lbe1 = FNX - FBW;
        const int ube1 = FNX - 1;
        const int hbn0 = FNX - 1 - FBW;
        const int lbn0 = FBW;
        const int hbn1 = FNX - 1 - FBW;
        const int lbn1 = FBW;

        MPI_Type_contiguous(sizeof(multipole_t), MPI_BYTE, &MPI_multipole_t);
        MPI_Type_commit(&MPI_multipole_t);

        MPI_Type_contiguous(sizeof(multipole_dot_t), MPI_BYTE, &MPI_multipole_dot_t);
        MPI_Type_commit(&MPI_multipole_dot_t);

        MPI_Type_contiguous(sizeof(_4force_t), MPI_BYTE, &MPI_4force_t);
        MPI_Type_commit(&MPI_4force_t);
        MPI_Type_contiguous(sizeof(taylor_t), MPI_BYTE, &MPI_taylor_t);
        MPI_Type_commit(&MPI_taylor_t);

        MPI_Type_contiguous(sizeof(taylor_dot_t), MPI_BYTE, &MPI_taylor_dot_t);
        MPI_Type_commit(&MPI_taylor_dot_t);

        Array3d<multipole_dot_t, FNX, FNX, FNX>::mpi_datatype(MPI_send_bnd_dot_t + XL, lbi0, ubi0, lbn0, hbn0, lbn0, hbn0, MPI_multipole_dot_t);
        Array3d<multipole_dot_t, FNX, FNX, FNX>::mpi_datatype(MPI_send_bnd_dot_t + YL, lbn1, hbn1, lbi0, ubi0, lbn0, hbn0, MPI_multipole_dot_t);
        Array3d<multipole_dot_t, FNX, FNX, FNX>::mpi_datatype(MPI_send_bnd_dot_t + ZL, lbn1, hbn1, lbn1, hbn1, lbi0, ubi0, MPI_multipole_dot_t);

        Array3d<multipole_dot_t, FNX, FNX, FNX>::mpi_datatype(MPI_send_bnd_dot_t + XU, lbi1, ubi1, lbn0, hbn0, lbn0, hbn0, MPI_multipole_dot_t);
        Array3d<multipole_dot_t, FNX, FNX, FNX>::mpi_datatype(MPI_send_bnd_dot_t + YU, lbn1, hbn1, lbi1, ubi1, lbn0, hbn0, MPI_multipole_dot_t);
        Array3d<multipole_dot_t, FNX, FNX, FNX>::mpi_datatype(MPI_send_bnd_dot_t + ZU, lbn1, hbn1, lbn1, hbn1, lbi1, ubi1, MPI_multipole_dot_t);

        Array3d<multipole_dot_t, FNX, FNX, FNX>::mpi_datatype(MPI_recv_bnd_dot_t + XL, lbe0, ube0, lbn0, hbn0, lbn0, hbn0, MPI_multipole_dot_t);
        Array3d<multipole_dot_t, FNX, FNX, FNX>::mpi_datatype(MPI_recv_bnd_dot_t + YL, lbn1, hbn1, lbe0, ube0, lbn0, hbn0, MPI_multipole_dot_t);
        Array3d<multipole_dot_t, FNX, FNX, FNX>::mpi_datatype(MPI_recv_bnd_dot_t + ZL, lbn1, hbn1, lbn1, hbn1, lbe0, ube0, MPI_multipole_dot_t);

        Array3d<multipole_dot_t, FNX, FNX, FNX>::mpi_datatype(MPI_recv_bnd_dot_t + XU, lbe1, ube1, lbn0, hbn0, lbn0, hbn0, MPI_multipole_dot_t);
        Array3d<multipole_dot_t, FNX, FNX, FNX>::mpi_datatype(MPI_recv_bnd_dot_t + YU, lbn1, hbn1, lbe1, ube1, lbn0, hbn0, MPI_multipole_dot_t);
        Array3d<multipole_dot_t, FNX, FNX, FNX>::mpi_datatype(MPI_recv_bnd_dot_t + ZU, lbn1, hbn1, lbn1, hbn1, lbe1, ube1, MPI_multipole_dot_t);

        Array3d<multipole_dot_t, FNX, FNX, FNX>::mpi_datatype(MPI_send_bnd_dot_t + 6 + XLYL, lbi0, ubi0, lbi0, ubi0, lbn0, hbn0, MPI_multipole_dot_t);
        Array3d<multipole_dot_t, FNX, FNX, FNX>::mpi_datatype(MPI_send_bnd_dot_t + 6 + XLYU, lbi0, ubi0, lbi1, ubi1, lbn0, hbn0, MPI_multipole_dot_t);
        Array3d<multipole_dot_t, FNX, FNX, FNX>::mpi_datatype(MPI_send_bnd_dot_t + 6 + XUYL, lbi1, ubi1, lbi0, ubi0, lbn0, hbn0, MPI_multipole_dot_t);
        Array3d<multipole_dot_t, FNX, FNX, FNX>::mpi_datatype(MPI_send_bnd_dot_t + 6 + XUYU, lbi1, ubi1, lbi1, ubi1, lbn0, hbn0, MPI_multipole_dot_t);
        Array3d<multipole_dot_t, FNX, FNX, FNX>::mpi_datatype(MPI_send_bnd_dot_t + 6 + XLZL, lbi0, ubi0, lbn0, hbn0, lbi0, ubi0, MPI_multipole_dot_t);
        Array3d<multipole_dot_t, FNX, FNX, FNX>::mpi_datatype(MPI_send_bnd_dot_t + 6 + XLZU, lbi0, ubi0, lbn0, hbn0, lbi1, ubi1, MPI_multipole_dot_t);
        Array3d<multipole_dot_t, FNX, FNX, FNX>::mpi_datatype(MPI_send_bnd_dot_t + 6 + XUZL, lbi1, ubi1, lbn0, hbn0, lbi0, ubi0, MPI_multipole_dot_t);
        Array3d<multipole_dot_t, FNX, FNX, FNX>::mpi_datatype(MPI_send_bnd_dot_t + 6 + XUZU, lbi1, ubi1, lbn0, hbn0, lbi1, ubi1, MPI_multipole_dot_t);
        Array3d<multipole_dot_t, FNX, FNX, FNX>::mpi_datatype(MPI_send_bnd_dot_t + 6 + YLZL, lbn0, hbn0, lbi0, ubi0, lbi0, ubi0, MPI_multipole_dot_t);
        Array3d<multipole_dot_t, FNX, FNX, FNX>::mpi_datatype(MPI_send_bnd_dot_t + 6 + YLZU, lbn0, hbn0, lbi0, ubi0, lbi1, ubi1, MPI_multipole_dot_t);
        Array3d<multipole_dot_t, FNX, FNX, FNX>::mpi_datatype(MPI_send_bnd_dot_t + 6 + YUZL, lbn0, hbn0, lbi1, ubi1, lbi0, ubi0, MPI_multipole_dot_t);
        Array3d<multipole_dot_t, FNX, FNX, FNX>::mpi_datatype(MPI_send_bnd_dot_t + 6 + YUZU, lbn0, hbn0, lbi1, ubi1, lbi1, ubi1, MPI_multipole_dot_t);

        Array3d<multipole_dot_t, FNX, FNX, FNX>::mpi_datatype(MPI_recv_bnd_dot_t + 6 + XLYL, lbe0, ube0, lbe0, ube0, lbn0, hbn0, MPI_multipole_dot_t);
        Array3d<multipole_dot_t, FNX, FNX, FNX>::mpi_datatype(MPI_recv_bnd_dot_t + 6 + XLYU, lbe0, ube0, lbe1, ube1, lbn0, hbn0, MPI_multipole_dot_t);
        Array3d<multipole_dot_t, FNX, FNX, FNX>::mpi_datatype(MPI_recv_bnd_dot_t + 6 + XUYL, lbe1, ube1, lbe0, ube0, lbn0, hbn0, MPI_multipole_dot_t);
        Array3d<multipole_dot_t, FNX, FNX, FNX>::mpi_datatype(MPI_recv_bnd_dot_t + 6 + XUYU, lbe1, ube1, lbe1, ube1, lbn0, hbn0, MPI_multipole_dot_t);
        Array3d<multipole_dot_t, FNX, FNX, FNX>::mpi_datatype(MPI_recv_bnd_dot_t + 6 + XLZL, lbe0, ube0, lbn0, hbn0, lbe0, ube0, MPI_multipole_dot_t);
        Array3d<multipole_dot_t, FNX, FNX, FNX>::mpi_datatype(MPI_recv_bnd_dot_t + 6 + XLZU, lbe0, ube0, lbn0, hbn0, lbe1, ube1, MPI_multipole_dot_t);
        Array3d<multipole_dot_t, FNX, FNX, FNX>::mpi_datatype(MPI_recv_bnd_dot_t + 6 + XUZL, lbe1, ube1, lbn0, hbn0, lbe0, ube0, MPI_multipole_dot_t);
        Array3d<multipole_dot_t, FNX, FNX, FNX>::mpi_datatype(MPI_recv_bnd_dot_t + 6 + XUZU, lbe1, ube1, lbn0, hbn0, lbe1, ube1, MPI_multipole_dot_t);
        Array3d<multipole_dot_t, FNX, FNX, FNX>::mpi_datatype(MPI_recv_bnd_dot_t + 6 + YLZL, lbn0, hbn0, lbe0, ube0, lbe0, ube0, MPI_multipole_dot_t);
        Array3d<multipole_dot_t, FNX, FNX, FNX>::mpi_datatype(MPI_recv_bnd_dot_t + 6 + YLZU, lbn0, hbn0, lbe0, ube0, lbe1, ube1, MPI_multipole_dot_t);
        Array3d<multipole_dot_t, FNX, FNX, FNX>::mpi_datatype(MPI_recv_bnd_dot_t + 6 + YUZL, lbn0, hbn0, lbe1, ube1, lbe0, ube0, MPI_multipole_dot_t);
        Array3d<multipole_dot_t, FNX, FNX, FNX>::mpi_datatype(MPI_recv_bnd_dot_t + 6 + YUZU, lbn0, hbn0, lbe1, ube1, lbe1, ube1, MPI_multipole_dot_t);

        Array3d<multipole_dot_t, FNX, FNX, FNX>::mpi_datatype(MPI_send_bnd_dot_t + 18 + XLYLZL, lbi0, ubi0, lbi0, ubi0, lbi0, ubi0, MPI_multipole_dot_t);
        Array3d<multipole_dot_t, FNX, FNX, FNX>::mpi_datatype(MPI_send_bnd_dot_t + 18 + XLYLZU, lbi0, ubi0, lbi0, ubi0, lbi1, ubi1, MPI_multipole_dot_t);
        Array3d<multipole_dot_t, FNX, FNX, FNX>::mpi_datatype(MPI_send_bnd_dot_t + 18 + XLYUZL, lbi0, ubi0, lbi1, ubi1, lbi0, ubi0, MPI_multipole_dot_t);
        Array3d<multipole_dot_t, FNX, FNX, FNX>::mpi_datatype(MPI_send_bnd_dot_t + 18 + XLYUZU, lbi0, ubi0, lbi1, ubi1, lbi1, ubi1, MPI_multipole_dot_t);
        Array3d<multipole_dot_t, FNX, FNX, FNX>::mpi_datatype(MPI_send_bnd_dot_t + 18 + XUYLZL, lbi1, ubi1, lbi0, ubi0, lbi0, ubi0, MPI_multipole_dot_t);
        Array3d<multipole_dot_t, FNX, FNX, FNX>::mpi_datatype(MPI_send_bnd_dot_t + 18 + XUYLZU, lbi1, ubi1, lbi0, ubi0, lbi1, ubi1, MPI_multipole_dot_t);
        Array3d<multipole_dot_t, FNX, FNX, FNX>::mpi_datatype(MPI_send_bnd_dot_t + 18 + XUYUZL, lbi1, ubi1, lbi1, ubi1, lbi0, ubi0, MPI_multipole_dot_t);
        Array3d<multipole_dot_t, FNX, FNX, FNX>::mpi_datatype(MPI_send_bnd_dot_t + 18 + XUYUZU, lbi1, ubi1, lbi1, ubi1, lbi1, ubi1, MPI_multipole_dot_t);

        Array3d<multipole_dot_t, FNX, FNX, FNX>::mpi_datatype(MPI_recv_bnd_dot_t + 18 + XLYLZL, lbe0, ube0, lbe0, ube0, lbe0, ube0, MPI_multipole_dot_t);
        Array3d<multipole_dot_t, FNX, FNX, FNX>::mpi_datatype(MPI_recv_bnd_dot_t + 18 + XLYLZU, lbe0, ube0, lbe0, ube0, lbe1, ube1, MPI_multipole_dot_t);
        Array3d<multipole_dot_t, FNX, FNX, FNX>::mpi_datatype(MPI_recv_bnd_dot_t + 18 + XLYUZL, lbe0, ube0, lbe1, ube1, lbe0, ube0, MPI_multipole_dot_t);
        Array3d<multipole_dot_t, FNX, FNX, FNX>::mpi_datatype(MPI_recv_bnd_dot_t + 18 + XLYUZU, lbe0, ube0, lbe1, ube1, lbe1, ube1, MPI_multipole_dot_t);
        Array3d<multipole_dot_t, FNX, FNX, FNX>::mpi_datatype(MPI_recv_bnd_dot_t + 18 + XUYLZL, lbe1, ube1, lbe0, ube0, lbe0, ube0, MPI_multipole_dot_t);
        Array3d<multipole_dot_t, FNX, FNX, FNX>::mpi_datatype(MPI_recv_bnd_dot_t + 18 + XUYLZU, lbe1, ube1, lbe0, ube0, lbe1, ube1, MPI_multipole_dot_t);
        Array3d<multipole_dot_t, FNX, FNX, FNX>::mpi_datatype(MPI_recv_bnd_dot_t + 18 + XUYUZL, lbe1, ube1, lbe1, ube1, lbe0, ube0, MPI_multipole_dot_t);
        Array3d<multipole_dot_t, FNX, FNX, FNX>::mpi_datatype(MPI_recv_bnd_dot_t + 18 + XUYUZU, lbe1, ube1, lbe1, ube1, lbe1, ube1, MPI_multipole_dot_t);

        Array3d<multipole_t, FNX, FNX, FNX>::mpi_datatype(MPI_send_bnd_t + XL, lbi0, ubi0, lbn0, hbn0, lbn0, hbn0, MPI_multipole_t);
        Array3d<multipole_t, FNX, FNX, FNX>::mpi_datatype(MPI_send_bnd_t + YL, lbn1, hbn1, lbi0, ubi0, lbn0, hbn0, MPI_multipole_t);
        Array3d<multipole_t, FNX, FNX, FNX>::mpi_datatype(MPI_send_bnd_t + ZL, lbn1, hbn1, lbn1, hbn1, lbi0, ubi0, MPI_multipole_t);

        Array3d<multipole_t, FNX, FNX, FNX>::mpi_datatype(MPI_send_bnd_t + XU, lbi1, ubi1, lbn0, hbn0, lbn0, hbn0, MPI_multipole_t);
        Array3d<multipole_t, FNX, FNX, FNX>::mpi_datatype(MPI_send_bnd_t + YU, lbn1, hbn1, lbi1, ubi1, lbn0, hbn0, MPI_multipole_t);
        Array3d<multipole_t, FNX, FNX, FNX>::mpi_datatype(MPI_send_bnd_t + ZU, lbn1, hbn1, lbn1, hbn1, lbi1, ubi1, MPI_multipole_t);

        Array3d<multipole_t, FNX, FNX, FNX>::mpi_datatype(MPI_recv_bnd_t + XL, lbe0, ube0, lbn0, hbn0, lbn0, hbn0, MPI_multipole_t);
        Array3d<multipole_t, FNX, FNX, FNX>::mpi_datatype(MPI_recv_bnd_t + YL, lbn1, hbn1, lbe0, ube0, lbn0, hbn0, MPI_multipole_t);
        Array3d<multipole_t, FNX, FNX, FNX>::mpi_datatype(MPI_recv_bnd_t + ZL, lbn1, hbn1, lbn1, hbn1, lbe0, ube0, MPI_multipole_t);

        Array3d<multipole_t, FNX, FNX, FNX>::mpi_datatype(MPI_recv_bnd_t + XU, lbe1, ube1, lbn0, hbn0, lbn0, hbn0, MPI_multipole_t);
        Array3d<multipole_t, FNX, FNX, FNX>::mpi_datatype(MPI_recv_bnd_t + YU, lbn1, hbn1, lbe1, ube1, lbn0, hbn0, MPI_multipole_t);
        Array3d<multipole_t, FNX, FNX, FNX>::mpi_datatype(MPI_recv_bnd_t + ZU, lbn1, hbn1, lbn1, hbn1, lbe1, ube1, MPI_multipole_t);

        Array3d<multipole_t, FNX, FNX, FNX>::mpi_datatype(MPI_send_bnd_t + 6 + XLYL, lbi0, ubi0, lbi0, ubi0, lbn0, hbn0, MPI_multipole_t);
        Array3d<multipole_t, FNX, FNX, FNX>::mpi_datatype(MPI_send_bnd_t + 6 + XLYU, lbi0, ubi0, lbi1, ubi1, lbn0, hbn0, MPI_multipole_t);
        Array3d<multipole_t, FNX, FNX, FNX>::mpi_datatype(MPI_send_bnd_t + 6 + XUYL, lbi1, ubi1, lbi0, ubi0, lbn0, hbn0, MPI_multipole_t);
        Array3d<multipole_t, FNX, FNX, FNX>::mpi_datatype(MPI_send_bnd_t + 6 + XUYU, lbi1, ubi1, lbi1, ubi1, lbn0, hbn0, MPI_multipole_t);
        Array3d<multipole_t, FNX, FNX, FNX>::mpi_datatype(MPI_send_bnd_t + 6 + XLZL, lbi0, ubi0, lbn0, hbn0, lbi0, ubi0, MPI_multipole_t);
        Array3d<multipole_t, FNX, FNX, FNX>::mpi_datatype(MPI_send_bnd_t + 6 + XLZU, lbi0, ubi0, lbn0, hbn0, lbi1, ubi1, MPI_multipole_t);
        Array3d<multipole_t, FNX, FNX, FNX>::mpi_datatype(MPI_send_bnd_t + 6 + XUZL, lbi1, ubi1, lbn0, hbn0, lbi0, ubi0, MPI_multipole_t);
        Array3d<multipole_t, FNX, FNX, FNX>::mpi_datatype(MPI_send_bnd_t + 6 + XUZU, lbi1, ubi1, lbn0, hbn0, lbi1, ubi1, MPI_multipole_t);
        Array3d<multipole_t, FNX, FNX, FNX>::mpi_datatype(MPI_send_bnd_t + 6 + YLZL, lbn0, hbn0, lbi0, ubi0, lbi0, ubi0, MPI_multipole_t);
        Array3d<multipole_t, FNX, FNX, FNX>::mpi_datatype(MPI_send_bnd_t + 6 + YLZU, lbn0, hbn0, lbi0, ubi0, lbi1, ubi1, MPI_multipole_t);
        Array3d<multipole_t, FNX, FNX, FNX>::mpi_datatype(MPI_send_bnd_t + 6 + YUZL, lbn0, hbn0, lbi1, ubi1, lbi0, ubi0, MPI_multipole_t);
        Array3d<multipole_t, FNX, FNX, FNX>::mpi_datatype(MPI_send_bnd_t + 6 + YUZU, lbn0, hbn0, lbi1, ubi1, lbi1, ubi1, MPI_multipole_t);

        Array3d<multipole_t, FNX, FNX, FNX>::mpi_datatype(MPI_recv_bnd_t + 6 + XLYL, lbe0, ube0, lbe0, ube0, lbn0, hbn0, MPI_multipole_t);
        Array3d<multipole_t, FNX, FNX, FNX>::mpi_datatype(MPI_recv_bnd_t + 6 + XLYU, lbe0, ube0, lbe1, ube1, lbn0, hbn0, MPI_multipole_t);
        Array3d<multipole_t, FNX, FNX, FNX>::mpi_datatype(MPI_recv_bnd_t + 6 + XUYL, lbe1, ube1, lbe0, ube0, lbn0, hbn0, MPI_multipole_t);
        Array3d<multipole_t, FNX, FNX, FNX>::mpi_datatype(MPI_recv_bnd_t + 6 + XUYU, lbe1, ube1, lbe1, ube1, lbn0, hbn0, MPI_multipole_t);
        Array3d<multipole_t, FNX, FNX, FNX>::mpi_datatype(MPI_recv_bnd_t + 6 + XLZL, lbe0, ube0, lbn0, hbn0, lbe0, ube0, MPI_multipole_t);
        Array3d<multipole_t, FNX, FNX, FNX>::mpi_datatype(MPI_recv_bnd_t + 6 + XLZU, lbe0, ube0, lbn0, hbn0, lbe1, ube1, MPI_multipole_t);
        Array3d<multipole_t, FNX, FNX, FNX>::mpi_datatype(MPI_recv_bnd_t + 6 + XUZL, lbe1, ube1, lbn0, hbn0, lbe0, ube0, MPI_multipole_t);
        Array3d<multipole_t, FNX, FNX, FNX>::mpi_datatype(MPI_recv_bnd_t + 6 + XUZU, lbe1, ube1, lbn0, hbn0, lbe1, ube1, MPI_multipole_t);
        Array3d<multipole_t, FNX, FNX, FNX>::mpi_datatype(MPI_recv_bnd_t + 6 + YLZL, lbn0, hbn0, lbe0, ube0, lbe0, ube0, MPI_multipole_t);
        Array3d<multipole_t, FNX, FNX, FNX>::mpi_datatype(MPI_recv_bnd_t + 6 + YLZU, lbn0, hbn0, lbe0, ube0, lbe1, ube1, MPI_multipole_t);
        Array3d<multipole_t, FNX, FNX, FNX>::mpi_datatype(MPI_recv_bnd_t + 6 + YUZL, lbn0, hbn0, lbe1, ube1, lbe0, ube0, MPI_multipole_t);
        Array3d<multipole_t, FNX, FNX, FNX>::mpi_datatype(MPI_recv_bnd_t + 6 + YUZU, lbn0, hbn0, lbe1, ube1, lbe1, ube1, MPI_multipole_t);

        Array3d<multipole_t, FNX, FNX, FNX>::mpi_datatype(MPI_send_bnd_t + 18 + XLYLZL, lbi0, ubi0, lbi0, ubi0, lbi0, ubi0, MPI_multipole_t);
        Array3d<multipole_t, FNX, FNX, FNX>::mpi_datatype(MPI_send_bnd_t + 18 + XLYLZU, lbi0, ubi0, lbi0, ubi0, lbi1, ubi1, MPI_multipole_t);
        Array3d<multipole_t, FNX, FNX, FNX>::mpi_datatype(MPI_send_bnd_t + 18 + XLYUZL, lbi0, ubi0, lbi1, ubi1, lbi0, ubi0, MPI_multipole_t);
        Array3d<multipole_t, FNX, FNX, FNX>::mpi_datatype(MPI_send_bnd_t + 18 + XLYUZU, lbi0, ubi0, lbi1, ubi1, lbi1, ubi1, MPI_multipole_t);
        Array3d<multipole_t, FNX, FNX, FNX>::mpi_datatype(MPI_send_bnd_t + 18 + XUYLZL, lbi1, ubi1, lbi0, ubi0, lbi0, ubi0, MPI_multipole_t);
        Array3d<multipole_t, FNX, FNX, FNX>::mpi_datatype(MPI_send_bnd_t + 18 + XUYLZU, lbi1, ubi1, lbi0, ubi0, lbi1, ubi1, MPI_multipole_t);
        Array3d<multipole_t, FNX, FNX, FNX>::mpi_datatype(MPI_send_bnd_t + 18 + XUYUZL, lbi1, ubi1, lbi1, ubi1, lbi0, ubi0, MPI_multipole_t);
        Array3d<multipole_t, FNX, FNX, FNX>::mpi_datatype(MPI_send_bnd_t + 18 + XUYUZU, lbi1, ubi1, lbi1, ubi1, lbi1, ubi1, MPI_multipole_t);

        Array3d<multipole_t, FNX, FNX, FNX>::mpi_datatype(MPI_recv_bnd_t + 18 + XLYLZL, lbe0, ube0, lbe0, ube0, lbe0, ube0, MPI_multipole_t);
        Array3d<multipole_t, FNX, FNX, FNX>::mpi_datatype(MPI_recv_bnd_t + 18 + XLYLZU, lbe0, ube0, lbe0, ube0, lbe1, ube1, MPI_multipole_t);
        Array3d<multipole_t, FNX, FNX, FNX>::mpi_datatype(MPI_recv_bnd_t + 18 + XLYUZL, lbe0, ube0, lbe1, ube1, lbe0, ube0, MPI_multipole_t);
        Array3d<multipole_t, FNX, FNX, FNX>::mpi_datatype(MPI_recv_bnd_t + 18 + XLYUZU, lbe0, ube0, lbe1, ube1, lbe1, ube1, MPI_multipole_t);
        Array3d<multipole_t, FNX, FNX, FNX>::mpi_datatype(MPI_recv_bnd_t + 18 + XUYLZL, lbe1, ube1, lbe0, ube0, lbe0, ube0, MPI_multipole_t);
        Array3d<multipole_t, FNX, FNX, FNX>::mpi_datatype(MPI_recv_bnd_t + 18 + XUYLZU, lbe1, ube1, lbe0, ube0, lbe1, ube1, MPI_multipole_t);
        Array3d<multipole_t, FNX, FNX, FNX>::mpi_datatype(MPI_recv_bnd_t + 18 + XUYUZL, lbe1, ube1, lbe1, ube1, lbe0, ube0, MPI_multipole_t);
        Array3d<multipole_t, FNX, FNX, FNX>::mpi_datatype(MPI_recv_bnd_t + 18 + XUYUZU, lbe1, ube1, lbe1, ube1, lbe1, ube1, MPI_multipole_t);

        const int ds = INX / 2;
        int ci, xlb, xub, ylb, yub, zlb, zub;

        for (int k = 0; k < 2; k++) {
            zlb = FBW + k * (INX / 2);
            zub = zlb + ds - 1;
            for (int j = 0; j < 2; j++) {
                ylb = FBW + j * (INX / 2);
                yub = ylb + ds - 1;
                for (int i = 0; i < 2; i++) {
                    xlb = FBW + i * (INX / 2);
                    xub = xlb + ds - 1;
                    ci = 4 * k + 2 * j + i;
                    Array3d<multipole_dot_t, FNX, FNX, FNX>::mpi_datatype(MPI_comm_child_poles_dot_t + ci, xlb, xub, ylb, yub, zlb, zub, MPI_multipole_dot_t);
                    Array3d<taylor_dot_t, FNX, FNX, FNX>::mpi_datatype(MPI_comm_taylor_dot_t + ci, xlb, xub, ylb, yub, zlb, zub, MPI_taylor_dot_t);
                    Array3d<multipole_t, FNX, FNX, FNX>::mpi_datatype(MPI_comm_child_poles_t + ci, xlb, xub, ylb, yub, zlb, zub, MPI_multipole_t);
                    Array3d<taylor_t, FNX, FNX, FNX>::mpi_datatype(MPI_comm_taylor_t + ci, xlb, xub, ylb, yub, zlb, zub, MPI_taylor_t);
                    Array3d<_4force_t, FNX, FNX, FNX>::mpi_datatype(MPI_comm_child3_t + ci, xlb, xub, ylb, yub, zlb, zub, MPI_4force_t);
                }
            }
        }
        face_id[XL] = 0;
        face_id[XU] = 1;
        face_id[YL] = 2;
        face_id[YU] = 3;
        face_id[ZL] = 4;
        face_id[ZU] = 5;
        face_id[XLYL + 6] = 6;
        face_id[XLYU + 6] = 7;
        face_id[XLZL + 6] = 8;
        face_id[XLZU + 6] = 9;
        face_id[XUYL + 6] = 10;
        face_id[XUYU + 6] = 11;
        face_id[XUZL + 6] = 12;
        face_id[XUZU + 6] = 13;
        face_id[YLZL + 6] = 14;
        face_id[YLZU + 6] = 15;
        face_id[YUZL + 6] = 16;
        face_id[YUZU + 6] = 17;
        face_id[XLYLZL + 18] = 18;
        face_id[XLYLZU + 18] = 19;
        face_id[XLYUZL + 18] = 20;
        face_id[XLYUZU + 18] = 21;
        face_id[XUYLZL + 18] = 22;
        face_id[XUYLZU + 18] = 23;
        face_id[XUYUZL + 18] = 24;
        face_id[XUYUZU + 18] = 25;
        face_opp_id[XU] = 0;
        face_opp_id[XL] = 1;
        face_opp_id[YU] = 2;
        face_opp_id[YL] = 3;
        face_opp_id[ZU] = 4;
        face_opp_id[ZL] = 5;
        face_opp_id[XUYU + 6] = 6;
        face_opp_id[XUYL + 6] = 7;
        face_opp_id[XUZU + 6] = 8;
        face_opp_id[XUZL + 6] = 9;
        face_opp_id[XLYU + 6] = 10;
        face_opp_id[XLYL + 6] = 11;
        face_opp_id[XLZU + 6] = 12;
        face_opp_id[XLZL + 6] = 13;
        face_opp_id[YUZU + 6] = 14;
        face_opp_id[YUZL + 6] = 15;
        face_opp_id[YLZU + 6] = 16;
        face_opp_id[YLZL + 6] = 17;
        face_opp_id[XUYUZU + 18] = 18;
        face_opp_id[XUYUZL + 18] = 19;
        face_opp_id[XUYLZU + 18] = 20;
        face_opp_id[XUYLZL + 18] = 21;
        face_opp_id[XLYUZU + 18] = 22;
        face_opp_id[XLYUZL + 18] = 23;
        face_opp_id[XLYLZU + 18] = 24;
        face_opp_id[XLYLZL + 18] = 25;
        initialized = true;
    }
}

#endif
