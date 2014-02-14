#include <mpi.h>
#include <typeinfo>
#include "hydro_grid.h"
#include "../bits.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "../tag.h"
#include "../indexer3d.h"
#include "../reconstruct.h"

#ifdef USE_HYDRO_GRID

void HydroGrid::write(FILE* fp) const {
    fwrite(&dx, sizeof(Real), 1, fp);
    if (get_level() == 0) {
        fwrite(&h0, sizeof(Real), 1, fp);
        fwrite(&FO, sizeof(State), 1, fp);
    }
    if (proc() == MPI_rank()) {
        fwrite(U.ptr(), sizeof(State), GNX * GNX * GNX, fp);
    }
}

void HydroGrid::read(FILE* fp) {
    fread(&dx, sizeof(Real), 1, fp);
    if (get_level() == 0) {
        fread(&h0, sizeof(Real), 1, fp);
        fread(&FO, sizeof(State), 1, fp);
    }
    if (proc() == MPI_rank()) {
        fread(U.ptr(), sizeof(State), GNX * GNX * GNX, fp);
    }
}

void HydroGrid::inject_from_parent(ChildIndex c) {
    const HydroGrid* p = dynamic_cast<const HydroGrid*>(get_parent());
    Vector<Real, STATE_NF> s1, s2, s3;
    State u;
    int k, j, k0, j0, i, i0;
    for (k = BW; k < GNX - BW; k += 2) {
        for (j = BW; j < GNX - BW; j += 2) {
            k0 = (BW + k) / 2 + c.get_z() * (GNX / 2 - BW);
            j0 = (BW + j) / 2 + c.get_y() * (GNX / 2 - BW);
            for (i = BW, i0 = BW + c.get_x() * (GNX / 2 - BW); i < GNX - BW; i += 2, i0++) {u
                = (*p)(i0, j0, k0);
                U(i + 0, j + 0, k + 0) = u;
                U(i + 1, j + 0, k + 0) = u;
                U(i + 0, j + 1, k + 0) = u;
                U(i + 1, j + 1, k + 0) = u;
                U(i + 0, j + 0, k + 1) = u;
                U(i + 1, j + 0, k + 1) = u;
                U(i + 0, j + 1, k + 1) = u;
                U(i + 1, j + 1, k + 1) = u;
            }
        }
    }
}

Real HydroGrid::get_dx() const {
    return dx;
}

void HydroGrid::create_child(const ChildIndex& c) {
    Vector<int, 3> child_offset;
    HydroGrid* child;
    OctNode::create_child(c);
    if (MPI_rank() == get_child(c)->proc()) {
        dynamic_cast<HydroGrid*>(get_child(c))->inject_from_parent(c);
    }
    child = dynamic_cast<HydroGrid*>(get_child(c));
    child_offset = (offset * 2 + BW) + c.vector() * (GNX - 2 * BW);
    child->dx = h0 / Real(1 << (get_level() + 1));
    child->offset = child_offset;
}

Real HydroGrid::get_output_point(int i, int j, int k, int l) const {
    return (*this)(i, j, k)[l];
}

int HydroGrid::nvar_output() const {
    return STATE_NF;
}

const char* HydroGrid::output_field_names(int i) const {
    return State::field_name(i);
}

bool HydroGrid::zone_is_refined(int i, int j, int k) const {
    ChildIndex c;
    i = (2 * i) / GNX;
    j = (2 * j) / GNX;
    k = (2 * k) / GNX;
    c.set_index(i, j, k);
    return (get_child(c) != NULL);
}

HydroGrid::HydroGrid() :
        OctNode() {
    static bool init = false;
    if (!init) {
        origin[0] = 0.0;
        origin[1] = 0.0;
        origin[2] = 0.0;
        init = true;
    }
    mpi_datatypes_initialize();
    offset = 0;
    time = 0.0;

    for (int i = 0; i < OCT_NCHILD; i++) {
        send_request[i] = recv_request[i] = MPI_REQUEST_NULL;
        for (int j = 0; j < 3; j++) {
            amr_recv_request[2 * j][i] = MPI_REQUEST_NULL;
            amr_recv_request[2 * j + 1][i] = MPI_REQUEST_NULL;
            amr_send_request[j][i] = MPI_REQUEST_NULL;
        }
    }
}

void HydroGrid::deallocate_arrays() {
    for (int i = 0; i < 3; i++) {
        Fv[i].deallocate();
        Fs[i].deallocate();
    }
    U.deallocate();
    D.deallocate();
    U0.deallocate();
    if (shadow) {
        E0.deallocate();
        E.deallocate();
    }
}

void HydroGrid::allocate_arrays() {
    for (int i = 0; i < 3; i++) {
        Fv[i].allocate();
        Fs[i].allocate();
    }
    U.allocate();
    D.allocate();
    U0.allocate();
    if (shadow) {
        E.allocate();
        E0.allocate();
    }
}

void HydroGrid::init() {
    OctNode::init();
    dx = h0;
}

HydroGrid::~HydroGrid() {
}

Vector<State, 4> HydroGrid::state_sum() const {
    Vector<State, 4> s0, s1;
    Real h3 = pow(dx, 3);
    Real a0, a1, a2, a3;
    Real buffer[STATE_NF * 4];
    if (this->proc() == MPI_rank()) {
        for (int l = 0; l < STATE_NF; l++) {
            a0 = a1 = a2 = a3 = 0.0;
            Real ds;
            int k, j, i;
            for (k = BW; k < GNX - BW; k++) {
                for (j = BW; j < GNX - BW; j++) {
                    for (i = BW; i < GNX - BW; i++) {
                        //		if (!zone_is_refined(i, j, k)) {
                        ds = U(i, j, k)[l] * h3;
                        a0 += ds;
                        a1 += ds * HydroGrid::xc(i);
                        a2 += ds * HydroGrid::yc(j);
                        a3 += ds * HydroGrid::zc(k);
                        //	}
                    }
                }
            }
            s0[0][l] = a0;
            s0[1][l] = a1;
            s0[2][l] = a2;
            s0[3][l] = a3;
        }
        int cnt = 0;
        for (int i = 0; i < STATE_NF; i++) {
            for (int l = 0; l < 4; l++) {
                buffer[cnt] = s0[l][i];
                cnt++;
            }
        }
    }
    MPI_Bcast(buffer, 4 * STATE_NF, MPI_DOUBLE_PRECISION, this->proc(), MPI_COMM_WORLD );
    if (this->proc() != MPI_rank()) {
        int cnt = 0;
        for (int i = 0; i < STATE_NF; i++) {
            for (int l = 0; l < 4; l++) {
                s0[l][i] = buffer[cnt];
                cnt++;
            }
        }
    }

    return s0;
}

void HydroGrid::amr_bnd_send(int dir) {
    HydroGrid* g;
    amr_cnt[dir] = 0;
    State slp, u0, up, um;
    Vector<int, 3> ub, lb, o, v0, vp, vm, v1;
    int cnt, tag, i0, j0, k0;
    for (ChildIndex i = 0; i < 8; i++) {
        g = dynamic_cast<HydroGrid*>(get_child(i));
        if (g) {
            for (int f = 2 * dir; f < 2 * dir + 2; f++) {
                ub = GNX - BW - 1;
                lb = BW;
                if (f % 2 == 0) {
                    lb[f / 2] = 0;
                    ub[f / 2] = BW - 1;
                } else {
                    lb[f / 2] = GNX - BW;
                    ub[f / 2] = GNX - 1;
                }
                o = i.vec();
                o *= (GNX - 2 * BW);
                o += BW;
                if (g->is_amr_bound(f)) {
                    State tmp, tmp2;
                    const int sz = (GNX - 2 * BW) * (GNX - 2 * BW) * BW;
                    cnt = 0;
                    mpi_amr_buffer[2 * dir][amr_cnt[dir]] = new State[sz];
                    for (int k = lb[2]; k <= ub[2]; k++) {
                        k0 = (o[2] + k) / 2;
                        for (int j = lb[1]; j <= ub[1]; j++) {
                            j0 = (o[1] + j) / 2;
                            for (int i = lb[0]; i <= ub[0]; i++) {
                                i0 = (o[0] + i) / 2;
                                tmp = U(i0, j0, k0);
                              //  int di;
                               // if (dir == 0) {
                                  //  di = 2 * ((o[0] + i) % 2) - 1;
                                  //  tmp2 = U(i0 + di, j0, k0);
                                //    //    tmp2.to_prim(this->X(i0 + di, j0, k0));
                              //  } else if (dir == 1) {
                                   // di = 2 * ((o[1] + j) % 2) - 1;
                                 //   tmp2 = U(i0, j0 + di, k0);
                               //     //      tmp2.to_prim(this->X(i0, j0 + di, k0));
                             //   } else {
                           //         di = 2 * ((o[2] + k) % 2) - 1;
                         //           tmp2 = U(i0, j0, k0 + di);
                                    //        tmp2.to_prim(this->X(0, j0, k0 + di));
                               // }
                                //         tmp *= 0.75;
                                //        tmp2 *= 0.25;
                                //.       tmp += tmp2;
                                tmp.to_prim(this->X(i0, j0, k0));
                                tmp.from_prim(g->X(i, j, k));
                                mpi_amr_buffer[2 * dir][amr_cnt[dir]][cnt] = tmp;
                                cnt++;
                            }
                        }
                    }assert( cnt == sz);
                    tag = tag_gen(TAG_FLUX, g->get_id(), f);
                    MPI_Isend(mpi_amr_buffer[2 * dir][amr_cnt[dir]], cnt, MPI_state_t, g->proc(), tag, MPI_COMM_WORLD, &(amr_send_request[dir][amr_cnt[dir]]));
                    amr_cnt[dir]++;
                    if (amr_cnt[dir] > 8) {
                        printf("Internal amr buffer xceeded = %i\n", amr_cnt[dir]);
                        abort();
                    }
                }
            }
        }
    }
    inc_instruction_pointer(dir);
}

void HydroGrid::amr_bnd_send_wait(int dir) {
    int flag;
    MPI_Testall(amr_cnt[dir], amr_send_request[dir], &flag, MPI_STATUS_IGNORE );
    if (flag) {
        for (int i = 0; i < amr_cnt[dir]; i++) {
            delete[] mpi_amr_buffer[2 * dir][i];
        }
        inc_instruction_pointer(dir);
    }
}

void HydroGrid::inject_from_children_send(int dir) {
    if (dir == 0) {
        int cnt, tag;
        if (get_level() != 0) {
            mpi_buffer[0] = new State[(GNX / 2 - BW) * (GNX / 2 - BW) * (GNX / 2 - BW)];
            cnt = 0;
            for (int k = BW; k < GNX - BW; k += 2) {
                for (int j = BW; j < GNX - BW; j += 2) {
                    for (int i = BW; i < GNX - BW; i += 2) {
                        mpi_buffer[0][cnt] = U.oct_avg(i, j, k);
                        cnt++;
                    }
                }
            }
            tag = tag_gen(TAG_C2P, get_id(), 0);
            MPI_Isend(mpi_buffer[0], cnt, MPI_state_t, get_parent()->proc(), tag, MPI_COMM_WORLD, send_request);
        }
    }
    inc_instruction_pointer(dir);
}

void HydroGrid::inject_from_children_recv(int dir) {
    if (dir == 0) {
        int tag;
        OctNode* child;
        for (ChildIndex ci = 0; ci < 8; ci++) {
            child = get_child(ci);
            if (child != NULL) {
                tag = tag_gen(TAG_C2P, child->get_id(), 0);
                MPI_Irecv(U.ptr(), 1, MPI_child_t[ci], child->proc(), tag, MPI_COMM_WORLD, recv_request + ci);
            }
        }
    }
    inc_instruction_pointer(dir);
}

void HydroGrid::inject_from_children_send_wait(int dir) {
    int flag;
    if (dir == 0) {
        MPI_Test(send_request, &flag, MPI_STATUS_IGNORE );
        if (flag && get_level() != 0) {
            delete[] mpi_buffer[0];
        }
    } else {
        flag = true;
    }
    if (flag) {
        inc_instruction_pointer(dir);
    }
}

void HydroGrid::inject_from_children_recv_wait(int dir) {
    int flag;
    if (dir == 0) {
        MPI_Testall(8, recv_request, &flag, MPI_STATUS_IGNORE );
    } else {
        flag = true;
    }
    if (flag) {
        inc_instruction_pointer(dir);
    }

}

void HydroGrid::flux_bnd_comm(int dir) {
    int send_tag, recv_tag, f, gproc;
    for (f = 2 * dir; f < 2 * dir + 2; f++) {
        if (is_real_bound(f)) {
            send_tag = tag_gen(TAG_FLUX, get_sibling(f)->get_id(), f ^ 1);
            MPI_Isend(U.ptr(), 1, MPI_guard_send_t[f], get_sibling(f)->proc(), send_tag, MPI_COMM_WORLD, send_request + f);
        }
        if (!is_phys_bound(f)) {
            if (is_real_bound(f)) {
                gproc = get_sibling(f)->proc();
            } else {
                gproc = get_parent()->proc();
            }
            recv_tag = tag_gen(TAG_FLUX, get_id(), f);
            MPI_Irecv(U.ptr(), 1, MPI_guard_recv_t[f], gproc, recv_tag, MPI_COMM_WORLD, recv_request + f);
        }
    }
    inc_instruction_pointer(dir);
}

void HydroGrid::flux_bnd_recv_wait(int dir) {
    int flag;
    MPI_Testall(2, recv_request + 2 * dir, &flag, MPI_STATUS_IGNORE );
    if (flag) {
        inc_instruction_pointer(dir);
    }
}

void HydroGrid::flux_bnd_send_wait(int dir) {
    int flag;
    MPI_Testall(2, send_request + 2 * dir, &flag, MPI_STATUS_IGNORE );
    if (flag) {
        inc_instruction_pointer(dir);
    }
}

void HydroGrid::flux_compute(int dir) {
    State q0[GNX], ql[GNX], qr[GNX];
    int k, j, i;
    _3Vec x;
    Real a;
    for (k = BW; k < GNX - BW; k++) {
        for (j = BW; j < GNX - BW; j++) {
            if (dir == 0) {
                for (i = 0; i < GNX; i++) {
                    q0[i] = U(i, j, k);
                    q0[i].to_prim(HydroGrid::X(i, j, k));
                }
                reconstruct(q0, ql, qr);
                for (i = BW; i < GNX - BW + 1; i++) {
                    x = HydroGrid::Xfx(i, j, k);
                    ql[i].from_prim(x);
                    qr[i].from_prim(x);
                    a = max(ql[i].max_abs_x_eigen(x), qr[i].max_abs_x_eigen(x));
                    Fv[0](i, j, k) = ((ql[i].x_vector_flux(x) + qr[i].x_vector_flux(x)) - (qr[i] - ql[i]) * a) * 0.5;
                    Fs[0](i, j, k) = (ql[i].scalar_flux(x) + qr[i].scalar_flux(x)) * 0.5;

                }
            } else if (dir == 1) {
                for (i = 0; i < GNX; i++) {
                    q0[i] = U(j, i, k);
                    q0[i].to_prim(HydroGrid::X(j, i, k));
                }
                reconstruct(q0, ql, qr);
                for (i = BW; i < GNX - BW + 1; i++) {
                    x = HydroGrid::Xfy(j, i, k);
                    ql[i].from_prim(x);
                    qr[i].from_prim(x);
                    a = max(ql[i].max_abs_y_eigen(x), qr[i].max_abs_y_eigen(x));
                    Fv[1](j, i, k) = ((ql[i].y_vector_flux(x) + qr[i].y_vector_flux(x)) - (qr[i] - ql[i]) * a) * 0.5;
                    Fs[1](j, i, k) = (ql[i].scalar_flux(x) + qr[i].scalar_flux(x)) * 0.5;
                }
            } else {
                for (i = 0; i < GNX; i++) {
                    q0[i] = U(j, k, i);
                    q0[i].to_prim(HydroGrid::X(j, k, i));
                }
                reconstruct(q0, ql, qr);
                for (i = BW; i < GNX - BW + 1; i++) {
                    x = HydroGrid::Xfz(j, k, i);
                    ql[i].from_prim(x);
                    qr[i].from_prim(x);
                    a = max(ql[i].max_abs_z_eigen(x), qr[i].max_abs_z_eigen(x));
                    Fv[2](j, k, i) = ((ql[i].z_vector_flux(x) + qr[i].z_vector_flux(x)) - (qr[i] - ql[i]) * a) * 0.5;
                    Fs[2](j, k, i) = (ql[i].scalar_flux(x) + qr[i].scalar_flux(x)) * 0.5;
                }
            }
        }
    }
    flux_physical_bounds(dir);
    inc_instruction_pointer(dir);
}

void HydroGrid::flux_physical_bounds(int dir) {
    Vector<int, 3> lb, ub;
    if (is_phys_bound(2 * dir)) {
        lb = BW;
        ub = GNX - BW - 1;
        lb[dir] = BW;
        ub[dir] = BW;
        for (Indexer3d i(lb, ub); !i.end(); i++) {
            for (int l = 0; l < STATE_NF; l++) {
                if (U(i)[l] > 0.0) {
                    (Fv[dir](i))[l] = min(0.0, (Fv[dir](i))[l]);
                } else {
                    (Fv[dir](i))[l] = max(0.0, (Fv[dir](i))[l]);
                }
                            (Fv[dir](i))[l] = 0.0;
            }
        }
    }
    if (is_phys_bound(2 * dir + 1)) {
        lb = BW;
        ub = GNX - BW - 1;
        lb[dir] = GNX - BW;
        ub[dir] = GNX - BW;
        for (Indexer3d i(lb, ub); !i.end(); i++) {
            for (int l = 0; l < STATE_NF; l++) {
                if (U(i)[l] > 0.0) {
                    (Fv[dir](i))[l] = max(0.0, (Fv[dir](i))[l]);
                } else {
                    (Fv[dir](i))[l] = min(0.0, (Fv[dir](i))[l]);
                }
                        (Fv[dir](i))[l] = 0.0;
            }
        }
    }
}

void HydroGrid::flux_cf_adjust_send(int dir) {
    const int sz = (GNX - 2 * BW) * (GNX - 2 * BW) / 4;
    State v;
    HydroGrid* p;
    int cnt, i, tag;
    if (dir == 0) {
        for (int f = 0; f < 2; f++) {
            if (is_amr_bound(f)) {
                mpi_buffer[f] = new State[sz];
                i = BW + f * (GNX - 2 * BW);
                cnt = 0;
                for (int k = BW; k < GNX - BW; k += 2) {
                    for (int j = BW; j < GNX - BW; j += 2) {
                        v = +Fv[0](i, j + 0, k + 0);
                        v += Fv[0](i, j + 1, k + 0);
                        v += Fv[0](i, j + 0, k + 1);
                        v += Fv[0](i, j + 1, k + 1);
                        v *= 0.25;
                        mpi_buffer[f][cnt] = v;
                        cnt++;
                    }
                }
                tag = tag_gen(TAG_CF, get_id(), f);
                assert(cnt==sz);
                p = dynamic_cast<HydroGrid*>(get_parent());
                if (f % 2 == my_child_index().get_x()) {
                    p = dynamic_cast<HydroGrid*>(p->get_sibling(f));
                }
                MPI_Isend(mpi_buffer[f], cnt, MPI_state_t, p->proc(), tag, MPI_COMM_WORLD, send_request + f);
            }
        }
    } else if (dir == 1) {
        for (int f = 2; f < 4; f++) {
            if (is_amr_bound(f)) {
                mpi_buffer[f] = new State[sz];
                i = BW + (f - 2) * (GNX - 2 * BW);
                cnt = 0;
                for (int k = BW; k < GNX - BW; k += 2) {
                    for (int j = BW; j < GNX - BW; j += 2) {
                        v = +Fv[1](j + 0, i, k + 0);
                        v += Fv[1](j + 1, i, k + 0);
                        v += Fv[1](j + 0, i, k + 1);
                        v += Fv[1](j + 1, i, k + 1);
                        v *= 0.25;
                        mpi_buffer[f][cnt] = v;
                        cnt++;
                    }
                }
                tag = tag_gen(TAG_CF, get_id(), f);
                assert(cnt==sz);
                p = dynamic_cast<HydroGrid*>(get_parent());
                if (f % 2 == my_child_index().get_y()) {
                    p = dynamic_cast<HydroGrid*>(p->get_sibling(f));
                }
                MPI_Isend(mpi_buffer[f], cnt, MPI_state_t, p->proc(), tag, MPI_COMM_WORLD, send_request + f);
            }
        }
    } else {
        for (int f = 4; f < 6; f++) {
            if (is_amr_bound(f)) {
                mpi_buffer[f] = new State[sz];
                i = BW + (f - 4) * (GNX - 2 * BW);
                cnt = 0;
                for (int k = BW; k < GNX - BW; k += 2) {
                    for (int j = BW; j < GNX - BW; j += 2) {
                        v = +Fv[2](j + 0, k + 0, i);
                        v += Fv[2](j + 1, k + 0, i);
                        v += Fv[2](j + 0, k + 1, i);
                        v += Fv[2](j + 1, k + 1, i);
                        v *= 0.25;
                        mpi_buffer[f][cnt] = v;
                        cnt++;
                    }
                }
                tag = tag_gen(TAG_CF, get_id(), f);
                assert(cnt==sz);
                p = dynamic_cast<HydroGrid*>(get_parent());
                if (f % 2 == my_child_index().get_z()) {
                    p = dynamic_cast<HydroGrid*>(p->get_sibling(f));
                }
                MPI_Isend(mpi_buffer[f], cnt, MPI_state_t, p->proc(), tag, MPI_COMM_WORLD, send_request + f);
            }
        }
    }
    inc_instruction_pointer(dir);
}

void HydroGrid::flux_cf_adjust_recv(int dir) {
    const int sz = (GNX - 2 * BW) * (GNX - 2 * BW) / 4;
    int f, tag;
    ChildIndex ci;
    OctNode* child;
    amr_has[dir] = false;
    for (ci = 0; ci < 8; ci++) {
        child = (get_child(ci));
        if (child != NULL) {
            f = (2 * dir) ^ 1 ^ (ci.vec()[dir] % 2);
            mpi_amr_buffer[f][ci] = NULL;
            if (child->is_amr_bound(f)) {
                mpi_amr_buffer[f][ci] = new State[sz];
                tag = tag_gen(TAG_CF, child->get_id(), f);
                MPI_Irecv(mpi_amr_buffer[f][ci], sz, MPI_state_t, child->proc(), tag, MPI_COMM_WORLD, amr_recv_request[f] + ci);
                amr_has[dir] = true;
            }
        } else {
            f = (2 * dir) ^ (ci.vec()[dir] % 2);
            mpi_amr_buffer[f][ci] = NULL;
            if (is_real_bound(f)) {
                child = get_sibling(f)->get_child(ci ^ (1 << dir));
                if (child != NULL) {
                    if (child->is_amr_bound(f ^ 1)) {
                        mpi_amr_buffer[f][ci] = new State[sz];
                        tag = tag_gen(TAG_CF, child->get_id(), f ^ 1);
                        MPI_Irecv(mpi_amr_buffer[f][ci], sz, MPI_state_t, child->proc(), tag, MPI_COMM_WORLD, amr_recv_request[f] + ci);
                        amr_has[dir] = true;
                    }
                }
            }
        }
    }
    inc_instruction_pointer(dir);
}

void HydroGrid::flux_cf_adjust_send_wait(int dir) {
    int flag;
    MPI_Testall(2, send_request + 2 * dir, &flag, MPI_STATUS_IGNORE );
    if (flag) {
        for (int f = 2 * dir; f < 2 * dir + 2; f++) {
            if (is_amr_bound(f)) {
                delete[] mpi_buffer[f];
            }
        }
        inc_instruction_pointer(dir);
    }
}

void HydroGrid::flux_cf_adjust_recv_wait(int dir) {
    int flag, i, flag_odd, flag_even;
    HydroGrid* child;
    int xlb, xub, ylb, yub, zlb, zub, cnt, f;
    bool test;
    if (amr_has[dir]) {
        MPI_Testall(8, amr_recv_request[2 * dir], &flag_odd, MPI_STATUS_IGNORE );
        MPI_Testall(8, amr_recv_request[2 * dir + 1], &flag_even, MPI_STATUS_IGNORE );
        flag = flag_odd && flag_even;
        if (flag) {
            for (ChildIndex ci = 0; ci < 8; ci++) {
                child = dynamic_cast<HydroGrid*>(get_child(ci));
                test = false;
                if (child != NULL) {
                    f = (2 * dir) ^ 1 ^ (ci.vec()[dir] % 2);
                    if (child->is_amr_bound(f)) {
                        test = true;
                    }
                } else {
                    f = (2 * dir) ^ (ci.vec()[dir] % 2);
                    if (is_real_bound(f)) {
                        child = dynamic_cast<HydroGrid*>(get_sibling(f)->get_child(ci ^ (1 << dir)));
                        if (child != NULL) {
                            if (child->is_amr_bound(f ^ 1)) {
                                test = true;
                            }
                        }
                    }
                }
                if (test) {
                    zlb = BW + ci.get_z() * (GNX / 2 - BW);
                    zub = BW + (ci.get_z() + 1) * (GNX / 2 - BW) - 1;
                    ylb = BW + ci.get_y() * (GNX / 2 - BW);
                    yub = BW + (ci.get_y() + 1) * (GNX / 2 - BW) - 1;
                    xlb = BW + ci.get_x() * (GNX / 2 - BW);
                    xub = BW + (ci.get_x() + 1) * (GNX / 2 - BW) - 1;
                    cnt = 0;
                    if (dir == 0) {
                        i = BW + ((f % 2) + ci.get_x()) * (GNX / 2 - BW);
                        for (int k = zlb; k <= zub; k++) {
                            for (int j = ylb; j <= yub; j++) {
                                assert( mpi_amr_buffer[f][ci]);
                                Fv[dir](i, j, k) = mpi_amr_buffer[f][ci][cnt];
                                cnt++;
                            }
                        }
                    } else if (dir == 1) {
                        i = BW + ((f % 2) + ci.get_y()) * (GNX / 2 - BW);
                        for (int k = zlb; k <= zub; k++) {
                            for (int j = xlb; j <= xub; j++) {
                                Fv[dir](j, i, k) = mpi_amr_buffer[f][ci][cnt];
                                cnt++;
                            }
                        }
                    } else {
                        i = BW + ((f % 2) + ci.get_z()) * (GNX / 2 - BW);
                        for (int k = ylb; k <= yub; k++) {
                            for (int j = xlb; j <= xub; j++) {
                                Fv[dir](j, k, i) = mpi_amr_buffer[f][ci][cnt];
                                cnt++;
                            }
                        }
                    }assert( cnt == (GNX-2*BW)*(GNX-2*BW)/4);
                    delete[] mpi_amr_buffer[f][ci];
                }
            }
        }
    } else {
        flag = true;
    }
    if (flag) {
        inc_instruction_pointer(dir);
    }
}

void HydroGrid::redistribute_send() {
    int tag, j, k, l, cnt;
    if (proc() != next_proc()) {
        tag = tag_gen(TAG_MOVE, this->get_id(), 0);
        cnt = 0;
        mpi_buffer[0] = new State[(GNX - 2 * BW) * (GNX - 2 * BW) * (GNX - 2 * BW)];
        for (l = BW; l < GNX - BW; l++) {
            for (k = BW; k < GNX - BW; k++) {
                for (j = BW; j < GNX - BW; j++) {
                    mpi_buffer[0][cnt] = U(j, k, l);
                    cnt++;
                }
            }
        }
        MPI_Isend(mpi_buffer[0], cnt, MPI_state_t, next_proc(), tag, MPI_COMM_WORLD, send_request + 0);
    }
}

void HydroGrid::redistribute_recv() {
    if (proc() != last_proc()) {
        int tag = tag_gen(TAG_MOVE, this->get_id(), 0);
        MPI_Irecv(U.ptr(), 1, MPI_interior_t, last_proc(), tag, MPI_COMM_WORLD, recv_request + 0);
    }
}

void HydroGrid::error_from_parent_send(int dir) {
    if (dir == 0 && shadow) {
        OctNode* child;
        int cnt, tag, sz;
        Vector<int, 3> lb, ub, o;
        for (ChildIndex ci = 0; ci < 8; ci++) {
            child = get_child(ci);
            if (child != NULL) {
                cnt = 0;
                sz = (GNX / 2 - BW) * (GNX / 2 - BW) * (GNX / 2 - BW);
                mpi_buffer[ci] = new State[sz];
                lb = BW;
                ub = GNX / 2 - 1;
                o = ci.vec();
                o *= GNX / 2 - BW;
                lb += o;
                ub += o;
                for (int k = lb[2]; k <= ub[2]; k++) {
                    for (int j = lb[1]; j <= ub[1]; j++) {
                        for (int i = lb[0]; i <= ub[0]; i++) {
                            mpi_buffer[ci][cnt] = D(i, j, k);
                            cnt++;
                        }
                    }
                }assert(cnt==sz);
                tag = tag_gen(TAG_ERROR, child->get_id(), 0);
                MPI_Isend(mpi_buffer[ci], cnt, MPI_state_t, child->proc(), tag, MPI_COMM_WORLD, send_request + ci);
            }
        }
    }
    inc_instruction_pointer(dir);
}

void HydroGrid::error_from_parent_recv(int dir) {
    if (dir == 0 && shadow) {
        int tag, sz;
        State u;
        if (get_level() != 0) {
            tag = tag_gen(TAG_ERROR, get_id(), 0);
            sz = (GNX / 2 - BW) * (GNX / 2 - BW) * (GNX / 2 - BW);
            mpi_buffer[0] = new State[sz];
            MPI_Irecv(mpi_buffer[0], sz, MPI_state_t, get_parent()->proc(), tag, MPI_COMM_WORLD, recv_request);
        }
    }
    inc_instruction_pointer(dir);
}

void HydroGrid::error_from_parent_send_wait(int dir) {
    int flag;
    if (dir == 0 && shadow) {
        MPI_Testall(8, send_request, &flag, MPI_STATUS_IGNORE );
        if (flag) {
            for (int ci = 0; ci < 8; ci++) {
                if (get_child(ci) != NULL) {
                    delete[] mpi_buffer[ci];
                }
            }
        }
    } else {
        flag = true;
    }
    if (flag) {
        inc_instruction_pointer(dir);
    }

}

void HydroGrid::error_from_parent_recv_wait(int dir) {
    int flag, cnt;
    State u;
    if (dir == 0 && shadow) {
        MPI_Test(recv_request, &flag, MPI_STATUS_IGNORE );
        if (flag && get_level() != 0) {
            cnt = 0;
            for (int k = BW; k < GNX - BW; k += 2) {
                for (int j = BW; j < GNX - BW; j += 2) {
                    for (int i = BW; i < GNX - BW; i += 2) {
                        u = mpi_buffer[0][cnt++];
                        E(i + 0, j + 0, k + 0) = u - D(i + 0, j + 0, k + 0);
                        E(i + 1, j + 0, k + 0) = u - D(i + 1, j + 0, k + 0);
                        E(i + 0, j + 1, k + 0) = u - D(i + 0, j + 1, k + 0);
                        E(i + 1, j + 1, k + 0) = u - D(i + 1, j + 1, k + 0);
                        E(i + 0, j + 0, k + 1) = u - D(i + 0, j + 0, k + 1);
                        E(i + 1, j + 0, k + 1) = u - D(i + 1, j + 0, k + 1);
                        E(i + 0, j + 1, k + 1) = u - D(i + 0, j + 1, k + 1);
                        E(i + 1, j + 1, k + 1) = u - D(i + 1, j + 1, k + 1);

                    }
                }
            }assert(cnt==(GNX/2-BW)*(GNX/2-BW)*(GNX/2-BW));
            delete[] mpi_buffer[0];
        }
    } else {
        flag = true;
    }
    if (flag) {
        inc_instruction_pointer(dir);
    }

}

void HydroGrid::add_to_dudt(int i, int j, int k, const Vector<Real, STATE_NF>& s) {
    D(i, j, k) += s;
}

Vector<Real, STATE_NF> HydroGrid::get_dudt(int i, int j, int k) const {
    return D(i, j, k);
}

void HydroGrid::compute_dudt(int dir) {
    if (dir == 0) {
        const Real da = dx * dx;
        _3Vec x;
        int k, j, i;
        for (k = BW; k < GNX - BW; k++) {
            for (j = BW; j < GNX - BW; j++) {
                for (i = BW; i < GNX - BW; i++) {
                    x = HydroGrid::X(i, j, k);
                    D(i, j, k) = 0.0;
                    D(i, j, k) += -(Fv[0](i + 1, j, k) - Fv[0](i, j, k)) / dx;
                    D(i, j, k) += -(Fv[1](i, j + 1, k) - Fv[1](i, j, k)) / dx;
                    D(i, j, k) += -(Fv[2](i, j, k + 1) - Fv[2](i, j, k)) / dx;
                    D(i, j, k) += -State::x_scalar_flux_coeff(x) * (Fs[0](i + 1, j, k) - Fs[0](i, j, k)) / dx;
                    D(i, j, k) += -State::y_scalar_flux_coeff(x) * (Fs[1](i, j + 1, k) - Fs[1](i, j, k)) / dx;
                    D(i, j, k) += -State::z_scalar_flux_coeff(x) * (Fs[2](i, j, k + 1) - Fs[2](i, j, k)) / dx;
                    D(i, j, k) += U(i, j, k).source(this->X(i, j, k), get_time());
                }
            }
        }
        this->compute_flow_off();
    }
    inc_instruction_pointer(dir);
}

void HydroGrid::compute_flow_off() {
    const Real da = dx * dx;
    int k, j, i;
    for (k = BW; k < GNX - BW; k++) {
        for (j = BW; j < GNX - BW; j++) {
            for (i = BW; i < GNX - BW; i++) {
                if (!zone_is_refined(i, j, k)) {
                    _3Vec x = HydroGrid::X(i, j, k);
                    if (is_phys_bound(XL) && i == BW) {
                        DFO -= (Fv[0](i, j, k)) * da;
                    }
                    if (is_phys_bound(XU) && i == GNX - BW - 1) {
                        DFO += (Fv[0](i + 1, j, k)) * da;
                    }
                    if (is_phys_bound(YL) && j == BW) {
                        DFO -= (Fv[1](i, j, k)) * da;
                    }
                    if (is_phys_bound(YU) && j == GNX - BW - 1) {
                        DFO += (Fv[1](i, j + 1, k)) * da;
                    }
                    if (is_phys_bound(ZL) && k == BW) {
                        DFO -= (Fv[2](i, j, k)) * da;
                    }
                    if (is_phys_bound(ZU) && k == GNX - BW - 1) {
                        DFO += (Fv[2](i, j, k + 1)) * da;
                    }
                    if (is_phys_bound(XL) && i == BW) {
                        DFO -= State::x_scalar_flux_coeff(x) * (Fs[0](i, j, k)) * da;
                    }
                    if (is_phys_bound(XU) && i == GNX - BW - 1) {
                        DFO += State::x_scalar_flux_coeff(x) * (Fs[0](i + 1, j, k)) * da;
                    }
                    if (is_phys_bound(YL) && j == BW) {
                        DFO -= State::y_scalar_flux_coeff(x) * (Fs[1](i, j, k)) * da;
                    }
                    if (is_phys_bound(YU) && j == GNX - BW - 1) {
                        DFO += State::y_scalar_flux_coeff(x) * (Fs[1](i, j + 1, k)) * da;
                    }
                    if (is_phys_bound(ZL) && k == BW) {
                        DFO -= State::z_scalar_flux_coeff(x) * (Fs[2](i, j, k)) * da;
                    }
                    if (is_phys_bound(ZU) && k == GNX - BW - 1) {
                        DFO += State::z_scalar_flux_coeff(x) * (Fs[2](i, j, k + 1)) * da;
                    }
                }
            }
        }
    }
}

void HydroGrid::compute_update(int dir) {
    if (dir == 0) {
        int k, j, i;
        for (k = BW; k < GNX - BW; k++) {
            for (j = BW; j < GNX - BW; j++) {
                for (i = BW; i < GNX - BW; i++) {
                    U(i, j, k) = (U(i, j, k) + D(i, j, k) * _dt) * _beta + U0(i, j, k) * (1.0 - _beta);
                    if (shadow) {
                        E0(i, j, k) = E0(i, j, k) * (1.0 - _beta) + E(i, j, k) * _beta;
                    }
                    U(i, j, k).floor(X(i, j, k));

                }
            }
        }
    }
    inc_instruction_pointer(dir);
}

void HydroGrid::enforce_dual_energy_formalism(int dir) {
    if (dir == 0) {
        HydroGrid* g = this;
        for (int k = BW; k < GNX - BW; k++) {
            for (int j = BW; j < GNX - BW; j++) {
                for (int i = BW; i < GNX - BW; i++) {
                    g->U(i, j, k).enforce_dual_energy_formalism(g->X(i, j, k), g->U(i + 1, j, k), g->U(i - 1, j, k), g->U(i, j + 1, k), g->U(i, j - 1, k),
                            g->U(i, j, k + 1), g->U(i, j, k - 1));
                }
            }
        }
    }
    inc_instruction_pointer(dir);
}

void HydroGrid::sync(int) {
    if (threads_are_synced(3)) {
        inc_instruction_pointer(0);
        inc_instruction_pointer(1);
        inc_instruction_pointer(2);
    }
}

void HydroGrid::set_refine_flags() {
    if (!shadow) {
        printf("set_refine_flags not defined\n");
        abort();
    }
    ChildIndex c;
    if (get_level() < 1) {
        for (int i = 0; i < OCT_NCHILD; i++) {
            set_refine_flag(i, true);
        }
    } else if (get_level() < get_max_level_allowed()) {
        for (int k = BW; k < GNX - BW; k++) {
            c.set_z(2 * k / GNX);
            for (int j = BW; j < GNX - BW; j++) {
                c.set_y(2 * j / GNX);
                for (int i = BW; i < GNX - BW; i++) {
                    c.set_x(2 * i / GNX);
                    if (!get_refine_flag(c)) {
                        if (E0(i, j, k).mag() > 10000) {
                            set_refine_flag(c, true);
                        }
                    }
                }
            }
        }
    }
}

_3Vec HydroGrid::X(int i, int j, int k) const {
    _3Vec x;
    x[0] = HydroGrid::xc(i);
    x[1] = HydroGrid::yc(j);
    x[2] = HydroGrid::zc(k);
    return x;

}

_3Vec HydroGrid::Xfx(int i, int j, int k) const {
    _3Vec x;
    x[0] = HydroGrid::xf(i);
    x[1] = HydroGrid::yc(j);
    x[2] = HydroGrid::zc(k);
    return x;

}

_3Vec HydroGrid::Xfy(int i, int j, int k) const {
    _3Vec x;
    x[0] = HydroGrid::xc(i);
    x[1] = HydroGrid::yf(j);
    x[2] = HydroGrid::zc(k);
    return x;

}

_3Vec HydroGrid::Xfz(int i, int j, int k) const {
    _3Vec x;
    x[0] = HydroGrid::xc(i);
    x[1] = HydroGrid::yc(j);
    x[2] = HydroGrid::zf(k);
    return x;

}

Reconstruct HydroGrid::reconstruct;

State HydroGrid::operator()(int i, int j, int k) const {
    return U(i, j, k);
}

State& HydroGrid::operator()(int i, int j, int k) {
    return U(i, j, k);
}

Real HydroGrid::xc(int i) const {
    return HydroGrid::xf(i) + 0.5 * dx;
}

Real HydroGrid::yc(int j) const {
    return HydroGrid::yf(j) + 0.5 * dx;
}

Real HydroGrid::zc(int k) const {
    return HydroGrid::zf(k) + 0.5 * dx;
}

Real HydroGrid::xf(int i) const {
    return Real(offset[0] + i) * dx - (GNX / 2) * h0 - origin[0];
}

Real HydroGrid::yf(int i) const {
    return Real(offset[1] + i) * dx - (GNX / 2) * h0 - origin[1];
}

Real HydroGrid::zf(int i) const {
    return Real(offset[2] + i) * dx - (GNX / 2) * h0 - origin[2];
}

void HydroGrid::max_dt_compute(int dir) {
    Real this_dt;
    Real dtinv = 0.0;
    _3Vec x;
    State q0[GNX], ql[GNX], qr[GNX];
    int k, j, i;
    Real this_dtinv;
    for (k = BW; k < GNX - BW; k++) {
        for (j = BW; j < GNX - BW; j++) {
            this_dtinv = 0.0;
            if (dir == 0) {
                for (i = 0; i < GNX; i++) {
                    q0[i] = U(i, j, k);
                    q0[i].to_prim(HydroGrid::X(i, j, k));
                }
                reconstruct(q0, ql, qr);
                for (i = BW; i < GNX - BW + 1; i++) {
                    x = HydroGrid::Xfx(i, j, k);
                    ql[i].from_prim(x);
                    qr[i].from_prim(x);
                    this_dtinv = max(this_dtinv, ql[i].max_abs_x_eigen(x), qr[i].max_abs_x_eigen(x));
                }
            } else if (dir == 1) {
                for (i = 0; i < GNX; i++) {
                    q0[i] = U(j, i, k);
                    q0[i].to_prim(HydroGrid::X(j, i, k));
                }
                reconstruct(q0, ql, qr);
                for (i = BW; i < GNX - BW + 1; i++) {
                    x = HydroGrid::Xfy(j, i, k);
                    ql[i].from_prim(x);
                    qr[i].from_prim(x);
                    this_dtinv = max(this_dtinv, ql[i].max_abs_y_eigen(x), qr[i].max_abs_y_eigen(x));
                }
            } else {
                for (i = 0; i < GNX; i++) {
                    q0[i] = U(j, k, i);
                    q0[i].to_prim(HydroGrid::X(j, k, i));
                }
                reconstruct(q0, ql, qr);
                for (i = BW; i < GNX - BW + 1; i++) {
                    x = HydroGrid::Xfz(j, k, i);
                    ql[i].from_prim(x);
                    qr[i].from_prim(x);
                    this_dtinv = max(this_dtinv, ql[i].max_abs_z_eigen(x), qr[i].max_abs_z_eigen(x));
                }
            }
            dtinv = max(this_dtinv, dtinv);
        }
    }
    dtinv /= dx;
    if (dtinv == 0.0) {
        dtinv = 1.0E-10;
    }
    this_dt = 1.0 / dtinv;
    this_dt *= GRID_CFL_FACTOR;
    eax = min(this_dt, eax);
    inc_instruction_pointer(dir);
}

void HydroGrid::reduce_dt(Real dt) {
    eax = min(eax, dt);
}

void HydroGrid::physical_boundary(int dir) {
    Vector<int, 3> lb, ub, k;
    if (is_phys_bound(2 * dir + 0)) {
        lb = BW;
        ub = GNX - BW - 1;
        lb[dir] = 0;
        ub[dir] = BW - 1;
        for (Indexer3d i(lb, ub); !i.end(); i++) {
            k = i;
            k[dir] = BW;
            U(i[0], i[1], i[2]) = U(k[0], k[1], k[2]);
        }
    }
    if (is_phys_bound(2 * dir + 1)) {
        lb = BW;
        ub = GNX - BW - 1;
        lb[dir] = GNX - BW;
        ub[dir] = GNX - 1;
        for (Indexer3d i(lb, ub); !i.end(); i++) {
            k = i;
            k[dir] = GNX - BW - 1;
            U(i[0], i[1], i[2]) = U(k[0], k[1], k[2]);
        }
    }
    inc_instruction_pointer(dir);
}

#endif
