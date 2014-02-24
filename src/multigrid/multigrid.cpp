#include "multigrid.h"
#include <math.h>
#include <stdlib.h>

#define RELAX_RESID 0.5
#define VDOWN_RESID 1.0
#define VUP_RESID 0.05

#include "../indexer3d.h"
#include "../tag.h"
#include "../legendre.h"


void MultiGrid::write(FILE* fp) const {
	fwrite(&dx, sizeof(Real), 1, fp);
	if (get_level() == 0) {
		fwrite(&h0, sizeof(Real), 1, fp);
	}
	if (proc() == MPI_rank()) {
		fwrite(phi.ptr(), sizeof(Real), PNX * PNX * PNX, fp);
	}
}

void MultiGrid::read(FILE* fp) {
	fread(&dx, sizeof(Real), 1, fp);
	if (get_level() == 0) {
		fread(&h0, sizeof(Real), 1, fp);
	}
	if (proc() == MPI_rank()) {
		fread(phi.ptr(), sizeof(Real), PNX * PNX * PNX, fp);
	}
}

Real MultiGrid::get_dx() const {
	return dx;
}

Real MultiGrid::get_residual(int i, int j, int k) const {
	return dphi(i, j, k);
}

Real MultiGrid::get_phi(int i, int j, int k) const {
	return phi(i, j, k);
}

Real* MultiGrid::get_phi_ptr(int i, int j, int k)  {
	return phi.ptr(i, j, k);
}

Real MultiGrid::get_source(int i, int j, int k) const {
	return S(i, j, k);
}

Real MultiGrid::get_fx(int i, int j, int k) const {
	return fx(i, j, k);
}

Real MultiGrid::get_fy(int i, int j, int k) const {
	return fy(i, j, k);
}

Real MultiGrid::get_fz(int i, int j, int k) const {
	return fz(i, j, k);
}

void MultiGrid::set_source(int i, int j, int k, Real s) {
	S(i, j, k) = s;
}

void MultiGrid::create_multigrid_child(const ChildIndex& c) {
	Vector<int, 3> child_offset;
	MultiGrid* child;
	child = dynamic_cast<MultiGrid*>(get_child(c));
	child_offset = (offset * 2 + 1) + c.vector() * (PNX - 2);
	child->dx = h0 / Real(1 << (get_level() + 1));
	child->offset = child_offset;
	child->inject_from_parent(c);
}

void MultiGrid::create_child(const ChildIndex& c) {
	OctNode::create_child(c);
	this->create_multigrid_child(c);
}

int MultiGrid::guard_proc(OctFace f) const {
	if (is_real_bound(f)) {
		return get_sibling(f)->proc();
	} else if (is_phys_bound(f)) {
		return proc();
	} else {
		if (((my_child_index().vec()[f / 2]) ^ (f % 2)) != 1) {
			return get_parent()->get_sibling(f)->proc();
		} else {
			return get_parent()->proc();
		}
	}
}

void MultiGrid::dphi_amr_boundary_communicate(int) {
	if (amr_cnt > 0) {
		int cnt, tag_send;
		for (int i = 0; i < amr_cnt; i++) {
			int sz = (PNX - 2) * (PNX - 2);
			mpi_buffer[i] = new Real[sz];
			cnt = 0;
			for (int l = amr_lb[i][2]; l <= amr_ub[i][2]; l++) {
				for (int k = amr_lb[i][1]; k <= amr_ub[i][1]; k++) {
					for (int j = amr_lb[i][0]; j <= amr_ub[i][0]; j++) {
						mpi_buffer[i][cnt] = dphi(j / 2, k / 2, l / 2);
						++cnt;
					}
				}
			}
			tag_send = tag_gen(TAG_RELAX, amr_id[i], amr_face[i], 0);
			MPI_Isend(mpi_buffer[i], cnt, MPI_DOUBLE_PRECISION, amr_child_proc[i], tag_send, MPI_COMM_WORLD, amr_request + i);
		}
	}
	inc_instruction_pointer();
}

void MultiGrid::dphi_amr_boundary_wait(int) {
	int rc;
	MPI_Testall(amr_cnt, amr_request, &rc, MPI_STATUS_IGNORE);
	if (rc) {
		for (int i = 0; i < amr_cnt; i++) {
			delete[] mpi_buffer[i];
		}
		inc_instruction_pointer();
	}
}

void MultiGrid::forces_compute(int) {
	for (int k = 1; k < PNX; k++) {
		for (int j = 1; j < PNX; j++) {
			for (int i = 1; i < PNX; i++) {
				fx(i, j, k) = -(phi(i, j, k) - phi(i - 1, j, k)) / dx;
				fy(i, j, k) = -(phi(i, j, k) - phi(i, j - 1, k)) / dx;
				fz(i, j, k) = -(phi(i, j, k) - phi(i, j, k - 1)) / dx;
			}
		}
	}
	inc_instruction_pointer();
}

void MultiGrid::forces_adjust_send(int) {
	Vector<int, 3> lb, ub;
	Real v;
	int i, cnt, tag, dir;
	for (int face = 0; face < 6; face++) {
		dir = face / 2;
		if (is_amr_bound(face)) {
			mpi_buffer[face] = new Real[PNX * PNX / 4];
			i = 1 + (PNX - 2) * (face & 1);
			cnt = 0;
			for (int k = 1; k < PNX - 1; k += 2) {
				for (int j = 1; j < PNX - 1; j += 2) {
					if (dir == 0) {
						v = +fx(i, j + 0, k + 0);
						v += fx(i, j + 0, k + 1);
						v += fx(i, j + 1, k + 0);
						v += fx(i, j + 1, k + 1);
					} else if (dir == 1) {
						v = +fy(j + 0, i, k + 0);
						v += fy(j + 0, i, k + 1);
						v += fy(j + 1, i, k + 0);
						v += fy(j + 1, i, k + 1);
					} else {
						v = +fz(j + 0, k + 0, i);
						v += fz(j + 0, k + 1, i);
						v += fz(j + 1, k + 0, i);
						v += fz(j + 1, k + 1, i);
					}
					v *= 0.25;
					mpi_buffer[face][cnt] = v;
					cnt++;
				}
			}
			assert(cnt == (PNX-2)*(PNX-2)/4);
			tag = tag_gen(TAG_FORCE, get_id(), face);
			MPI_Isend(mpi_buffer[face], cnt, MPI_DOUBLE_PRECISION, guard_proc(face), tag, MPI_COMM_WORLD, send_request + face);
		} else {
			send_request[face] = MPI_REQUEST_NULL;
		}
	}
	inc_instruction_pointer();
}

void MultiGrid::forces_adjust_send_wait(int) {
	int flag;
	MPI_Testall(6, send_request, &flag, MPI_STATUS_IGNORE);
	if (flag) {
		for (int i = 0; i < 6; i++) {
			if (is_amr_bound(i)) {
				delete[] mpi_buffer[i];
			}
		}
	}
	if (flag) {
		inc_instruction_pointer();
	}
}

void MultiGrid::forces_adjust_recv(int) {
	if (amr_cnt > 0) {
		int tag;
		for (int i = 0; i < amr_cnt; i++) {
			const int sz = (PNX - 2) * (PNX - 2) / 4;
			mpi_buffer[i] = new Real[sz];
			tag = tag_gen(TAG_FORCE, amr_id[i], amr_face[i]);
			MPI_Irecv(mpi_buffer[i], sz, MPI_DOUBLE_PRECISION, amr_child_proc[i], tag, MPI_COMM_WORLD, amr_request + i);
		}
	}
	inc_instruction_pointer();
}

void MultiGrid::forces_adjust_recv_wait(int) {
	int flag;
	if (amr_cnt > 0) {
		MPI_Testall(amr_cnt, amr_request, &flag, MPI_STATUS_IGNORE);
		if (flag) {
			int cnt, dir;
			int i;
			for (int ci = 0; ci < amr_cnt; ci++) {
				dir = amr_face[ci] / 2;
				cnt = 0;
				i = amr_lb[ci][dir] / 2;
				if (i == PNX / 2 - 1) {
					i++;
				} else if (i == PNX - 2) {
					i = PNX - 1;
				}
				if (dir == 0) {
					for (int k = amr_lb[ci][2] / 2; k <= amr_ub[ci][2] / 2; k++) {
						for (int j = amr_lb[ci][1] / 2; j <= amr_ub[ci][1] / 2; j++) {
							fx(i, j, k) = mpi_buffer[ci][cnt];
							++cnt;
						}
					}
				} else if (dir == 1) {
					for (int k = amr_lb[ci][2] / 2; k <= amr_ub[ci][2] / 2; k++) {
						for (int j = amr_lb[ci][0] / 2; j <= amr_ub[ci][0] / 2; j++) {
							fy(j, i, k) = mpi_buffer[ci][cnt];
							++cnt;
						}
					}
				} else {
					for (int k = amr_lb[ci][1] / 2; k <= amr_ub[ci][1] / 2; k++) {
						for (int j = amr_lb[ci][0] / 2; j <= amr_ub[ci][0] / 2; j++) {
							fz(j, k, i) = mpi_buffer[ci][cnt];
							++cnt;
						}
					}
				}
				delete[] mpi_buffer[ci];
			}
		}
	} else {
		flag = true;
	}
	if (flag) {
		inc_instruction_pointer();
	}
}
void MultiGrid::phi_calc_amr_bounds() {
	ChildIndex cl, cr;
	const OctNode *child;
	int cj, ck;
	int i;
	amr_cnt = 0;
	Vector<int, 3> lb, ub;
	for (int l = 0; l < 4; l++) {
		for (cj = 0; cj < 2; cj++) {
			for (ck = 0; ck < 2; ck++) {
				cr.set_index(1, cj, ck);
				cl.set_index(0, cj, ck);
				child = NULL;
				switch (l) {
				case 0:
					if (get_child(cl) == NULL && get_sibling(XL)) {
						if (get_sibling(XL)->get_child(cr) != NULL) {
							child = get_sibling(XL)->get_child(cr);
						}
					}
					i = 2;
					amr_face[amr_cnt] = XU;
					break;
				case 1:
					if (get_child(cl) != NULL && get_child(cr) == NULL) {
						child = get_child(cl);
					}
					amr_face[amr_cnt] = XU;
					i = PNX;
					break;
				case 2:
					if (get_child(cl) == NULL && get_child(cr) != NULL) {
						child = get_child(cr);
					}
					amr_face[amr_cnt] = XL;
					i = PNX - 1;
					break;
				case 3:
					if (get_child(cr) == NULL && get_sibling(XU)) {
						if (get_sibling(XU)->get_child(cl) != NULL) {
							child = get_sibling(XU)->get_child(cl);
						}
					}
					i = 2 * PNX - 3;
					amr_face[amr_cnt] = XL;
					break;
				}
				if (child != NULL) {
					lb[0] = ub[0] = i;
					lb[1] = 2 + cj * (PNX - 2);
					ub[1] = lb[1] + PNX - 3;
					lb[2] = 2 + ck * (PNX - 2);
					ub[2] = lb[2] + PNX - 3;
					amr_lb[amr_cnt] = lb;
					amr_ub[amr_cnt] = ub;
					amr_child_proc[amr_cnt] = child->proc();
					amr_id[amr_cnt] = child->get_id();
					amr_cnt++;
				}
			}
		}
	}
	for (int l = 0; l < 4; l++) {
		for (cj = 0; cj < 2; cj++) {
			for (ck = 0; ck < 2; ck++) {
				cr.set_index(cj, 1, ck);
				cl.set_index(cj, 0, ck);
				child = NULL;
				switch (l) {
				case 0:
					if (get_child(cl) == NULL && get_sibling(YL)) {
						if (get_sibling(YL)->get_child(cr) != NULL) {
							child = get_sibling(YL)->get_child(cr);
						}
					}
					i = 2;
					amr_face[amr_cnt] = YU;
					break;
				case 1:
					if (get_child(cl) != NULL && get_child(cr) == NULL) {
						child = get_child(cl);
					}
					i = PNX;
					amr_face[amr_cnt] = YU;
					break;
				case 2:
					if (get_child(cl) == NULL && get_child(cr) != NULL) {
						child = get_child(cr);
					}
					i = PNX - 1;
					amr_face[amr_cnt] = YL;
					break;
				case 3:
					if (get_child(cr) == NULL && get_sibling(YU)) {
						if (get_sibling(YU)->get_child(cl) != NULL) {
							child = get_sibling(YU)->get_child(cl);
						}
					}
					i = 2 * PNX - 3;
					amr_face[amr_cnt] = YL;
					break;
				}
				if (child != NULL) {
					lb[1] = ub[1] = i;
					lb[0] = 2 + cj * (PNX - 2);
					ub[0] = lb[0] + PNX - 3;
					lb[2] = 2 + ck * (PNX - 2);
					ub[2] = lb[2] + PNX - 3;
					amr_lb[amr_cnt] = lb;
					amr_ub[amr_cnt] = ub;
					amr_child_proc[amr_cnt] = child->proc();
					amr_id[amr_cnt] = child->get_id();
					amr_cnt++;
				}
			}
		}
	}
	for (int l = 0; l < 4; l++) {
		for (cj = 0; cj < 2; cj++) {
			for (ck = 0; ck < 2; ck++) {
				cr.set_index(cj, ck, 1);
				cl.set_index(cj, ck, 0);
				child = NULL;
				switch (l) {
				case 0:
					if (get_child(cl) == NULL && get_sibling(ZL) != NULL) {
						if (get_sibling(ZL)->get_child(cr) != NULL) {
							child = get_sibling(ZL)->get_child(cr);
						}
					}
					i = 2;
					amr_face[amr_cnt] = ZU;
					break;
				case 1:
					if (get_child(cl) != NULL && get_child(cr) == NULL) {
						child = get_child(cl);
					}
					amr_face[amr_cnt] = ZU;
					i = PNX;
					break;
				case 2:
					if (get_child(cl) == NULL && get_child(cr) != NULL) {
						child = get_child(cr);
					}
					amr_face[amr_cnt] = ZL;
					i = PNX - 1;
					break;
				case 3:
					if (get_child(cr) == NULL && get_sibling(ZU) != NULL) {
						if (get_sibling(ZU)->get_child(cl) != NULL) {
							child = get_sibling(ZU)->get_child(cl);
						}
					}
					amr_face[amr_cnt] = ZL;
					i = 2 * PNX - 3;
					break;
				}
				if (child != NULL) {
					lb[0] = 2 + cj * (PNX - 2);
					ub[0] = lb[0] + PNX - 3;
					lb[1] = 2 + ck * (PNX - 2);
					ub[1] = lb[1] + PNX - 3;
					lb[2] = ub[2] = i;
					amr_lb[amr_cnt] = lb;
					amr_ub[amr_cnt] = ub;
					amr_child_proc[amr_cnt] = child->proc();
					amr_id[amr_cnt] = child->get_id();
					amr_cnt++;
				}
			}
		}
	}

}

void MultiGrid::phi_children_recv(int) {
	int tag;
	MultiGrid* child;
	for (int i = 0; i < 8; i++) {
		child = dynamic_cast<MultiGrid*>(get_child(i));
		if (child != NULL) {
			tag = tag_gen(TAG_CHILD, get_id(), i);
			MPI_Irecv(phi.ptr(), 1, MPI_recv_amr_child_t[i], get_child(i)->proc(), tag, MPI_COMM_WORLD, recv_request + i);
		} else {
			recv_request[i] = MPI_REQUEST_NULL;
		}
	}
	inc_instruction_pointer();
}

void MultiGrid::phi_children_recv_wait(int) {
	int rc;
	MPI_Testall(8, recv_request, &rc, MPI_STATUS_IGNORE);
	if (rc) {
		inc_instruction_pointer();
	}
}

void MultiGrid::phi_children_send(int) {
	if (get_level() != 0) {
		int tag;
		int cnt = 0;
		/*	mpi_buffer[0] = new Real[3 * (PNX - 2) * (PNX - 2) / 2];
		 for (int k = 1; k < PNX - 1; k += 2) {
		 for (int j = 1; j < PNX - 1; j += 2) {
		 for (int i = 1; i < PNX - 1; i += 2) {*/
		mpi_buffer[0] = new Real[(PNX - 2) * (PNX - 2) * (PNX - 2) / 8];
		for (int k = 1; k < PNX - 1; k += 2) {
			for (int j = 1; j < PNX - 1; j += 2) {
				for (int i = 1; i < PNX - 1; i += 2) {
					//		if ((k == 1) || (k == PNX - 3) || (j == 1) || (j == PNX - 3) || (i == 1) || (i == PNX - 3)) {
					mpi_buffer[0][cnt] = phi.oct_avg(i, j, k);
					cnt++;
					//	}
				}
			}

		}
		tag = tag_gen(TAG_CHILD, get_parent()->get_id(), my_child_index());
		if (cnt != 0) {
			MPI_Isend(mpi_buffer[0], cnt, MPI_DOUBLE_PRECISION, get_parent()->proc(), tag, MPI_COMM_WORLD, send_request);
		} else {
			send_request[0] = MPI_REQUEST_NULL;
		}
	}
	inc_instruction_pointer();
}

void MultiGrid::phi_children_send_wait(int) {
	int rc;
	if (get_level() != 0) {
		MPI_Test(send_request, &rc, MPI_STATUS_IGNORE);
		if (rc) {
			delete[] mpi_buffer[0];
		}
	} else {
		rc = true;
	}
	if (rc) {
		inc_instruction_pointer();
	}
}

void MultiGrid::phi_amr_boundary_communicate(int) {
	if (amr_cnt > 0) {
		int cnt, tag_send;
		for (int i = 0; i < amr_cnt; i++) {
			int sz = (PNX - 2) * (PNX - 2);
			mpi_buffer[i] = new Real[sz];
			cnt = 0;
			for (int l = amr_lb[i][2]; l <= amr_ub[i][2]; l++) {
				for (int k = amr_lb[i][1]; k <= amr_ub[i][1]; k++) {
					for (int j = amr_lb[i][0]; j <= amr_ub[i][0]; j++) {
						mpi_buffer[i][cnt] = phi.interp(j, k, l);
						++cnt;
					}
				}
			}
			tag_send = tag_gen(TAG_PHI, amr_id[i], amr_face[i]);
			MPI_Isend(mpi_buffer[i], cnt, MPI_DOUBLE_PRECISION, amr_child_proc[i], tag_send, MPI_COMM_WORLD, amr_request + i);
		}
	}
	inc_instruction_pointer();
}

void MultiGrid::phi_amr_boundary_wait(int) {
	int flag;
	bool rc;
	MPI_Testall(amr_cnt, amr_request, &flag, MPI_STATUS_IGNORE);
	rc = flag;
	if (rc) {
		for (int i = 0; i < amr_cnt; i++) {
			delete[] mpi_buffer[i];
		}
	}
	if (rc) {
		inc_instruction_pointer();
	}
}


void MultiGrid::phi_real_boundary_communicate(int) {
	const int dir = cx;
	int tag_send, tag_recv;
	MultiGrid* sibling;
	MPI_Datatype send_type, recv_type;
	for (int face = 2 * dir; face < 2 * dir + 2; face++) {
		sibling = dynamic_cast<MultiGrid*>(get_sibling(face));
		tag_recv = tag_gen(TAG_PHI, get_id(), face);
		if (is_real_bound(face)) {
			tag_send = tag_gen(TAG_PHI, sibling->get_id(), (face ^ 1));
			if (get_level() == OctNode::get_max_level()) {
				send_type = MPI_send_bnd1_t[face];
				recv_type = MPI_recv_bnd1_t[face];
			} else {
				send_type = MPI_send_bnd2_t[face];
				recv_type = MPI_recv_bnd2_t[face];
			}
			MPI_Isend(phi.ptr(), 1, send_type, sibling->proc(), tag_send, MPI_COMM_WORLD, send_request + face);
			MPI_Irecv(phi.ptr(), 1, recv_type, guard_proc(face), tag_recv, MPI_COMM_WORLD, recv_request + face);
		} else if (!is_phys_bound(face)) {
			MPI_Irecv(phi.ptr(), 1, MPI_recv_bnd1_t[face], guard_proc(face), tag_recv, MPI_COMM_WORLD, recv_request + face);
			send_request[face] = MPI_REQUEST_NULL;
		} else {
			send_request[face] = MPI_REQUEST_NULL;
			recv_request[face] = MPI_REQUEST_NULL;
		}
	}
	inc_instruction_pointer();
}

void MultiGrid::phi_real_boundary_wait(int) {
	const int dir = cx;
	bool rc;
	int flag_recv, flag_send;
	MPI_Testall(2, recv_request + 2 * dir, &flag_recv, MPI_STATUS_IGNORE);
	MPI_Testall(2, send_request + 2 * dir, &flag_send, MPI_STATUS_IGNORE);
	rc = flag_recv && flag_send;
	if (rc) {
		inc_instruction_pointer();
	}
}

void MultiGrid::phi_real_boundary_begin_loop(int) {
	cx = 0;
	inc_instruction_pointer();
	ax = get_instruction_pointer();
}

void MultiGrid::phi_real_boundary_end_loop(int) {
	cx++;
	if (cx < 3) {
		set_instruction_pointer(ax);
	} else {
		inc_instruction_pointer();
	}
}

void MultiGrid::init() {
	OctNode::init();
	dx = h0;
}

void MultiGrid::set_phi(int i, int j, int k, Real a) {
	phi(i, j, k) = a;
}

void MultiGrid::set_refine_flags() {
	ChildIndex c;
	if (get_level() < 1) {
		for (int i = 0; i < OCT_NCHILD; i++) {
			set_refine_flag(i, true);
		}
	} else if (get_level() < OctNode::get_max_level_allowed()) {
		for (int k = 1; k < PNX - 1; k++) {
			c.set_z(2 * k / PNX);
			for (int j = 1; j < PNX - 1; j++) {
				c.set_y(2 * j / PNX);
				for (int i = 1; i < PNX - 1; i++) {
					c.set_x(2 * i / PNX);
					if (!get_refine_flag(c)) {
						if (S(i, j, k) > 1.0e-3) {
							set_refine_flag(c, true);
						}
					}
				}
			}
		}
	}
}

Real MultiGrid::get_output_point(int i, int j, int k, int l) const {
	if (l == 0) {
		return S(i, j, k);
	} else {
		return phi(i, j, k);
	}
}

int MultiGrid::nvar_output() const {
	return 2;
}

const char* MultiGrid::output_field_names(int i) const {
	if (i == 0) {
		return "rho";
	} else {
		return "phi";
	}
}

void MultiGrid::deallocate_arrays() {
	fx.deallocate();
	fy.deallocate();
	fz.deallocate();
//	phi_fx.deallocate();
//	phi_fy.deallocate();
//	phi_fz.deallocate();
	S.deallocate();
	phi0.deallocate();
	phi.deallocate();
	dphi.deallocate();
	dphi1.deallocate();
}

void MultiGrid::allocate_arrays() {
	fx.allocate();
	fy.allocate();
	fz.allocate();
//	phi_fx.allocate();
//	phi_fy.allocate();
//	phi_fz.allocate();
	S.allocate();
	phi0.allocate();
	phi.allocate();
	dphi.allocate();
	dphi1.allocate();
	int i, j, k;
	for (k = 0; k < PNX; k++) {
		for (j = 0; j < PNX; j++) {
			for (i = 0; i < PNX; i++) {
				phi(i, j, k) = 0.0;
				dphi1(i, j, k) = 0.0;
				dphi(i, j, k) = 0.0;
				phi0(i, j, k) = 0.0;
			}
		}
	}
}

void MultiGrid::inject_from_parent(ChildIndex c) {
	if (proc() == MPI_rank()) {
		const MultiGrid* p = dynamic_cast<const MultiGrid*>(get_parent());
		Real u;
		int k, j, k0, j0, i, i0;
		for (k = 1; k < PNX - 1; k += 2) {
			for (j = 1; j < PNX - 1; j += 2) {
				k0 = (1 + k) / 2 + c.get_z() * (PNX / 2 - 1);
				j0 = (1 + j) / 2 + c.get_y() * (PNX / 2 - 1);
				for (i = 1, i0 = 1 + c.get_x() * (PNX / 2 - 1); i < PNX - 1; i += 2, i0++) {
					u = p->phi(i0, j0, k0);
					phi(i + 0, j + 0, k + 0) = u;
					phi(i + 1, j + 0, k + 0) = u;
					phi(i + 0, j + 1, k + 0) = u;
					phi(i + 1, j + 1, k + 0) = u;
					phi(i + 0, j + 0, k + 1) = u;
					phi(i + 1, j + 0, k + 1) = u;
					phi(i + 0, j + 1, k + 1) = u;
					phi(i + 1, j + 1, k + 1) = u;
					phi0(i + 0, j + 0, k + 0) = u;
					phi0(i + 1, j + 0, k + 0) = u;
					phi0(i + 0, j + 1, k + 0) = u;
					phi0(i + 1, j + 1, k + 0) = u;
					phi0(i + 0, j + 0, k + 1) = u;
					phi0(i + 1, j + 0, k + 1) = u;
					phi0(i + 0, j + 1, k + 1) = u;
					phi0(i + 1, j + 1, k + 1) = u;
				}
			}
		}
	}
}

MultiGrid::MultiGrid() :
		OctNode() {
	static bool init = false;
	if (!init) {
		origin[0] = origin[1] = origin[2] = 0.0;
		init = true;
	}
	amr_cnt = 0;
	dx = h0;
	offset = 0;

}

MultiGrid::~MultiGrid() {
}

bool MultiGrid::poisson_zone_is_refined(int i, int j, int k) const {
	ChildIndex c;
	i = (2 * i) / PNX;
	j = (2 * j) / PNX;
	k = (2 * k) / PNX;
	c.set_index(i, j, k);
	return (get_child(c) != NULL);
}

Real MultiGrid::xc(int i) const {
	return MultiGrid::xf(i) + 0.5 * dx;
}

Real MultiGrid::yc(int j) const {
	return MultiGrid::yf(j) + 0.5 * dx;
}

Real MultiGrid::zc(int k) const {
	return MultiGrid::zf(k) + 0.5 * dx;
}


Real MultiGrid::xf(int i) const {
	return Real(offset[0] + i) * dx - (PNX / 2) * h0 - origin[0];
}

Real MultiGrid::yf(int i) const {
	return Real(offset[1] + i) * dx - (PNX / 2) * h0 - origin[1];
}

Real MultiGrid::zf(int i) const {
	return Real(offset[2] + i) * dx - (PNX / 2) * h0 - origin[2];
}

void MultiGrid::relax_compute(int) {
	for (int k = 1; k < PNX - 1; k++) {
		for (int j = 1; j < PNX - 1; j++) {
			for (int i = (j + k) % 2 + 1; i < PNX - 1; i += 2) {
				dphi(i, j, k) += dphi1(i, j, k) + dphi.divergence(i, j, k) / 6.0;
			}
		}
	}
	for (int k = 1; k < PNX - 1; k++) {
		for (int j = 1; j < PNX - 1; j++) {
			for (int i = (j + k + 1) % 2 + 1; i < PNX - 1; i += 2) {
				dphi(i, j, k) += dphi1(i, j, k) + dphi.divergence(i, j, k) / 6.0;
			}
		}
	}
	inc_instruction_pointer();
}

void MultiGrid::relax_begin_down_loop(int) {
	cx = 1;
	inc_instruction_pointer();
	ax = get_instruction_pointer();
}

void MultiGrid::relax_end_down_loop(int) {
	cx++;
	if (cx <= 2) {
		set_instruction_pointer(ax);
	} else {
		inc_instruction_pointer();
	}
}
void MultiGrid::relax_begin_up_loop(int) {
	cx = 0;
	inc_instruction_pointer();
	ax = get_instruction_pointer();
}

void MultiGrid::relax_end_up_loop(int) {
	cx++;
	if (cx < 4) {
		set_instruction_pointer(ax);
	} else {
		inc_instruction_pointer();
	}
}

void MultiGrid::relax_communicate(int) {
	int face, tag_send, tag_recv;
	OctNode* sibling;
	static Array3d<Real, PNX, PNX, PNX> tmp0(true);
	for (face = 0; face < 6; face++) {
		sibling = get_sibling(face);
		recv_request[face] = MPI_REQUEST_NULL;
		send_request[face] = MPI_REQUEST_NULL;
		if (!is_phys_bound(face)) {
			if (is_real_bound(face)) {
				tag_send = tag_gen(TAG_RELAX, sibling->get_id(), (face ^ 1), cx);
				MPI_Isend(dphi.ptr(), 1, MPI_send_bnd1_t[face], sibling->proc(), tag_send, MPI_COMM_WORLD, send_request + face);
			}
			if (is_real_bound(face) || cx == 0) {
				tag_recv = tag_gen(TAG_RELAX, get_id(), face, cx);
				MPI_Irecv(dphi.ptr(), 1, MPI_recv_bnd1_t[face], guard_proc(face), tag_recv, MPI_COMM_WORLD, recv_request + face);
			}
		}
	}
	inc_instruction_pointer();
}

void MultiGrid::relax_wait(int) {
	int flag_recv, flag_send;
	bool ready;
	MPI_Testall(6, recv_request, &flag_recv, MPI_STATUS_IGNORE);
	MPI_Testall(6, send_request, &flag_send, MPI_STATUS_IGNORE);
	ready = flag_recv && flag_send;
	if (ready) {
		inc_instruction_pointer();
	}
}

void MultiGrid::residual_error_compute(int) {
	Real sum;
	sum = 0.0;
	const Real h3 = pow(dx, 3);
	const Real h1 = dx;
	int i, j, k;
	Real df;
	const Real h1inv = 1.0 / h1;
	for (k = 1; k < PNX - 1; k++) {
		for (j = 1; j < PNX - 1; j++) {
			for (i = 1; i < PNX - 1; i++) {
				if (!poisson_zone_is_refined(i, j, k)) {
					df = -(fx(i + 1, j, k) - fx(i, j, k)) * h1inv;
					df -= (fy(i, j + 1, k) - fy(i, j, k)) * h1inv;
					df -= (fz(i, j, k + 1) - fz(i, j, k)) * h1inv;
					dphi(i, j, k) = (df - S(i, j, k));
					sum += fabs(dphi(i, j, k)) * h3;
				} else {
					dphi(i, j, k) = 0.0;
				}
			}
		}
	}
	inc_instruction_pointer();
//	printf( "\t \t \t%e\n", sum);
	st0 = sum;
}

void MultiGrid::null(int) {
}

void MultiGrid::vdown_init_recv(int) {
	int tag;
	MultiGrid* child;
	for (int j = 0; j < PNX; j++) {
		for (int k = 0; k < PNX; k++) {
			dphi(0, j, k) = 0.0;
			dphi(j, 0, k) = 0.0;
			dphi(j, k, 0) = 0.0;
			dphi(PNX - 1, j, k) = 0.0;
			dphi(j, PNX - 1, k) = 0.0;
			dphi(j, k, PNX - 1) = 0.0;
		}
	}
	for (int i = 0; i < 8; i++) {
		child = dynamic_cast<MultiGrid*>(get_child(i));
		if (child != NULL) {
			tag = tag_gen(TAG_VDOWN, get_id(), i);
			MPI_Irecv(dphi.ptr(), 1, MPI_comm_child_t[i], get_child(i)->proc(), tag, MPI_COMM_WORLD, recv_request + i);
		} else {
			recv_request[i] = MPI_REQUEST_NULL;
		}
	}
	inc_instruction_pointer();
}

void MultiGrid::vdown_init_recv_wait(int) {
	int flag;
	MPI_Testall(8, recv_request, &flag, MPI_STATUS_IGNORE);
	if (flag) {
		inc_instruction_pointer();
	}
}

void MultiGrid::vdown_init_compute(int) {
	const Real dx2 = dx * dx;
	for (int i = 1; i < PNX - 1; i++) {
		for (int j = 1; j < PNX - 1; j++) {
			for (int k = 1; k < PNX - 1; k++) {
				if (!poisson_zone_is_refined(i, j, k)) {
					dphi(i, j, k) = (phi.divergence(i, j, k) - dx2 * S(i, j, k)) * (1.0 / 6.0);
				}
			}
		}
	}
	inc_instruction_pointer();
}

void MultiGrid::vdown_init_adjust_send(int) {
	Vector<int, 3> lb, ub;
	Real v;
	int i, cnt, tag, dir;
	for (int face = 0; face < 6; face++) {
		dir = face / 2;
		if (is_amr_bound(face)) {
			mpi_buffer[face] = new Real[PNX * PNX / 4];
			i = 1 + (PNX - 2) * (face & 1);
			cnt = 0;
			for (int k = 1; k < PNX - 1; k += 2) {
				for (int j = 1; j < PNX - 1; j += 2) {
					if (dir == 0) {
						v = +phi(i, j + 0, k + 0) - phi(i - 1, j + 0, k + 0);
						v += phi(i, j + 1, k + 0) - phi(i - 1, j + 1, k + 0);
						v += phi(i, j + 0, k + 1) - phi(i - 1, j + 0, k + 1);
						v += phi(i, j + 1, k + 1) - phi(i - 1, j + 1, k + 1);
					} else if (dir == 1) {
						v = +phi(j + 0, i, k + 0) - phi(j + 0, i - 1, k + 0);
						v += phi(j + 1, i, k + 0) - phi(j + 1, i - 1, k + 0);
						v += phi(j + 0, i, k + 1) - phi(j + 0, i - 1, k + 1);
						v += phi(j + 1, i, k + 1) - phi(j + 1, i - 1, k + 1);
					} else {
						v = +phi(j + 0, k + 0, i) - phi(j + 0, k + 0, i - 1);
						v += phi(j + 1, k + 0, i) - phi(j + 1, k + 0, i - 1);
						v += phi(j + 0, k + 1, i) - phi(j + 0, k + 1, i - 1);
						v += phi(j + 1, k + 1, i) - phi(j + 1, k + 1, i - 1);
					}
					v *= 0.5;
					mpi_buffer[face][cnt] = v;
					cnt++;
				}
			}
			assert(cnt == (PNX-2)*(PNX-2)/4);
			tag = tag_gen(TAG_DPHI, get_id(), face);
			MPI_Isend(mpi_buffer[face], cnt, MPI_DOUBLE_PRECISION, guard_proc(face), tag, MPI_COMM_WORLD, send_request + face);
		} else {
			send_request[face] = MPI_REQUEST_NULL;
		}
	}
	inc_instruction_pointer();
}

void MultiGrid::vdown_init_adjust_send_wait(int) {
	int flag;
	MPI_Testall(6, send_request, &flag, MPI_STATUS_IGNORE);
	if (flag) {
		for (int i = 0; i < 6; i++) {
			if (is_amr_bound(i)) {
				delete[] mpi_buffer[i];
			}
		}
		inc_instruction_pointer();
	}
}

void MultiGrid::vdown_init_adjust_recv(int) {
	if (amr_cnt > 0) {
		int tag;
		for (int i = 0; i < amr_cnt; i++) {
			const int sz = (PNX - 2) * (PNX - 2) / 4;
			mpi_buffer[i] = new Real[sz];
			tag = tag_gen(TAG_DPHI, amr_id[i], amr_face[i]);
			MPI_Irecv(mpi_buffer[i], sz, MPI_DOUBLE_PRECISION, amr_child_proc[i], tag, MPI_COMM_WORLD, amr_request + i);
		}
	}
	inc_instruction_pointer();
}

void MultiGrid::vdown_init_adjust_recv_wait(int) {
	int flag;
	if (amr_cnt > 0) {
		MPI_Testall(amr_cnt, amr_request, &flag, MPI_STATUS_IGNORE);
		if (flag) {
			int cnt, dir;
			Real v, d;
			int i;
			for (int ci = 0; ci < amr_cnt; ci++) {
				dir = amr_face[ci] / 2;
				cnt = 0;
				i = amr_lb[ci][dir] / 2;
				if (i == PNX / 2 - 1) {
					i++;
				} else if (i == PNX - 2) {
					i = PNX - 1;
				}
				if (dir == 0) {
					for (int k = amr_lb[ci][2] / 2; k <= amr_ub[ci][2] / 2; k++) {
						for (int j = amr_lb[ci][1] / 2; j <= amr_ub[ci][1] / 2; j++) {
							v = mpi_buffer[ci][cnt];
							d = (phi(i, j, k) - phi(i - 1, j, k) - v) * (1.0 / 6.0);
							if (!poisson_zone_is_refined(i, j, k)) {
								dphi(i, j, k) += d;
							}
							if (!poisson_zone_is_refined(i - 1, j, k)) {
								dphi(i - 1, j, k) -= d;
							}
							++cnt;
						}
					}
				} else if (dir == 1) {
					for (int k = amr_lb[ci][2] / 2; k <= amr_ub[ci][2] / 2; k++) {
						for (int j = amr_lb[ci][0] / 2; j <= amr_ub[ci][0] / 2; j++) {
							v = mpi_buffer[ci][cnt];
							d = (phi(j, i, k) - phi(j, i - 1, k) - v) * (1.0 / 6.0);
							if (!poisson_zone_is_refined(j, i, k)) {
								dphi(j, i, k) += d;
							}
							if (!poisson_zone_is_refined(j, i - 1, k)) {
								dphi(j, i - 1, k) -= d;
							}
							++cnt;
						}
					}
				} else {
					for (int k = amr_lb[ci][1] / 2; k <= amr_ub[ci][1] / 2; k++) {
						for (int j = amr_lb[ci][0] / 2; j <= amr_ub[ci][0] / 2; j++) {
							v = mpi_buffer[ci][cnt];
							d = (phi(j, k, i) - phi(j, k, i - 1) - v) * (1.0 / 6.0);
							if (!poisson_zone_is_refined(j, k, i)) {
								dphi(j, k, i) += d;
							}
							if (!poisson_zone_is_refined(j, k, i - 1)) {
								dphi(j, k, i - 1) -= d;
							}
							++cnt;
						}
					}
				}
				delete[] mpi_buffer[ci];
			}
			inc_instruction_pointer();
		}
	} else {
		inc_instruction_pointer();
	}
}

void MultiGrid::vdown_init_store_dphi(int) {
	for (int i = 1; i < PNX - 1; i++) {
		for (int j = 1; j < PNX - 1; j++) {
			for (int k = 1; k < PNX - 1; k++) {
				dphi1(i, j, k) = dphi(i, j, k);
			}
		}
	}
	inc_instruction_pointer();
}

void MultiGrid::vdown_init_send(int) {
	if (get_level() != 0) {
		int tag;
		int cnt = 0;
		mpi_buffer[0] = new Real[(PNX - 2) * (PNX - 2) * (PNX - 2)];
		for (int k = 1; k < PNX - 1; k += 2) {
			for (int j = 1; j < PNX - 1; j += 2) {
				for (int i = 1; i < PNX - 1; i += 2) {
					mpi_buffer[0][cnt] = dphi.oct_avg(i, j, k);
					cnt++;
				}
			}

		}
		tag = tag_gen(TAG_VDOWN, get_parent()->get_id(), my_child_index());
		if (cnt != 0) {
			//	printf("%i %i %i %i \n", get_parent()->proc(), get_location()[0],get_location()[1],get_location()[2]);
			MPI_Isend(mpi_buffer[0], cnt, MPI_DOUBLE_PRECISION, get_parent()->proc(), tag, MPI_COMM_WORLD, send_request);
		} else {
			send_request[0] = MPI_REQUEST_NULL;
		}
	}
	inc_instruction_pointer();
}

void MultiGrid::vdown_init_send_wait(int) {
	int flag;
	if (get_level() != 0) {
		MPI_Test(send_request, &flag, MPI_STATUS_IGNORE);
		if (flag) {
			delete[] mpi_buffer[0];
			inc_instruction_pointer();
		}
	} else {
		inc_instruction_pointer();
	}
}

void MultiGrid::retire_dphi(int) {
	if (proc() == MPI_rank()) {
		int i, j, k;
		for (k = 1; k < PNX - 1; k++) {
			for (j = 1; j < PNX - 1; j++) {
				for (i = 1; i < PNX - 1; i++) {
					phi(i, j, k) += dphi(i, j, k);
				}
			}
		}
	}
	inc_instruction_pointer();
}

void MultiGrid::vup_init_recv(int) {
	if (get_level() != 0) {
		const ChildIndex ci = my_child_index();
		const MultiGrid* p = dynamic_cast<const MultiGrid*>(get_parent());
		int i, j, k, sz, tag;
		Real a;
		MPI_Type_size(MPI_comm_child_t[ci], &sz);
		mpi_buffer[0] = new Real[sz / sizeof(Real)];
//	printf("%i\n", sz / (int) sizeof(double));
		tag = tag_gen(TAG_COARSE, get_id(), ci);
		MPI_Irecv(mpi_buffer[0], sz / sizeof(Real), MPI_DOUBLE_PRECISION, p->proc(), tag, MPI_COMM_WORLD, recv_request + 0);
		for (k = 1; k < PNX - 1; k += 2) {
			for (j = 1; j < PNX - 1; j += 2) {
				for (i = 1; i < PNX - 1; i += 2) {
					a = dphi.oct_avg(i, j, k);
					dphi(i + 0, j + 0, k + 0) -= a;
					dphi(i + 1, j + 0, k + 0) -= a;
					dphi(i + 0, j + 1, k + 0) -= a;
					dphi(i + 1, j + 1, k + 0) -= a;
					dphi(i + 0, j + 0, k + 1) -= a;
					dphi(i + 1, j + 0, k + 1) -= a;
					dphi(i + 0, j + 1, k + 1) -= a;
					dphi(i + 1, j + 1, k + 1) -= a;
				}
			}
		}
	}
	inc_instruction_pointer();
}

void MultiGrid::vup_init_recv_wait(int) {
	bool rc;
	if (get_level() != 0) {
		int flag;
		MPI_Test(recv_request, &flag, MPI_STATUS_IGNORE);
		rc = flag;
		if (rc) {
			int cnt = 0;
			for (int k = 1; k <= PNX - 2; k += 2) {
				for (int j = 1; j <= PNX - 2; j += 2) {
					for (int i = 1; i <= PNX - 2; i += 2) {
						dphi(i, j, k) += mpi_buffer[0][cnt];
						dphi(i, j, k + 1) += mpi_buffer[0][cnt];
						dphi(i, j + 1, k) += mpi_buffer[0][cnt];
						dphi(i, j + 1, k + 1) += mpi_buffer[0][cnt];
						dphi(i + 1, j, k) += mpi_buffer[0][cnt];
						dphi(i + 1, j, k + 1) += mpi_buffer[0][cnt];
						dphi(i + 1, j + 1, k) += mpi_buffer[0][cnt];
						dphi(i + 1, j + 1, k + 1) += mpi_buffer[0][cnt];
						cnt++;
					}
				}
			}
			int sz;
			MPI_Type_size(MPI_comm_child_t[my_child_index()], &sz);
			delete[] mpi_buffer[0];
		}
	} else {
		rc = true;
	}
	if (rc) {
		inc_instruction_pointer();
	}
}

void MultiGrid::vup_init_send(int) {
	MultiGrid* child;
	int tag;
	for (int ci = 0; ci < 8; ci++) {
		child = dynamic_cast<MultiGrid*>(get_child(ci));
		if (child == NULL) {
			send_request[ci] = MPI_REQUEST_NULL;
		} else {
			tag = tag_gen(TAG_COARSE, child->get_id(), ci);
			MPI_Isend(dphi.ptr(), 1, MPI_comm_child_t[ci], child->proc(), tag, MPI_COMM_WORLD, send_request + ci);
		}
	}
	inc_instruction_pointer();
}

void MultiGrid::vup_init_send_wait(int) {
	int flag;
	MPI_Testall(8, send_request, &flag, MPI_STATUS_IGNORE);
	if (flag) {
		inc_instruction_pointer();
	}

}

void MultiGrid::accumulate_com() {
    const int DX = (PNX / 2 - 1);
    const Real dv = pow(get_dx(), 3);
    Real dm;
    for (ChildIndex ci = 0; ci < 8; ci++) {
        if (get_child(ci) == NULL) {
            for (int k = 1 + ci.get_z() * DX; k < PNX / 2 + ci.get_z() * DX; k++) {
                for (int j = 1 + ci.get_y() * DX; j < PNX / 2 + ci.get_y() * DX; j++) {
                    for (int i = 1 + ci.get_x() * DX; i < PNX / 2 + ci.get_x() * DX; i++) {
                        dm = get_source(i, j, k) * dv;
                        com[0] += dm;
                        com[1] += MultiGrid::xc(i)*dm;
                        com[2] += MultiGrid::yc(j)*dm;
                        com[3] += MultiGrid::zc(k)*dm;
                    }
                }
            }
        }
    }
}

Real MultiGrid::compute_phi(Real x, Real y, Real z) {
    static AssociatedLegendrePolynomial P(LMAX);
    Real real = 0.0;
    Real r, theta, phi, prpow, mphi, rpow;
    x -= com[1];
    y -= com[2];
    z -= com[3];
    r = sqrt(x * x + y * y + z * z);
    theta = acos(z / r);
    phi = atan2(y, x);
    P.generate(cos(theta));

    rpow = 1.0;
    for (int l = 0; l <= LMAX; l++) {
        rpow /= r;
        for (int m = -l; m <= l; m++) {
            prpow = P.get(l, m) * rpow;
            mphi = Real(m) * phi;
            real -= prpow * (r_poles[l][m] * cos(mphi) - i_poles[l][m] * sin(mphi));
        }
    }
    return real;

}

void MultiGrid::compute_local_physical_boundaries() {
    Vector<int, 3> ub, lb;
    for (int f = 0; f < OCT_NSIB; f++) {
        if (is_phys_bound(f)) {
            lb = 0;
            ub = PNX - 1;
            lb[f / 2] = ub[f / 2] = (f % 2) * (PNX - 1);
            for (Indexer3d i(lb, ub); !i.end(); i++) {
                set_phi(i[0], i[1], i[2], MultiGrid::compute_phi(MultiGrid::xc(i[0]), MultiGrid::yc(i[1]), MultiGrid::zc(i[2])));
            }
        }
    }
}

void MultiGrid::poles_compute() {
    static AssociatedLegendrePolynomial P(LMAX);
    int l, m, xlb, xub, ylb, yub, zlb, zub;
    Real theta, phi, r, x, y, z, prpow, mphi;
    Real dv = pow(get_dx(), 3) / (4.0 * M_PI);
    Real x2, z2, y2;
    for (ChildIndex ci = 0; ci < 8; ci++) {
        if (get_child(ci) == NULL) {
            xlb = 1 + ci.get_x() * (PNX / 2 - 1);
            ylb = 1 + ci.get_y() * (PNX / 2 - 1);
            zlb = 1 + ci.get_z() * (PNX / 2 - 1);
            xub = xlb + (PNX / 2 - 1) - 1;
            yub = ylb + (PNX / 2 - 1) - 1;
            zub = zlb + (PNX / 2 - 1) - 1;
            for (int j = ylb; j <= yub; j++) {
                y = MultiGrid::yc(j) - com[2];
                y2 = y * y;
                for (int i = xlb; i <= xub; i++) {
                    x = MultiGrid::xc(i) - com[1];
                    x2 = x * x;
                    phi = atan2(y, x);
                    for (int k = zlb; k <= zub; k++) {
                        z = MultiGrid::zc(k) - com[3];
                        z2 = z * z;
                        r = sqrt(x2 + y2 + z2);
                        theta = acos(z / r);
                        P.generate(cos(theta));
                        for (l = 0; l <= LMAX; l++) {
                            for (m = -l; m <= l; m++) {
                                prpow = P.get(l, m) * pow(r, l) * dv * get_source(i, j, k);
                                mphi = Real(m) * phi;
                                r_poles[l][m] += prpow * cos(mphi);
                                i_poles[l][m] -= prpow * sin(mphi);
                            }
                        }
                    }
                }
            }
        }
    }
}

