#include "../defs.h"
#include "oct_node.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "../comm.h"
#include "../real.h"
#include "grid_output.h"
#include <mpi.h>
#ifndef USING_MIC
#include <silo.h>
#endif

void OctNode::set_max_level_allowed(int l) {
	max_refine_level = l;
}

int OctNode::get_max_level_allowed() {
	return max_refine_level;
}

bool OctNode::get_refine_flag(const ChildIndex& i) const {
	return refine[i];
}

void OctNode::set_refine_flag(const ChildIndex& i, bool f) {
	refine[i] = f;
}

int OctNode::next_proc() const {
	return next_processor;
}

int OctNode::last_proc() const {
	return last_processor;
}

int OctNode::proc() const {
	return processor;
}

void OctNode::set_to_next_proc() {
	last_processor = processor;
	processor = next_processor;
}

void OctNode::set_proc(int p) {
	processor = p;
}

int OctNode::get_id() const {
	return myid;
}

void OctNode::allocate_nodes() {
	bool flag;
	flag = proc() != next_proc();
	OctNode::set_to_next_proc();
	if (flag && last_processor != -1) {
		if (next_proc() == MPI_rank()) {
			this->allocate_arrays();
		} else {
			this->deallocate_arrays();
		}
	}
	for (int i = 0; i < OCT_NCHILD; i++) {
		if (get_child(i) != NULL) {
			get_child(i)->allocate_nodes();
		}
	}
}
void OctNode::init() {
	assert( get_level() == 0);
	root = this;
	set_proc(proc());
	if (proc() == MPI_rank()) {
		this->allocate_arrays();
	}
	for (int i = 0; i < OCT_NSIB; i++) {
		set_sibling(OctFace(i), NULL);
	}
	find_local_nodes();
}

void OctNode::find_local_nodes() {

	if (level == 0) {
		if (initialized) {
			delete[] local_node;
		}
		local_node = new OctNode*[get_node_cnt()];
		local_node_cnt = 0;
		initialized = true;
	}

	if (processor == MPI_rank()) {
		local_node[local_node_cnt] = this;
		local_node_cnt++;
	}

	for (int i = 0; i < 8; i++) {
		if (get_child(i)) {
			get_child(i)->find_local_nodes();
		}
	}

}

ChildIndex OctNode::my_child_index() const {
	assert( get_parent() != NULL);
	ChildIndex c;
	c.set_x(location[0] % 2);
	c.set_y(location[1] % 2);
	c.set_z(location[2] % 2);
	return c;
}

void OctNode::compute_distribution() {
	if (get_level() == 0) {
		node_counter = 0;
		//	printf( "%i %i\n", MPI_rank(), get_node_cnt() );
	}
//	next_processor = node_dist[get_level()] * MPI_size() / node_sums[get_level()];
//	node_dist[get_level()]++;
#ifdef RANK_ZERO_HAS_ONE_GRID
	if( node_counter == 0 ) {
		next_processor = 0;
	} else {
		next_processor = (node_counter-1) * (MPI_size()-1) / (get_node_cnt()-1)+1;
	}
#else
	next_processor = node_counter * MPI_size() / get_node_cnt();
#endif
	this->myid = node_counter;
	node_counter++;
	for (int i = 0; i < OCT_NCHILD; i++) {
		if (get_child(i) != NULL) {
			get_child(i)->compute_distribution();
		}
	}
	if (get_level() == 0) {
		id = node_counter;
	}
}

bool OctNode::is_real_bound(OctFace f) const {
	return get_sibling(f) != NULL;
}

void OctNode::write_checkpoint(FILE* fp) const {
	if (get_level() == 0) {
		fwrite(&grid_time, sizeof(Real), 1, fp);
	}
	bool child_flags[OCT_NCHILD];
	fwrite(&processor, sizeof(int), 1, fp);
	fwrite(&myid, sizeof(int), 1, fp);
	this->write(fp);
	for (int i = 0; i < OCT_NCHILD; i++) {
		child_flags[i] = (get_child(i) != NULL);
	}
	fwrite(child_flags, sizeof(bool), OCT_NCHILD, fp);
	for (int i = 0; i < OCT_NCHILD; i++) {
		if (get_child(i) != NULL) {
			get_child(i)->write_checkpoint(fp);
		}
	}
}

void OctNode::read_checkpoint(FILE* fp) {
	if (get_level() == 0) {
		fread(&grid_time, sizeof(Real), 1, fp);
	}
	bool child_flags[OCT_NCHILD];
	fread(&processor, sizeof(int), 1, fp);
	fread(&myid, sizeof(int), 1, fp);
	if (processor == MPI_rank()) {
		this->allocate_arrays();
	} else {
		this->deallocate_arrays();
	}
	this->read(fp);
	fread(child_flags, sizeof(bool), OCT_NCHILD, fp);
	for (int i = 0; i < OCT_NCHILD; i++) {
		if (child_flags[i]) {
			if (!get_child(i)) {
				this->create_child(i);
			}
		} else {
			if (get_child(i)) {
				this->destroy_child(i);
			}
		}
	}
	for (int i = 0; i < OCT_NCHILD; i++) {
		if (child_flags[i]) {
			get_child(i)->read_checkpoint(fp);
		}
	}
	if (get_level() == 0) {
		find_local_nodes();
	}
}

bool OctNode::is_amr_bound(OctFace f) const {
	return get_sibling(f) == NULL && !is_phys_bound(f);
}

OctNode::OctNode() :
		OctConnections() {

	for (int i = 0; i < 8; i++) {
		refine[i] = false;
	}
	processor = 0;
	next_processor = 0;
	last_processor = -1;
	level = 0;

	location = 0;
	if (node_sums[0] == 0) {
		node_sums[0] = 1;
		node_total = 1;
		myid = 0;
	}
	myid = id;
	id++;
}

bool OctNode::is_phys_bound(OctFace f) const {
	bool rc;
	switch (f) {
	case XL:
		rc = location[0] == 0;
		break;
	case XU:
		rc = location[0] == ((1 << level) - 1);
		break;
	case YL:
		rc = location[1] == 0;
		break;
	case YU:
		rc = location[1] == ((1 << level) - 1);
		break;
	case ZL:
		rc = location[2] == 0;
		break;
	default:
		//case ZU:
		rc = location[2] == ((1 << level) - 1);
		break;
	}
	return rc;
}

int OctNode::get_level() const {
	return level;
}

void OctNode::create_child(const ChildIndex& c) {
	assert( get_child(c) == NULL);
	ChildIndex s;
	set_child(c, this->new_octnode());
	get_child(c)->level = level + 1;
	//printf( "%i\n", level);
	get_child(c)->processor = processor;
	get_child(c)->set_parent(this);
	get_child(c)->location = location * 2 + c.vector();
	OctFace f1;
	OctFace f2;
	s = c;
	if (c.get_x() == 0) {
		s.set_x(1);
		f1 = XL;
		f2 = XU;
	} else {
		s.set_x(0);
		f2 = XL;
		f1 = XU;
	}
	if (get_child(s) != NULL) {
		get_child(c)->set_sibling(f2, get_child(s));
		get_child(s)->set_sibling(f1, get_child(c));
	}
	if (get_sibling(f1) != NULL) {
		get_child(c)->set_sibling(f1, get_sibling(f1)->get_child(s));
		if (get_sibling(f1)->get_child(s) != NULL) {
			(get_sibling(f1)->get_child(s))->set_sibling(f2, get_child(c));
		}
	}
	s = c;
	if (c.get_y() == 0) {
		s.set_y(1);
		f1 = YL;
		f2 = YU;
	} else {
		s.set_y(0);
		f2 = YL;
		f1 = YU;
	}
	if (get_child(s) != NULL) {
		get_child(c)->set_sibling(f2, get_child(s));
		get_child(s)->set_sibling(f1, get_child(c));
	}
	if (get_sibling(f1) != NULL) {
		get_child(c)->set_sibling(f1, get_sibling(f1)->get_child(s));
		if (get_sibling(f1)->get_child(s) != NULL) {
			(get_sibling(f1)->get_child(s))->set_sibling(f2, get_child(c));
		}
	}
	s = c;
	if (c.get_z() == 0) {
		s.set_z(1);
		f1 = ZL;
		f2 = ZU;
	} else {
		s.set_z(0);
		f2 = ZL;
		f1 = ZU;
	}
	if (get_child(s) != NULL) {
		get_child(c)->set_sibling(f2, get_child(s));
		get_child(s)->set_sibling(f1, get_child(c));
	}
	if (get_sibling(f1) != NULL) {
		get_child(c)->set_sibling(f1, get_sibling(f1)->get_child(s));
		if (get_sibling(f1)->get_child(s) != NULL) {
			(get_sibling(f1)->get_child(s))->set_sibling(f2, get_child(c));
		}
	}
	if (node_sums[get_level() + 1] == 0) {
		maxlev++;
	}
	node_total++;
	node_sums[get_level() + 1]++;
	get_child(c)->set_proc(OctNode::proc());
	if (proc() == MPI_rank()) {
		get_child(c)->allocate_arrays();
	}

}

void OctNode::destroy_child(const ChildIndex& c) {
	assert( get_child(c) != NULL);
	ChildIndex s;
	OctFace f1;
	OctFace f2;
	s = c;
	if (c.get_x() == 0) {
		s.set_x(1);
		f1 = XL;
		f2 = XU;
	} else {
		s.set_x(0);
		f2 = XL;
		f1 = XU;
	}
	if (get_child(s) != NULL) {
		get_child(s)->set_sibling(f1, NULL);
	}
	if (get_sibling(f1) != NULL) {
		if (get_sibling(f1)->get_child(s) != NULL) {
			(get_sibling(f1)->get_child(s))->set_sibling(f2, NULL);
		}
	}
	s = c;
	if (c.get_y() == 0) {
		s.set_y(1);
		f1 = YL;
		f2 = YU;
	} else {
		s.set_y(0);
		f2 = YL;
		f1 = YU;
	}
	if (get_child(s) != NULL) {
		get_child(s)->set_sibling(f1, NULL);
	}
	if (get_sibling(f1) != NULL) {
		if (get_sibling(f1)->get_child(s) != NULL) {
			(get_sibling(f1)->get_child(s))->set_sibling(f2, NULL);
		}
	}
	s = c;
	if (c.get_z() == 0) {
		s.set_z(1);
		f1 = ZL;
		f2 = ZU;
	} else {
		s.set_z(0);
		f2 = ZL;
		f1 = ZU;
	}
	if (get_child(s) != NULL) {
		get_child(s)->set_sibling(f1, NULL);
	}
	if (get_sibling(f1) != NULL) {
		if (get_sibling(f1)->get_child(s) != NULL) {
			(get_sibling(f1)->get_child(s))->set_sibling(f2, NULL);
		}
	}
	delete get_child(c);
	set_child(c, NULL);
	node_sums[get_level() + 1]--;
	if (node_sums[get_level() + 1] == 0) {
		maxlev--;
	}
	node_total--;
}

OctNode::~OctNode() {
	for (int i = 0; i < OCT_NCHILD; i++) {
		if (get_child(i) != NULL) {
			this->destroy_child(i);
		}
	}
}

void OctNode::clear_refine_flags() {
	for (int i = 0; i < OCT_NCHILD; i++) {
		refine[i] = false;
	}
	if (get_level() < max_refine_level) {
		for (int i = 0; i < OCT_NCHILD; i++) {
			if (this->get_child(i) != NULL) {
				(get_child(i))->clear_refine_flags();
			}
		}
	}
}

void OctNode::enforce_proper_nesting() {
	if (get_level() > 0) {
		const ChildIndex mi = my_child_index();
		OctNode* parent = get_parent();
		ChildIndex ci, i;
		OctNode* g;
		OctNode* tmp;
		for (int ci0 = 0; ci0 < OCT_NCHILD; ci0++) {
			ci = ci0;
			if (refine[ci]) {
				i = mi;
				i.flip_x();
				g = NULL;
				if ((mi.get_x() == ci.get_x()) && !is_phys_bound(ci.x_face())) {
					g = (parent->get_sibling(ci.x_face()));
					assert(g!=NULL);
				} else {
					g = (parent);
					assert(g!=NULL);
				}
				if (g != NULL) {
					(g)->refine[i] = true;
				}
				i = mi;
				i.flip_y();
				g = NULL;
				if ((mi.get_y() == ci.get_y()) && !is_phys_bound(ci.y_face())) {
					g = (parent->get_sibling(ci.y_face()));
					assert(g!=NULL);
				} else {
					g = (parent);
					assert(g!=NULL);
				}
				if (g != NULL) {
					(g)->refine[i] = true;
				}
				i = mi;
				i.flip_z();
				g = NULL;
				if ((mi.get_z() == ci.get_z()) && !is_phys_bound(ci.z_face())) {
					g = (parent->get_sibling(ci.z_face()));
					assert(g!=NULL);
				} else {
					g = (parent);
					assert(g!=NULL);
				}
				if (g != NULL) {
					(g)->refine[i] = true;
				}
				i = mi;
				i.flip_x();
				i.flip_y();
				g = NULL;
				if ((mi.get_x() == ci.get_x()) && !is_phys_bound(ci.x_face())) {
					if ((mi.get_y() == ci.get_y()) && !is_phys_bound(ci.y_face())) {
						tmp = parent->get_sibling(ci.x_face());
						assert(tmp);
						tmp = tmp->get_sibling(ci.y_face());
						g = (tmp);
						assert(g!=NULL);
					} else {
						g = (parent->get_sibling(ci.x_face()));
						assert(g!=NULL);
					}
				} else {
					if ((mi.get_y() == ci.get_y()) && !is_phys_bound(ci.y_face())) {
						g = (parent->get_sibling(ci.y_face()));
						assert(g!=NULL);
					} else {
						g = (parent);
						assert(g!=NULL);
					}
				}
				if (g != NULL) {
					(g)->refine[i] = true;
				}
				i = mi;
				i.flip_x();
				i.flip_z();
				g = NULL;
				if ((mi.get_x() == ci.get_x()) && !is_phys_bound(ci.x_face())) {
					if ((mi.get_z() == ci.get_z()) && !is_phys_bound(ci.z_face())) {
						tmp = parent->get_sibling(ci.x_face());
						tmp = tmp->get_sibling(ci.z_face());
						g = (tmp);
						assert(g!=NULL);
					} else {
						g = (parent->get_sibling(ci.x_face()));
						assert(g!=NULL);
					}
				} else {
					if ((mi.get_z() == ci.get_z()) && !is_phys_bound(ci.z_face())) {
						g = (parent->get_sibling(ci.z_face()));
						assert(g!=NULL);
					} else {
						g = (parent);
						assert(g!=NULL);
					}
				}
				if (g != NULL) {
					(g)->refine[i] = true;
				}
				i = mi;
				i.flip_y();
				i.flip_z();
				g = NULL;
				if ((mi.get_y() == ci.get_y()) && !is_phys_bound(ci.y_face())) {
					if ((mi.get_z() == ci.get_z()) && !is_phys_bound(ci.z_face())) {
						tmp = parent->get_sibling(ci.y_face());
						tmp = tmp->get_sibling(ci.z_face());
						g = (tmp);
						assert(g!=NULL);
					} else {
						g = (parent->get_sibling(ci.y_face()));
						assert(g!=NULL);
					}
				} else {
					if ((mi.get_z() == ci.get_z()) && !is_phys_bound(ci.z_face())) {
						g = (parent->get_sibling(ci.z_face()));
						assert(g!=NULL);
					} else {
						g = (parent);
						assert(g!=NULL);
					}
				}
				if (g != NULL) {
					(g)->refine[i] = true;
				}
				i = mi;
				i.flip_x();
				i.flip_y();
				i.flip_z();
				g = NULL;
				if ((mi.get_x() == ci.get_x()) && !is_phys_bound(ci.x_face())) {
					if ((mi.get_y() == ci.get_y()) && !is_phys_bound(ci.y_face())) {
						if ((mi.get_z() == ci.get_z()) && !is_phys_bound(ci.z_face())) {
							tmp = parent->get_sibling(ci.x_face());
							assert(tmp!=NULL);
							tmp = tmp->get_sibling(ci.y_face());
							assert(tmp!=NULL);
							tmp = tmp->get_sibling(ci.z_face());
							assert(tmp!=NULL);
							g = (tmp);
						} else {
							tmp = parent->get_sibling(ci.x_face());
							tmp = tmp->get_sibling(ci.y_face());
							g = (tmp);
							assert(g!=NULL);
						}
					} else {
						if ((mi.get_z() == ci.get_z()) && !is_phys_bound(ci.z_face())) {
							tmp = parent->get_sibling(ci.x_face());
							tmp = tmp->get_sibling(ci.z_face());
							g = (tmp);
							assert(g!=NULL);
						} else {
							tmp = parent->get_sibling(ci.x_face());
							g = (tmp);
							assert(g!=NULL);
						}
					}
				} else {
					if ((mi.get_y() == ci.get_y()) && !is_phys_bound(ci.y_face())) {
						if ((mi.get_z() == ci.get_z()) && !is_phys_bound(ci.z_face())) {
							tmp = parent->get_sibling(ci.y_face());
							tmp = tmp->get_sibling(ci.z_face());
							g = (tmp);
							assert(g!=NULL);
						} else {
							tmp = parent->get_sibling(ci.y_face());
							g = (tmp);
							assert(g!=NULL);
						}
					} else {
						if ((mi.get_z() == ci.get_z()) && !is_phys_bound(ci.z_face())) {
							tmp = parent->get_sibling(ci.z_face());
							g = (tmp);
							assert(g!=NULL);
						} else {
							g = (parent);
							assert(g!=NULL);
						}
					}
				}
				if (g != NULL) {
					(g)->refine[i] = true;
				}
			}
		}
	}
}

void OctNode::propagate_refine_flags_up() {
	OctNode *c, *p;
	c = this;
	for (int j = 0; j < 8; j++) {
		if (c->refine[j] == true || c->get_child(j) != NULL) {
			while (c->get_level() > 0) {
				p = (c->get_parent());
				p->refine[c->my_child_index()] = true;
				c = p;
			}
			break;
		}
	}
}

void OctNode::use_refine_bits(Bits* bits, int* n) {
	for (int i = 0; i < OCT_NCHILD; i++) {
		refine[i] = bits->get(*n);
		++(*n);
		if (get_child(i) != NULL && get_child(i)->get_level() != max_refine_level) {
			(get_child(i))->use_refine_bits(bits, n);
		}
	}
}

void OctNode::set_refine_bits(Bits* bits, int* n) {
	for (int i = 0; i < OCT_NCHILD; i++) {
		bits->set(*n, refine[i]);
		++(*n);
		if (get_child(i) != NULL && get_child(i)->get_level() != max_refine_level) {
			(get_child(i))->set_refine_bits(bits, n);
		}
	}
}

bool OctNode::use_refine_flags() {
	bool rc = false;
	for (int i = 0; i < OCT_NCHILD; i++) {
		if (refine[i] && (get_child(i) == NULL) && (get_level() + 1 <= max_refine_level)) {
			this->create_child(i);
			rc = true;
		} else if (!refine[i] && (get_child(i) != NULL)) {
			this->destroy_child(i);
			rc = true;
		}
	}
	for (int i = 0; i < OCT_NCHILD; i++) {
		if (this->get_child(i) != NULL) {
			bool test = (get_child(i))->use_refine_flags();
			rc = rc || test;
		}
	}
	return rc;
}

void OctNode::output(grid_output_t* ptr, int nx0, int bw0, Real dtheta) const {
	ChildIndex ci;
	const int nvar = this->nvar_output();
	int cnt;
	if (MPI_rank() == 0) {
		if (proc() != 0) {

			cnt = (nx0 - 2 * bw0) * (nx0 - 2 * bw0) * (nx0 - 2 * bw0);
			for (int i = 0; i < OCT_NCHILD; i++) {
				if (get_child(i) != NULL) {
					cnt -= (nx0 / 2 - bw0) * (nx0 / 2 - bw0) * (nx0 / 2 - bw0);
				}
			}
			if (cnt != 0) {
				MPI_Recv(output_buffer, cnt * nvar, MPI_DOUBLE_PRECISION, proc(), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
			}
		}
		for (int k = bw0; k < nx0 + 1 - bw0; k++) {
			for (int j = bw0; j < nx0 + 1 - bw0; j++) {
				for (int i = bw0; i < nx0 + 1 - bw0; i++) {
				    Real xtmp = this->xf(i);
					Real ytmp  = this->yf(j);
                   *(ptr->z) = this->zf(k);
                   *(ptr->x) = cos(dtheta)*xtmp - sin(dtheta)*ytmp;
                   *(ptr->y) = sin(dtheta)*xtmp + cos(dtheta)*ytmp;
					ptr->x++;
					ptr->y++;
					ptr->z++;
				}
			}
		}
		cnt = 0;
		for (int k = bw0; k < nx0 - bw0; k++) {
			ci.set_z(2 * k / nx0);
			for (int j = bw0; j < nx0 - bw0; j++) {
				ci.set_y(2 * j / nx0);
				for (int i = bw0; i < nx0 - bw0; i++) {
					ci.set_x(2 * i / nx0);
					if (get_child(ci) == NULL) {
						const int nxi = nx0 - 2 * bw0;
						const int i0 = i - bw0;
						const int j0 = j - bw0;
						const int k0 = k - bw0;
						(ptr->nodelist)[ptr->ni++] = (i0 + 0) + (nxi + 1) * ((j0 + 0) + (nxi + 1) * (k0 + 0)) + ptr->pi;
						(ptr->nodelist)[ptr->ni++] = (i0 + 1) + (nxi + 1) * ((j0 + 0) + (nxi + 1) * (k0 + 0)) + ptr->pi;
						(ptr->nodelist)[ptr->ni++] = (i0 + 1) + (nxi + 1) * ((j0 + 1) + (nxi + 1) * (k0 + 0)) + ptr->pi;
						(ptr->nodelist)[ptr->ni++] = (i0 + 0) + (nxi + 1) * ((j0 + 1) + (nxi + 1) * (k0 + 0)) + ptr->pi;
						(ptr->nodelist)[ptr->ni++] = (i0 + 0) + (nxi + 1) * ((j0 + 0) + (nxi + 1) * (k0 + 1)) + ptr->pi;
						(ptr->nodelist)[ptr->ni++] = (i0 + 1) + (nxi + 1) * ((j0 + 0) + (nxi + 1) * (k0 + 1)) + ptr->pi;
						(ptr->nodelist)[ptr->ni++] = (i0 + 1) + (nxi + 1) * ((j0 + 1) + (nxi + 1) * (k0 + 1)) + ptr->pi;
						(ptr->nodelist)[ptr->ni++] = (i0 + 0) + (nxi + 1) * ((j0 + 1) + (nxi + 1) * (k0 + 1)) + ptr->pi;
						if (proc() == 0) {
							for (int l = 0; l < nvar; l++) {
								ptr->ele[l][ptr->ei] = this->get_output_point(i, j, k, l);
							}
						} else {
							for (int l = 0; l < nvar; l++) {
								ptr->ele[l][ptr->ei] = output_buffer[nvar * cnt + l];
							}
							cnt++;
						}
						ptr->ei++;
					}
				}
			}
		}
		ptr->pi += (nx0 - 2 * bw0 + 1) * (nx0 - 2 * bw0 + 1) * (nx0 - 2 * bw0 + 1);
	} else {
		cnt = 0;
		if (proc() == MPI_rank()) {
			for (int k = bw0; k < nx0 - bw0; k++) {
				ci.set_z(2 * k / nx0);
				for (int j = bw0; j < nx0 - bw0; j++) {
					ci.set_y(2 * j / nx0);
					for (int i = bw0; i < nx0 - bw0; i++) {
						ci.set_x(2 * i / nx0);
						if (get_child(ci) == NULL) {
							for (int l = 0; l < nvar; l++) {
								output_buffer[nvar * cnt + l] = this->get_output_point(i, j, k, l);
							}
							cnt++;
						}
					}
				}
			}
			if (cnt != 0) {
				MPI_Send(output_buffer, cnt * nvar, MPI_DOUBLE_PRECISION, 0, 0, MPI_COMM_WORLD );
			}
		}
	}
	for (int i = 0; i < OCT_NCHILD; i++) {
		if (get_child(i) != NULL) {
			get_child(i)->output(ptr, nx0, bw0);
		}
	}
}

void OctNode::output(const char* prefix, int counter, int nx0, int bw0, double dtheta) const {
		output_buffer = new Real[this->nvar_output() * (nx0 - 2 * bw0) * (nx0 - 2 * bw0) * (nx0 - 2 * bw0)];
#ifdef USING_MIC
		printf( "Called output on MIC - no silo!\n");
		abort();
#else
	assert( get_level() == 0);
	char filename[32];
	int cnt, nzones, nnodes;
	grid_output_t go;
	float* coords[3];
	DBfile* db;
	DBoptlist* olist;
	char* coordnames[3];
	int shapesize[1];
	int shapecnt[1];
	int shapetype[1];
	int nshapes = 1;
	float ftime = float(get_time());
	const int nf = this->nvar_output();
	if (MPI_rank() == 0) {
		sprintf(filename, "%s.%i.silo", prefix, counter);
		olist = DBMakeOptlist(1);
		DBAddOption(olist, DBOPT_TIME, &ftime);
		db = DBCreate(filename, DB_CLOBBER, DB_LOCAL, "Euler Mesh", DB_PDB);
		cnt = get_node_cnt();
		nnodes = cnt * (nx0 - 2 * bw0 + 1) * (nx0 - 2 * bw0 + 1) * (nx0 - 2 * bw0 + 1);
		nzones = get_leaf_cnt() * (nx0 - 2 * bw0) * (nx0 - 2 * bw0) * (nx0 - 2 * bw0) / OCT_NCHILD;
		for (int i = 0; i < 3; i++) {
			coordnames[i] = new char[2];
			coords[i] = reinterpret_cast<float*>(new OReal[nnodes]);
		}
		go.nodelist = new int[8 * nzones];
		go.x = reinterpret_cast<OReal*>(coords[0]);
		go.y = reinterpret_cast<OReal*>(coords[1]);
		go.z = reinterpret_cast<OReal*>(coords[2]);
		go.ele = new OReal*[nf];
		for (int i = 0; i < nf; i++) {
			go.ele[i] = new OReal[nzones];
		}
		go.pi = 0;
		go.ni = 0;
		go.ei = 0;
		output(&go, nx0, bw0,dtheta);
		shapesize[0] = 8;
		shapecnt[0] = nzones;
		shapetype[0] = DB_ZONETYPE_HEX;
		strcpy(coordnames[0], "x");
		strcpy(coordnames[1], "y");
		strcpy(coordnames[2], "z");
		DBPutZonelist2(db, "zones", nzones, 3, go.nodelist, 8 * nzones, 0, 0, 0, shapetype, shapesize, shapecnt, nshapes, olist);
		DBPutUcdmesh(db, "mesh", 3, coordnames, coords, nnodes, nzones, "zones", NULL, DB_OREAL, olist);
		for (int i = 0; i < nf; i++) {
			DBPutUcdvar1(db, output_field_names(i), "mesh", (float*) go.ele[i], nzones, 0, 0, DB_OREAL, DB_ZONECENT, olist);
		}
		DBClose(db);
		char str[80];

		delete[] go.nodelist;
		for (int i = 0; i < nf; i++) {
			delete[] go.ele[i];
		}
		delete[] go.ele;
		for (int i = 0; i < 3; i++) {
			delete[] coords[i];
			delete[] coordnames[i];
		}
	} else {
		output(NULL, nx0, bw0);

	}
	delete[] output_buffer;
#endif
}

int OctNode::get_leaf_cnt() const {
	int leaf_cnt = 0;
	for (int i = 0; i < OCT_NCHILD; i++) {
		if (get_child(i) != NULL) {
			leaf_cnt += get_child(i)->get_leaf_cnt();
		} else {
			leaf_cnt++;
		}
	}
	return leaf_cnt;
}

