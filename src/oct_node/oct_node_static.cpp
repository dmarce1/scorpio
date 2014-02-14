#include <mpi.h>
#include "oct_node.h"

int OctNode::max_refine_level = 5;
Real* OctNode::output_buffer;
Real OctNode::grid_time = 0.0;
OctNode* OctNode::root;
int OctNode::node_sums[MAX_LEVEL] = { 0 };
int OctNode::node_total = 0;
int OctNode::node_counter;
int OctNode::maxlev = 0;
int OctNode::id = 0;
bool OctNode::initialized = false;
int OctNode::local_node_cnt;
OctNode** OctNode::local_node;


int OctNode::get_node_cnt() {
	return node_total;
}

int OctNode::get_node_cnt(int lev) {
	return node_sums[lev];
}

void OctNode::initialize_grids() {
	for (int i = 0; i < OctNode::get_local_node_cnt(); i++) {
		OctNode::get_local_node(i)->initialize();

	}
}

Real OctNode::get_time() {
	return grid_time;
}

void OctNode::set_time(Real a) {
	grid_time = a;
}

OctNode* OctNode::get_root() {
	return root;
}

int OctNode::get_max_level() {
	return maxlev;
}

OctNode* OctNode::get_local_node(int i) {
	return local_node[i];
}

OctNode** OctNode::get_local_node_ptr() {
	return local_node;
}

int OctNode::get_local_node_cnt() {
	return local_node_cnt;
}

bool OctNode::check_for_refine() {
	OctNode* g;
	bool rc = false;
	if (max_refine_level == 0) {
		return false;
	}
	root->clear_refine_flags();
	for (int i = 0; i < get_local_node_cnt(); i++) {
		g = (get_local_node(i));
		g->set_refine_flags();
		g->propagate_refine_flags_up();
	}
	for (int lev = get_max_level(); lev >= 0; lev--) {
		for (int i = 0; i < get_local_node_cnt(); i++) {
			g = (get_local_node(i));
			if (g->get_level() == lev) {
				g->enforce_proper_nesting();
				g->propagate_refine_flags_up();
			}
		}
		reduce_refine_flags();
	}
	rc = root->use_refine_flags();
	if (rc) {
		root->find_local_nodes();
	}
	return rc;
}

void OctNode::reduce_refine_flags() {
	const int N = OCT_NCHILD * (get_node_cnt() - get_node_cnt(max_refine_level));
	Bits bits(N);
	Bits all_bits(N);
	int counter = 0;
	root->set_refine_bits(&bits, &counter);
	MPI_Allreduce(bits.ptr(), all_bits.ptr(), bits.size(), MPI_BYTE, MPI_BOR, MPI_COMM_WORLD);
	counter = 0;
	root->use_refine_bits(&all_bits, &counter);
}
