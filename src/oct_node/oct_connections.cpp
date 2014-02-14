#include "oct_connections.h"
#include <stdlib.h>

OctConnections::OctConnections() {
	parent = NULL;
	for (int i = 0; i < OCT_NCHILD; i++) {
		children[i] = NULL;
	}
	for (int i = 0; i < OCT_NSIB; i++) {
		siblings[i] = NULL;
	}
}

OctConnections::~OctConnections() {
	return;
}

void OctConnections::set_sibling(OctFace f, OctNode* n) {
	siblings[f] = n;
}

const OctNode* OctConnections::get_sibling(int i) const {
	return siblings[i];
}

OctNode* OctConnections::get_sibling(int i) {
	return siblings[i];
}

OctNode* OctConnections::get_parent() const {
	return parent;
}

OctNode* OctConnections::get_child(const ChildIndex& i) {
	return children[i];
}

const OctNode* OctConnections::get_child(const ChildIndex& i) const {
	return children[i];
}

void OctConnections::set_child(const ChildIndex& i, OctNode* n) {
	children[i] = n;
}

void OctConnections::set_parent(OctNode* n) {
	parent = n;
}


