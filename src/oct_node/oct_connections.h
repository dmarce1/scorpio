/*
 * oct_connections.h
 *
 *  Created on: Mar 20, 2013
 *      Author: dmarce1
 */

#ifndef OCT_CONNECTIONS_H_
#define OCT_CONNECTIONS_H_

#include "child_index.h"

class OctNode;

class OctConnections {
private:
	OctNode* children[OCT_NCHILD];
	OctNode* siblings[OCT_NSIB];
	OctNode* parent;
public:
	OctConnections();
	virtual ~OctConnections();
	void set_sibling(OctFace f, OctNode* n);
	const OctNode* get_sibling(int i) const;
	OctNode* get_sibling(int i);
	OctNode* get_parent() const;
	OctNode* get_child(const ChildIndex& i);
	const OctNode* get_child(const ChildIndex& i) const;
	void set_child(const ChildIndex& i, OctNode* n);
	void set_parent(OctNode* n);
};

#endif /* OCT_CONNECTIONS_H_ */
