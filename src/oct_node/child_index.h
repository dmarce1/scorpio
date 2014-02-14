#ifndef CHILD_INDEX_H_
#define CHILD_INDEX_H_

#define OCT_NCHILD 8

#include "../real.h"
#include "../vector.h"
#include "oct_face.h"

class ChildIndex {
private:
	int index;
public:
	Vector<int, 3> vec() const {
		Vector<int, 3> v;
		v[0] = get_x();
		v[1] = get_y();
		v[2] = get_z();
		return v;
	}
	ChildIndex();
	ChildIndex(int, int, int);
	ChildIndex(const ChildIndex&);
	ChildIndex(int);
	ChildIndex& operator=(const ChildIndex&);
	ChildIndex& operator=(int);
	virtual ~ChildIndex();
	Vector<int, 3> vector() const;
	int get_x() const;
	int get_y() const;
	int get_z() const;
	operator int() const;
	ChildIndex operator++();
	ChildIndex operator++(int);
	void set_index(int, int, int);
	void flip_x();
	void flip_y();
	void flip_z();
	void set_x(int);
	void set_y(int);
	void set_z(int);
	OctFace x_face() const;
	OctFace y_face() const;
	OctFace z_face() const;
};

#endif /* CHILD_INDEX_H_ */
