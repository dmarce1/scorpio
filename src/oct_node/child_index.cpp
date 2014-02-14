#include "assert.h"
#include "child_index.h"

ChildIndex::ChildIndex() {
	index = 0;
	return;
}

OctFace ChildIndex::x_face() const {
	return 0x1 & index;
}

ChildIndex::ChildIndex(int a, int b, int c) {
	index = a | (b << 1) | (c << 2);
}

OctFace ChildIndex::y_face() const {
	return ((0x2 & index) >> 1) + 2;
}

OctFace ChildIndex::z_face() const {
	return ((0x4 & index) >> 2) + 4;
}

ChildIndex::~ChildIndex() {
	return;
}

ChildIndex ChildIndex::operator++() {
	index++;
	return *this;
}

Vector<int, 3> ChildIndex::vector() const {
	Vector<int, 3> v;
	v[0] = get_x();
	v[1] = get_y();
	v[2] = get_z();
	return v;
}

ChildIndex ChildIndex::operator++(int) {
	ChildIndex prev(*this);
	index++;
	return prev;
}

ChildIndex::ChildIndex(const ChildIndex& a) {
	index = a.index;
}

ChildIndex& ChildIndex::operator=(const ChildIndex& a) {
	index = a.index;
	return *this;
}

ChildIndex::ChildIndex(int a) {
	index = a;
}

ChildIndex& ChildIndex::operator=(int a) {
	index = a;
	return *this;
}

int ChildIndex::get_x() const {
	return index & 0x1;
}

int ChildIndex::get_y() const {
	return (index >> 1) & 0x1;
}

int ChildIndex::get_z() const {
	return (index >> 2) & 0x1;
}

ChildIndex::operator int() const {
	return index;
}

void ChildIndex::set_index(int a, int b, int c) {
	index = a + 2 * b + 4 * c;
}

void ChildIndex::flip_x() {
	index ^= 0x1;
}

void ChildIndex::flip_y() {
	index ^= 0x2;
}

void ChildIndex::flip_z() {
	index ^= 0x4;
}

void ChildIndex::set_x(int a) {
	assert( a == 0 || a == 1 );
	if (a == 0) {
		index &= ~(0x1);
	} else {
		index |= 0x1;
	}
}

void ChildIndex::set_y(int a) {
	assert( a == 0 || a == 1 );
	if (a == 0) {
		index &= ~(0x2);
	} else {
		index |= 0x2;
	}
}

void ChildIndex::set_z(int a) {
	assert( a == 0 || a == 1 );
	if (a == 0) {
		index &= ~(0x4);
	} else {
		index |= 0x4;
	}
}
