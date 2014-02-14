#ifndef INDEXER3D_H_
#define INDEXER3D_H_

#include "vector.h"

class Indexer3d: public Vector<int, 3> {
private:
	Vector<int, 3> lb, ub;
public:
	Indexer3d(Vector<int, 3> _lb, Vector<int, 3> _ub) {
		lb = _lb;
		ub = _ub;
		*static_cast<Vector<int, 3>*>(this) = lb;
	}
	void operator++() {
		operator++(0);
	}
	void operator++(int) {
		(*this)[0]++;
		if ((*this)[0] > ub[0]) {
			(*this)[0] = lb[0];
			(*this)[1]++;
			if ((*this)[1] > ub[1]) {
				(*this)[1] = lb[1];
				(*this)[2]++;
			}
		}
	}
	bool end() const {
		if ((*this)[2] > ub[2]) {
			return true;
		} else {
			return false;
		}
	}

};

#endif /* INDEXER3D_H_ */
