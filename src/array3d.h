#ifndef ARRAY3D_H_
#define ARRAY3D_H_

#include <mpi.h>
#include "real.h"
#include "vector.h"
#include <stdio.h>

template<class T, int NX, int NY, int NZ>
class Array3d {
protected:
	T* data;
	int index(int i, int j, int k) const {
		assert(data!= NULL);
		return (i + NX * (j + k * NY));
	}
public:
	bool is_allocated() const {
		return data != NULL;
	}
	T interp(int i2, int j2, int k2) const {
		static const double interp[3][3][3] = { { { -8.239746093750E-004, +8.239746093750E-003, +1.373291015625E-003 }, { +8.239746093750E-003,
				-8.239746093750E-002, -1.373291015625E-002 }, { +1.373291015625E-003, -1.373291015625E-002, -2.288818359375E-003 } }, { { +8.239746093750E-003,
				-8.239746093750E-002, -1.373291015625E-002 }, { -8.239746093750E-002, +8.239746093750E-001, +1.373291015625E-001 }, { -1.373291015625E-002,
				+1.373291015625E-001, +2.288818359375E-002 } }, { { +1.373291015625E-003, -1.373291015625E-002, -2.288818359375E-003 }, { -1.373291015625E-002,
				+1.373291015625E-001, +2.288818359375E-002 }, { -2.288818359375E-003, +2.288818359375E-002, +3.814697265625E-003 } } };
		T sum = 0.0;
		for (int l = 0; l < 3; l++) {
			int i = (i2 >> 1) + (l - 1) * (((i2 & 1) << 1) - 1);
			for (int n = 0; n < 3; n++) {
				int j = (j2 >> 1) + (n - 1) * (((j2 & 1) << 1) - 1);
				for (int m = 0; m < 3; m++) {
					int k = (k2 >> 1) + (m - 1) * (((k2 & 1) << 1) - 1);
					sum += (*this)(i, j, k) * interp[l][n][m];

				}
			}
		}
		return sum;
	}
	static void mpi_datatype(MPI_Datatype* newtype, int xlb, int xub, int ylb, int yub, int zlb, int zub, MPI_Datatype oldtype) {
		int sizes[3] = { NX, NY, NZ };
		int starts[3];
		int subsizes[3];
		starts[0] = xlb;
		starts[1] = ylb;
		starts[2] = zlb;
		subsizes[0] = xub - xlb + 1;
		subsizes[1] = yub - ylb + 1;
		subsizes[2] = zub - zlb + 1;
		MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, oldtype, newtype);
		MPI_Type_commit(newtype);
	}
	Array3d(bool alloc_now = false) {
		data = NULL;
		if (alloc_now) {
			allocate();
		}
	}
	void allocate() {
		if (data == NULL) {
			data = new T[NX * NY * NZ];
		}
	}
	void deallocate() {
		if (data != NULL) {
			delete[] data;
		}
		data = NULL;
	}
	Array3d(const Array3d<T, NX, NY, NZ>& a) {
		*this = a;
	}
	Array3d<T, NX, NY, NZ>& operator=(const Array3d<T, NX, NY, NZ>& a) {
		if (&a != this) {
			assert( data != NULL);
			for (int k = 0; k < NZ; k++) {
				for (int j = 0; j < NY; j++) {
					for (int i = 0; i < NX; i++) {
						data[index(i, j, k)] = a.data[index(i, j, k)];
					}
				}
			}
		}
		return *this;
	}
	Array3d<T, NX, NY, NZ>& operator/=(const Array3d<T, NX, NY, NZ>& a) {
		for (int k = 0; k < NZ; k++) {
			for (int j = 0; j < NY; j++) {
				for (int i = 0; i < NX; i++) {
					data[index(i, j, k)] /= a.data[index(i, j, k)];
				}
			}
		}
		return *this;
	}
	Array3d<T, NX, NY, NZ>& operator*=(const Array3d<T, NX, NY, NZ>& a) {
		for (int k = 0; k < NZ; k++) {
			for (int j = 0; j < NY; j++) {
				for (int i = 0; i < NX; i++) {
					data[index(i, j, k)] *= a.data[index(i, j, k)];
				}
			}
		}
		return *this;
	}
	Array3d<T, NX, NY, NZ>& operator+=(const Array3d<T, NX, NY, NZ>& a) {
		for (int k = 0; k < NZ; k++) {
			for (int j = 0; j < NY; j++) {
				for (int i = 0; i < NX; i++) {
					data[index(i, j, k)] += a.data[index(i, j, k)];
				}
			}
		}
		return *this;
	}
	Array3d<T, NX, NY, NZ>& operator-=(const Array3d<T, NX, NY, NZ>& a) {
		for (int k = 0; k < NZ; k++) {
			for (int j = 0; j < NY; j++) {
				for (int i = 0; i < NX; i++) {
					data[index(i, j, k)] -= a.data[index(i, j, k)];
				}
			}
		}
		return *this;
	}
	Array3d<T, NX, NY, NZ>& operator+=(const T a) {
		for (int k = 0; k < NZ; k++) {
			for (int j = 0; j < NY; j++) {
				for (int i = 0; i < NX; i++) {
					data[index(i, j, k)] += a;
				}
			}
		}
		return *this;
	}
	Array3d<T, NX, NY, NZ>& operator-=(const T a) {
		for (int k = 0; k < NZ; k++) {
			for (int j = 0; j < NY; j++) {
				for (int i = 0; i < NX; i++) {
					data[index(i, j, k)] -= a;
				}
			}
		}
		return *this;
	}
	Array3d<T, NX, NY, NZ>& operator=(const T a) {
		for (int k = 0; k < NZ; k++) {
			for (int j = 0; j < NY; j++) {
				for (int i = 0; i < NX; i++) {
					data[index(i, j, k)] = a;
				}
			}
		}
		return *this;
	}
	Array3d<T, NX, NY, NZ>& operator*=(const T a) {
		for (int k = 0; k < NZ; k++) {
			for (int j = 0; j < NY; j++) {
				for (int i = 0; i < NX; i++) {
					data[index(i, j, k)] *= a;
				}
			}
		}
		return *this;
	}
	Array3d<T, NX, NY, NZ>& operator/=(const T a) {
		for (int k = 0; k < NZ; k++) {
			for (int j = 0; j < NY; j++) {
				for (int i = 0; i < NX; i++) {
					data[index(i, j, k)] /= a;
				}
			}
		}
		return *this;
	}
	T operator()(const Vector<int, 3>& v) const {
		return (*this)(v[0], v[1], v[2]);
	}
	T& operator()(const Vector<int, 3>& v) {
		return (*this)(v[0], v[1], v[2]);
	}
	T& operator()(int i, int j, int k) {
		assert( i < NX);assert( j < NY);assert( k < NZ);assert( i >= 0);assert( j >= 0);assert( k >= 0);
		return data[index(i, j, k)];
	}
	const T operator()(int i, int j, int k) const {
		assert( i < NX);assert( j < NY);assert( k < NZ);assert( i >= 0);assert( j >= 0);assert( k >= 0);
		return data[index(i, j, k)];
	}
	T divergence(int i, int j, int k) const {
		assert( i < NX-1);assert( j < NY-1);assert( k < NZ-1);assert( i >= 1);assert( j >= 1);assert( k >= 1);
		return (data[index(i + 1, j, k)] + data[index(i - 1, j, k)] + data[index(i, j + 1, k)] + data[index(i, j - 1, k)] + data[index(i, j, k + 1)]
				+ data[index(i, j, k - 1)]) - (data[index(i, j, k)] * 6.0);
	}
	T oct_avg(int i, int j, int k) const {
		assert( i < NX-1);assert( j < NY-1);assert( k < NZ-1);assert( i >= 1);assert( j >= 1);assert( k >= 1);
		return (data[index(i, j, k)] + data[index(i + 1, j, k)] + data[index(i, j + 1, k)] + data[index(i + 1, j + 1, k)] + data[index(i, j, k + 1)]
				+ data[index(i + 1, j, k + 1)] + data[index(i, j + 1, k + 1)] + data[index(i + 1, j + 1, k + 1)]) * 0.125;
	}
	~Array3d() {
		if (data != NULL) {
			delete[] data;
		}
		data = NULL;
	}
	T* ptr(int i, int j, int k) {
		assert( data != NULL);
		return &(data[index(i, j, k)]);
	}
	T* ptr() {
		assert( data != NULL);
		return data;
	}
	const T* ptr() const {
		assert( data != NULL);
		return data;
	}
	int size() const {
		return NX * NY * NZ;
	}
};

template<int N>
class LinOp {
private:
	Real Jp[3], Jm[3];
	Real K0;
public:
	LinOp() {
		for (int i = 0; i < 3; i++) {
			Jp[i] = Jm[i] = 1.0;
		}
		K0 = 0.0;
	}
	~LinOp() {
	}
	Real& D(int i) {
		assert(i!=0);
		if (i < 0) {
			return Jm[-i - 1];
		} else {
			return Jp[i - 1];
		}
	}
	Real D(int i) const {
		assert(i!=0);
		if (i < 0) {
			return Jm[-i - 1];
		} else {
			return Jp[i - 1];
		}
	}
	Real K() const {
		return K0;
	}
	Real& K() {
		return K0;
	}
	Real J0() const {
		Real sum = K0;
		for (int i = 0; i < 3; i++) {
			sum -= (Jp[i] + Jm[i]);
		}
		return sum;
	}
	Real op(const Array3d<Real, N, N, N>& A, int i, int j, int k) const {
		Real sum = J0() * A(i, j, k);
		sum += A(i - 1, j, k) * Jm[0];
		sum += A(i, j - 1, k) * Jm[1];
		sum += A(i, j, k - 1) * Jm[2];
		sum += A(i + 1, j, k) * Jp[0];
		sum += A(i, j + 1, k) * Jp[1];
		sum += A(i, j, k + 1) * Jp[2];
		return sum;
	}
};

#endif /* ARRAY3D_H_ */
