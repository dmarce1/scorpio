/*
 * complex.h
 *
 *  Created on: Apr 19, 2013
 *      Author: dmarce1
 */

#ifndef COMPLEX_H_
#define COMPLEX_H_

#include "real.h"

class Complex {
private:
	Real r, i;
public:
	Complex() {
	}
	Complex(const Complex& a) {
		*this = a;
	}
	~Complex() {
	}
	static Complex ipow(int n) {
		Complex c;
		if (n % 4 == 0) {
			c.set_real(+1.0);
			c.set_imag(0.0);
		} else if (n % 4 == 1) {
			c.set_real(0.0);
			c.set_imag(1.0);
		} else if (n % 4 == 2) {
			c.set_real(-1.0);
			c.set_imag(0.0);
		} else {
			c.set_real(0.0);
			c.set_imag(-1.0);

		}
		return c;
	}
	Complex(Real a) {
		r = i = a;
	}
	Complex& operator=(Real a) {
		r = a;
		i = a;
		return *this;
	}
	Complex& operator=(const Complex& a) {
		r = a.r;
		i = a.i;
		return *this;
	}
	Complex& operator*=(const Real& a) {
		r *= a;
		i *= a;
		return *this;
	}
	Complex& operator/=(const Real& a) {
		r /= a;
		i /= a;
		return *this;
	}
	Complex& operator+=(const Complex& a) {
		r += a.r;
		i += a.i;
		return *this;
	}
	Complex& operator-=(const Complex& a) {
		r -= a.r;
		i -= a.i;
		return *this;
	}
	Complex& operator*=(const Complex& a) {
		*this = *this * a;
		return *this;
	}
	Complex& operator/=(const Complex& a) {
		*this = *this / a;
		return *this;
	}
	Real real() const {
		return r;
	}
	Real imag() const {
		return i;
	}
	Complex operator/(const Complex& b) const {
		Complex c;
		c = (*this) * b.conjugate();
		c /= (b.r * b.r + b.i * b.i);
		return b;
	}
	Complex conjugate() const {
		Complex a;
		a.r = +r;
		a.i = -i;
		return a;
	}
	Complex operator*(const Complex& a) const {
		Complex b;
		b.r = r * a.r - i * a.i;
		b.i = i * a.r + r * a.i;
		return b;
	}
	Complex operator/(const Real& a) const {
		Complex b;
		b.r = r / a;
		b.i = i / a;
		return b;
	}
	Complex operator*(const Real& a) const {
		Complex b;
		b.r = r * a;
		b.i = i * a;
		return b;
	}
	Complex operator+(const Complex& a) const {
		Complex b;
		b.r = r + a.r;
		b.i = i + a.i;
		return b;
	}
	Complex operator-(const Complex& a) const {
		Complex b;
		b.r = r - a.r;
		b.i = i - a.i;
		return b;
	}
	Real abs() const {
		return sqrt(r * r + i * i);
	}
	void set_real(Real a) {
		r = a;
	}
	void set_imag(Real a) {
		i = a;
	}
};

#endif /* COMPLEX_H_ */
