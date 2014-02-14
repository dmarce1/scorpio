#ifndef __MYBITS_H
#define __MYBITS_H

#include <stdio.h>

typedef unsigned char byte_t;

class Bits {
private:
	byte_t* a;
	int cnt;
public:
	Bits(int size) {
		cnt = (size - 1) / 8 + 1;
		a = new byte_t[cnt];
		for (int i = 0; i < cnt; i++) {
			a[i] = 0;
		}
	}
	void print() const {
		for (int i = 0; i < 8 * cnt; i++) {
			if (get(i)) {
				printf("1");
			} else {
				printf("0");
			}
			if ((i+1) % 8 == 0) {
				printf("|");
			}
		}
		printf("\n");
	}
	int size() const {
		return cnt;
	}
	byte_t* ptr() {
		return a;
	}
	bool get(int i) const {
		if (i >= 8 * cnt) {
			printf("%i %i\n", i, cnt);
			assert(i < 8 * cnt);
		}
		return (((a[i / 8] >> (i % 8)) & 0x1) != 0x0);
	}
	void set(int i, bool b) const {
		assert(i < 8 * cnt);
		if (b) {
			a[i / 8] |= (0x1 << (i % 8));
		} else {
			a[i / 8] &= ~(0x1 << (i % 8));
		}
	}
	~Bits() {
		delete[] a;
	}
};

#endif
