#include "reconstruct.h"
#include <stdlib.h>


Reconstruct::Reconstruct() {
	return;
}

Reconstruct::~Reconstruct() {
	return;
}


void Reconstruct::operator()(State q0[GNX], State ql[GNX], State qr[GNX]) const {
    State s[GNX];
    int i, l;
    Real tmp1, tmp2;
    for (i = 1; i < GNX - 1; i++) {
        State up = q0[i + 1] - q0[i];
        State um = q0[i] - q0[i - 1];
        s[i] = minmod_theta(up, um, 2.0);
    }
    for (i = 1; i < GNX - 2; i++) {
        qr[i] = (q0[i] + q0[i + 1]) * 0.5;
        for (l = 0; l < STATE_NF; l++) {
            if (!State::low_order_variable(l)) {
                qr[i][l] += (s[i][l] - s[i + 1][l]) * (1.0 / 6.0);
            }
        }
        ql[i] = qr[i - 1];
    }
    for (i = 2; i < GNX - 2; i++) {
        for (l = 0; l < STATE_NF; l++) {
            if (!State::smooth_variable(l)) {
                tmp1 = qr[i][l] - ql[i][l];
                tmp2 = qr[i][l] + ql[i][l];
                if ((qr[i][l] - q0[i][l]) * (q0[i][l] - ql[i][l]) <= 0.0) {
                    qr[i][l] = ql[i][l] = q0[i][l];
                } else if (tmp1 * (q0[i][l] - 0.5 * tmp2) > (1.0 / 6.0) * tmp1 * tmp1) {
                    ql[i][l] = 3.0 * q0[i][l] - 2.0 * qr[i][l];
                } else if (-(1.0 / 6.0) * tmp1 * tmp1 > tmp1 * (q0[i][l] - 0.5 * tmp2)) {
                    qr[i][l] = 3.0 * q0[i][l] - 2.0 * ql[i][l];
                }
            }
        }
    }
    for (i = GNX - 3; i > 2; i--) {
        qr[i] = qr[i - 1];
    }
}

