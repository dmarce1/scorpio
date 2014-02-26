/*
 * parameter_reader.cpp
 *
 *  Created on: Feb 26, 2014
 *      Author: dmarce1
 */

#include <string.h>
#include <stdlib.h>

bool read_parameter(int argc, char* argv[], const char* name, double* param) {
    bool rc;
    int namelen;
    namelen = strlen(name);
    rc = false;
    for (int i = 1; i < argc; i++) {
        if (argv[i][0] == '-') {
            if (strncmp(name, argv[i] + 1, namelen) == 0) {
                rc = true;
                if (param) {
                    *param = atof(argv[i] + namelen + 2);
                }
            }
        }
    }
    return rc;
}

bool read_parameter(int argc, char* argv[], const char* name, int* param) {
    bool rc;
    int namelen;
    namelen = strlen(name);
    rc = false;
    for (int i = 1; i < argc; i++) {
        if (argv[i][0] == '-') {
            if (strncmp(name, argv[i] + 1, namelen) == 0) {
                rc = true;
                if (param) {
                    *param = atoi(argv[i] + namelen + 2);
                }
            }
        }
    }
    return rc;
}

bool read_parameter(int argc, char* argv[], const char* name, char* param, int strl) {
    bool rc;
    int namelen;
    namelen = strlen(name);
    rc = false;
    for (int i = 1; i < argc; i++) {
        if (argv[i][0] == '-') {
            if (strncmp(name, argv[i] + 1, namelen) == 0) {
                rc = true;
                strncpy(param, argv[i] + namelen + 2, strl);
            }
        }
    }
    return rc;
}
