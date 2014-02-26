/*
 * parameter_reader.h
 *
 *  Created on: Feb 26, 2014
 *      Author: dmarce1
 */

#ifndef PARAMETER_READER_H_
#define PARAMETER_READER_H_

bool read_parameter(int argc, char* argv[], const char* parameter_name, double* param = 0);
bool read_parameter(int argc, char* argv[], const char* parameter_name, int* param);
bool read_parameter(int argc, char* argv[], const char* parameter_name, char* param, int);

#endif /* PARAMETER_READER_H_ */
