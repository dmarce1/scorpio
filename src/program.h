#ifndef PROGRAM_H_
#define PROGRAM_H_

class Program {
public:
	Program();
	static void print_help();
	int run(int, char* a[]);
	virtual ~Program();
};

#endif /* PROGRAM_H_ */
