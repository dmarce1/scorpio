#ifndef VIRTUAL_PROCESS_H_
#define VIRTUAL_PROCESS_H_



template<class T>
class VirtualProcess {
protected:
	int* ip;
	typedef void (T::*ifunc_t)(int);
	int get_instruction_pointer(int dir = 0) const {
		return ip[dir];
	}
	void set_instruction_pointer(int i, int dir = 0) {
		ip[dir] = i;
	}
	static void run_program(T** procs, int nprocs, ifunc_t* cs, int i_cnt, int nthreads = 1) {
		int last_ip;
		bool done;
		for (int i = 0; i < nprocs; i++) {
			procs[i]->ip = new int[nthreads];
			for (int j = 0; j < nthreads; j++) {
				procs[i]->ip[j] = 0;
			}
		}
		do {
			done = true;
			for (int i = 0; i < nprocs; i++) {
				for (int j = 0; j < nthreads; j++) {
					if (procs[i]->ip[j] < i_cnt) {
						do {
							last_ip = procs[i]->ip[j];
							(procs[i]->*cs[procs[i]->ip[j]])(j);
						} while (last_ip != procs[i]->ip[j] && procs[i]->ip[j] < i_cnt);
						if (procs[i]->ip[j] < i_cnt) {
							done = false;
						}
					}
				}
			}
		} while (!done);
		for (int i = 0; i < nprocs; i++) {
			delete[] procs[i]->ip;
		}
	}
	bool threads_are_synced(int nthreads) {
		if (nthreads == 1) {
			return true;
		} else {
			bool rc = true;
			for (int i = 1; i < nthreads; i++) {
				if (ip[i] != ip[0]) {
					rc = false;
					break;
				}
			}
			return rc;
		}
	}
	void inc_instruction_pointer(int dir = 0) {
		ip[dir]++;
	}
	VirtualProcess() {
		return;
	}
	virtual ~VirtualProcess() {
		return;
	}
};

#endif /* VIRTUAL_PROCESS_H_ */
