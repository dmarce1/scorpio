#ifndef VIRTUAL_PROCESS_H_
#define VIRTUAL_PROCESS_H_

template<class T>
class VirtualProcess {
protected:
    int ip;
    typedef void (T::*ifunc_t)(int);
    int get_instruction_pointer() const {
        return ip;
    }
    void set_instruction_pointer(int i) {
        ip = i;
    }
    static void run_program(T** procs, int nprocs, ifunc_t* cs, int i_cnt) {
        int last_ip;
        bool done;
        for (int i = 0; i < nprocs; i++) {
            procs[i]->ip = 0;
        }
        do {
            done = true;
            for (int i = 0; i < nprocs; i++) {
                if (procs[i]->ip < i_cnt) {
                    do {
                        last_ip = procs[i]->ip;
                        (procs[i]->*cs[procs[i]->ip])(0);
                    } while (last_ip != procs[i]->ip && procs[i]->ip < i_cnt);
                    if (procs[i]->ip < i_cnt) {
                        done = false;
                    }
                }
            }
        } while (!done);
    }
    void inc_instruction_pointer() {
        ip++;
    }
    VirtualProcess() {
        return;
    }
    virtual ~VirtualProcess() {
        return;
    }
};

#endif /* VIRTUAL_PROCESS_H_ */
