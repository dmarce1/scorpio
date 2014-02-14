#ifndef ASSERT_H_
#define ASSERT_H_


#ifndef NDEBUG
#define assert(a) __assert(a,__FILE__,__LINE__)
void __assert(bool, const char*, int);
#else
#include <stdlib.h>
#undef assert
#define assert(a)
#endif

#endif /* ASSERT_H_ */
