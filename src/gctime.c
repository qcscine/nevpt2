#include <string.h>
#include <time.h>

#define CH_F2C(X) ((char *) (X))

int gctime(char *fstr, int lstr) {
    time_t t;
    t = time(NULL);
    strcpy(CH_F2C(fstr), ctime(&t));
    return 0;
}

int gctime_(char *fstr, int lstr) {
    time_t t;
    t = time(NULL);
    strcpy(CH_F2C(fstr), ctime(&t));
    return 0;
}
