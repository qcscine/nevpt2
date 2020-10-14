typedef long f77_int;     /* Fortran integer type */
typedef char * f77_char;    /* Fortran character argument */
#define CH_F2C(X) ((char *) (X))  /* How to get char ptr from F77 argument */
gctime (fstr, lstr) f77_char *fstr; int lstr; {
   long time(), t;
   char *ctime();
   t = time ( (long *) 0);
   strcpy(CH_F2C(fstr),ctime(&t));
   return (0);
   }
gctime_(fstr, lstr) f77_char *fstr; int lstr; {
   long time(), t;
   char *ctime();
   t = time ( (long *) 0);
   strcpy(CH_F2C(fstr),ctime(&t));
   return (0);
   }
