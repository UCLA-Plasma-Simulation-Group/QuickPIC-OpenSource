#include <stdio.h>
#include <sys/time.h>

void dtimer(double *time, struct timeval *itime, int icntrl) {
/* this subroutine performs timing
   input: icntrl, itime
   icntrl = (-1,0,1) = (initialize,ignore,read) clock
   clock should be initialized before it is read!
   time = elapsed time in seconds                   */
   long oss, usec;
   const double tick = 1.0/1000000.0;
   struct timeval jclock, nclock;
   if (!(icntrl)) 
      return;
   if (icntrl==1)
      goto L10;
/* initialize clock              
   calculate time elapsed in microseconds */
   oss = gettimeofday(&jclock,NULL);
   *itime = jclock;
   return;
/* read clock and write time difference from last clock initialization
   calculate time elapsed in microseconds                              */
L10: oss = gettimeofday(&nclock,NULL);
   jclock = *itime;
   usec = nclock.tv_usec - jclock.tv_usec;
   *time = (nclock.tv_sec - jclock.tv_sec) + usec*tick;
   return;
}

void dtimer_(double *time, unsigned long *itime, int *icntrl) {
/* in Fortran, itime is an array of two or four integers */
   dtimer(time,(struct timeval *)itime,*icntrl);
   return;
}
