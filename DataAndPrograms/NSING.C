/* ----------------------------------------------------------------------- */
/* file            NSING.C                                                 */
/* function(s)     Theory for N single Lorentzians                         */
/* author          Gerhard Grosse                                          */
/* last change     07.09.91                                                */
/*             (C) Copyright 1991 by Gerhard Grosse                        */
/* ----------------------------------------------------------------------- */


/* ----------------------------------------------------------------------- */
/* Vital: Include overhead for Lorentzian theories:                        */
/* ------------------------------------------------                        */
 # include  "lorentz.c"
/* ----------------------------------------------------------------------- */


/* ----------------------------------------------------------------------- */
   static int  nlrtz, npara;
/* Number of Lorentzians and parameters.                                   */
/* ----------------------------------------------------------------------- */


/* ----------------------------------------------------------------------- */
 # pragma argsused                      /* supress an annoying warning ... */
   static int  no_lrtz  (int     np,
                         double  p[])
/*                                                                         */
/* Determine number of Lorentzians.                                        */
/* ----------------------------------------------------------------------- */

{
   npara = np;
   nlrtz = (np-1) / 3;
   return nlrtz;
}


/* ----------------------------------------------------------------------- */
   static void  mk_lrtz  (double  p[],
                          lrtz_t  l[])
/*                                                                         */
/* Calculate the Lorentzians given the parameter array <p>.                */
/*                                                                         */
/* Theory for <nlrtz> single Lorentzians.                                  */
/*                                                                         */
/* General Parameters:                                                     */
/*                                                                         */
/* p[0]       = Baseline                                                   */
/* p[1]       = Total area                                                 */
/*                                                                         */
/* Parameters for the 1st Lorentzian:                                      */
/*                                                                         */
/* p[2+0]     = Line position (mm/s)                                       */
/* p[2+1]     = Linewidth (mm/s)                                           */
/*                                                                         */
/* Parameters for (k+1)'st Lorentzian (1 <= k < nlrtz):                    */
/*                                                                         */
/* p[1+3*k+0] = Line position (mm/s)                                       */
/* p[1+3*k+1] = Linewidth relative to first Lorentzian                     */
/* p[1+3*k+2] = Partial intensity                                          */
/* ----------------------------------------------------------------------- */

{
   int      k;
   lrtz_t  *_l;
   double  *_p;

   _l = l;
   _p = p + 2;                           /* parameters of first Lorentzian */
   _l[0].lpos = _p[0];
   _l[0].wdth = _p[1];
   _l[0].area = 1.0;
   for (k = 1; k < nlrtz; k++) {
     _l = l + k;                                    /* (k+1)'st Lorentzian */
     _p = p + 1 + 3*k;
     _l[0].lpos  = _p[0];
     _l[0].wdth  = _p[1] * l[0].wdth;
     _l[0].area  = _p[2];
     _l[0].dpth  = _l[0].area / _l[0].wdth;
     l[0].area  -= _l[0].area;
   }
   l[0].dpth = l[0].area / l[0].wdth;
   for (k = 0; k < nlrtz; k++)
     l[k].dpth *= p[1];
}


/* ----------------------------------------------------------------------- */
   static int  subspec  (int  n,
                         int  subl[])
/*                                                                         */
/* Return the no. of Lorentzians belonging to the n'th subspectrum and     */
/* write the indices of the Lorentzians to <subl>.                         */
/* ----------------------------------------------------------------------- */

{
  if (n > nlrtz)
    return 0;
  else {
    subl[0] = n-1;
    return 1;
  }
}


/* ----------------------------------------------------------------------- */
/* End of file NSING.C                                                     */
/* ----------------------------------------------------------------------- */