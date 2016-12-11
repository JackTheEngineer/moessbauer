/* ----------------------------------------------------------------------- */
/* file            NDOUB.C                                                 */
/* function(s)     Theory for N quadrupole doublets                        */
/* author          Gerhard Grosse                                          */
/* last change     07.10.91                                                */
/*             (C) Copyright 1991 by Gerhard Grosse                        */
/* ----------------------------------------------------------------------- */


/* ----------------------------------------------------------------------- */
/* Vital: Include overhead for Lorentzian theories:                        */
/* ------------------------------------------------                        */
 # include  "lorentz.c"
/* ----------------------------------------------------------------------- */


/* ----------------------------------------------------------------------- */
   static int  nlrtz, npara, ndoub;
/* Number of Lorentzians, parameters and doublets.                         */
/* ----------------------------------------------------------------------- */


/* ----------------------------------------------------------------------- */
 # pragma argsused                      /* supress an annoying warning ... */
   static int  no_lrtz  (int     np,
                         double  p[])
/*                                                                         */
/* Determine number of Lorentzians and doublets.                           */
/* ----------------------------------------------------------------------- */

{
   npara = np;
   ndoub = (np-1) / 4;
   return nlrtz = ndoub * 2;
}


/* ----------------------------------------------------------------------- */
   static void  mk_lrtz  (double  p[],
                          lrtz_t  l[])
/*                                                                         */
/* Calculate the Lorentzians given the parameter array <p>.                */
/*                                                                         */
/* Theory for "ndoub" quadrupole doublets.                                 */
/*                                                                         */
/* General Parameters:                                                     */
/*                                                                         */
/* p[0]       = Baseline                                                   */
/* p[1]       = Total area                                                 */
/*                                                                         */
/* Parameters for the 1st doublet:                                         */
/*                                                                         */
/* p[2+0]     = Quadrupole splitting (eQVzz/2 in mm/s)                     */
/* p[2+1]     = Isomer shift (mm/s)                                        */
/* p[2+2]     = Linewidth (mm/s)                                           */
/*                                                                         */
/* Parameters for (k+1)'st doublet (1 <= k < ndoub):                       */
/*                                                                         */
/* p[1+4*k+0] = Quadrupole splitting (eQVzz/2 in mm/s)                     */
/* p[1+4*k+1] = Isomer shift (mm/s)                                        */
/* p[1+4*k+2] = Linewidth relative to first doublet                        */
/* p[1+4*k+3] = Partial intensity                                          */
/* ----------------------------------------------------------------------- */

{
   int      k;
   lrtz_t  *_l;
   double  *_p;

   _l = l;
   _p = p + 2;                              /* parameters of first doublet */
   _l[0].lpos = _p[1] + 0.5 * _p[0];
   _l[1].lpos = _p[1] - 0.5 * _p[0];
   _l[0].wdth = _l[1].wdth = _p[2];
   _l[0].area = 0.5;
   for (k = 1; k < ndoub; k++) {
     _l = l + 2*k;                      /* Lorentzians of (k+1)-st doublet */
     _p = p + 1 + 4*k;
     _l[0].lpos  = _p[1] + 0.5 * _p[0];
     _l[1].lpos  = _p[1] - 0.5 * _p[0];
     _l[0].wdth  = _l[1].wdth =  _p[2] * l[0].wdth;
     _l[0].area  = _l[1].area =   0.5  * _p[3];
     _l[0].dpth  = _l[1].dpth = _l[0].area / _l[0].wdth;
     l[0].area  -= _l[0].area;
   }
   l[1].area = l[0].area;
   l[0].dpth = l[1].dpth = l[0].area / l[0].wdth;
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
  if (n > ndoub)
    return 0;
  else {
    subl[0] = 2*(n-1);
    subl[1] = 2*(n-1) + 1;
    return 2;
  }
}


/* ----------------------------------------------------------------------- */
/* End of file NDOUB.C                                                     */
/* ----------------------------------------------------------------------- */