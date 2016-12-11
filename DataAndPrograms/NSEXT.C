/* ----------------------------------------------------------------------- */
/* file            NSEXT.C                                                 */
/* function        Theory for N magnetic sextets and M doublets            */
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
   static int  nlrtz, npara, ndoub, nsex;
/* Number of Lorentzians, parameters, doublets and sextets.                */
/* ----------------------------------------------------------------------- */


/* ----------------------------------------------------------------------- */
   static int  no_lrtz  (int     np,
                         double  p[])
/*                                                                         */
/* Determine number of Lorentzians, doublets and sextets.                  */
/* ----------------------------------------------------------------------- */

{
   npara = np;
   nsex  = (int) p[2];
   ndoub = (int) p[3];
   if (npara != 3 + 7*nsex + 4*ndoub)
     return 0;
   else
     return nlrtz = ndoub * 2 + nsex * 6;
}


/* ----------------------------------------------------------------------- */
   static void  mk_lrtz  (double  p[],
                          lrtz_t  l[])
/*                                                                         */
/* Calculate the Lorentzians <l> given the parameter array <p>.            */
/* Theory for <nsex> magnetic sextets of Fe57 and <ndoub> quadrupole       */
/* doublets. Electrostatic quadrupole splitting of the sextets is treated  */
/* in 1st order pertubation theory. A small (assumed Lorentzian shaped)    */
/* distribution of magnetic hyperfine fields within one sextet is allowed  */
/* for by a correlated widening of the sextet lines.                       */
/*                                                                         */
/* General parameters:                                                     */
/*                                                                         */
/* p[0]         = Baseline                                                 */
/* p[1]         = Total area                                               */
/* p[2]         = Number of sextets  ( >= 1)                               */
/* p[3]         = Number of doublets ( >= 0)                               */
/*                                                                         */
/* Parameters for 1st sextet:                                              */
/*                                                                         */
/* p[4+0]       = Magnetic hyperfine field (Tesla)                         */
/* p[4+1]       = Quadrupole splitting (eQVzz/2 in mm/s)                   */
/* p[4+2]       = Isomer shift (mm/s)                                      */
/* p[4+3]       = Width of outer lines (mm/s)                              */
/* p[4+4]       = Sets intensity ratio to 3 : 2 * p[8] : 1 * (p[8])**2     */
/* p[4+5]       = Natural line width (mm/s)                                */
/*                                                                         */
/* Parameters for the (k+1)st sextet (k >= 1):                             */
/*                                                                         */
/* p[3+7k+0]    = Magnetic hyperfine field (Tesla)                         */
/* p[3+7k+1]    = Quadrupole splitting (eQVzz/2 in mm/s)                   */
/* p[3+7k+2]    = Isomer shift (mm/s)                                      */
/* p[3+7k+3]    = Width of outer lines (mm/s)                              */
/* p[3+7k+4]    = Sets intensity ratio to 3 : 2 * p : 1 * p**2             */
/* p[3+7k+5]    = Partial intensity of the sextet                          */
/* p[3+7k+6]    = Natural line width (mm/s)                                */
/*                                                                         */
/* Parameters for the 1st doublet (if one) (assume n sextets):             */
/*                                                                         */
/* p[3+7n+0]    = Quadrupole splitting (eQVzz/2 in mm/s)                   */
/* p[3+7n+1]    = Isomer shift (mm/s)                                      */
/* p[3+7n+2]    = Linewidth (mm/s)                                         */
/* p[3+7n+3]    = Partial intensity                                        */
/*                                                                         */
/* Parameters for the (i+1)st doublet (i >= i) (assume n sextets):         */
/*                                                                         */
/* p[3+7n+4i+0] = Quadrupole splitting (eQVzz/2 in mm/s)                   */
/* p[3+7n+4i+1] = Isomer shift (mm/s)                                      */
/* p[3+7n+4i+2] = Linewidth relativ to 1st doublet                         */
/* p[3+7n+4i+3] = Partial intensity                                        */
/* ----------------------------------------------------------------------- */

{
   lrtz_t  *_l;
   double  *_p;
   int      i;

   double  hf, qm, sm, amag, wnat;

   double  dg     = 0.059424;                 /* theoretical constants ... */
   double  de     = -0.57113 * dg;
   double  magou  = dg - 3 * de;
   double  magmi  = dg - de;
   double  magin  = dg + de;

   double  area   = 0.0;                 /* adds up fractional intensities */

   for (i = 1; i < nsex; i++) {              /* treat sextets (except 1st) */

     _l = l + 6*i;                   /* first Lorentzian of (i+1)st sextet */
     _p = p + 3 + 7*i;                     /* parameters of (i+1)st sextet */

     hf   = _p[0];
     qm   = _p[1] * 0.5;
     sm   = _p[2];
     amag = _p[5] / (1.0 + (0.66667 + 0.33333 * _p[4]) * _p[4]);
     wnat = _p[6];

     _l[0].lpos = sm + qm + hf * magou;
     _l[1].lpos = sm + qm - hf * magou;
     _l[2].lpos = sm - qm + hf * magmi;
     _l[3].lpos = sm - qm - hf * magmi;
     _l[4].lpos = sm - qm + hf * magin;
     _l[5].lpos = sm - qm - hf * magin;

     _l[0].wdth = _l[1].wdth = _p[3];
     _l[2].wdth = _l[3].wdth = wnat + (_p[3] - wnat) * magmi / magou;
     _l[4].wdth = _l[5].wdth = wnat + (_p[3] - wnat) * magin / magou;

     _l[0].area = _l[1].area = 0.5 * amag;
     _l[2].area = _l[3].area = 0.33333 * _p[4] * amag;
     _l[4].area = _l[5].area = 0.16667 * _p[4] * _p[4] * amag;

     area += _p[5];
   }

   for (i = 0; i < ndoub; i++) {                         /* treat doublets */

     _l = l + 6*nsex + 2*i;         /* first Lorentzian of (i+1)st doublet */
     _p = p + 3 + 7*nsex + 4*i;          /* parameters for (i+1)st doublet */

     _l[0].lpos = _p[1] + 0.5 * _p[0];
     _l[1].lpos = _p[1] - 0.5 * _p[0];

     _l[0].wdth = _l[1].wdth = (i == 0) ? _p[2] : p[3+7*nsex+2] * _p[2];
     _l[0].area = _l[1].area = 0.5 * _p[3];

     area += _p[3];
   }

   /* treat 1st sextet: */

   _p = p + 4;                          /* first parameter of first sextet */

   hf   = _p[0];
   qm   = _p[1] * 0.5;
   sm   = _p[2];
   amag = (1.0 - area) / (1.0 + (0.66667 + 0.33333 * _p[4]) * _p[4]);
   wnat = _p[5];

   l[0].lpos = sm + qm + hf * magou;
   l[1].lpos = sm + qm - hf * magou;
   l[2].lpos = sm - qm + hf * magmi;
   l[3].lpos = sm - qm - hf * magmi;
   l[4].lpos = sm - qm + hf * magin;
   l[5].lpos = sm - qm - hf * magin;

   l[0].wdth = l[1].wdth = _p[3];
   l[2].wdth = l[3].wdth = wnat + (_p[3] - wnat) * magmi / magou;
   l[4].wdth = l[5].wdth = wnat + (_p[3] - wnat) * magin / magou;

   l[0].area = l[1].area = 0.5 * amag;
   l[2].area = l[3].area = 0.33333 * _p[4] * amag;
   l[4].area = l[5].area = 0.16667 * _p[4] * _p[4] * amag;

   for (i = 0; i < nlrtz; i++)                 /* depth of all Lorentzians */
     l[i].dpth = p[1] * l[i].area / l[i].wdth;

}


/* ----------------------------------------------------------------------- */
   static int  subspec  (int   n,
                         int   subl[])
/*                                                                         */
/* Return the no. of Lorentzians belonging to the <n>'th subspectrum and   */
/* write the indices of the Lorentzians to <subl>.                         */
/* ----------------------------------------------------------------------- */

{
   int   i;

   if (n > nsex + ndoub)
     return 0;
   else if (n > nsex) {                          /* subspectrum is doublet */
     subl[0] = 6 * nsex + 2 * (n - nsex - 1);
     subl[1] = subl[0] + 1;
     return 2;
   }
   else {                                         /* subspectrum is sextet */
     subl[0] = 6 * (n - 1);
     for (i = 1; i < 6; i++) subl[i] = subl[0] + i;
     return 6;
   }
}


/* ----------------------------------------------------------------------- */
/* End of file NSEXT.C                                                     */
/* ----------------------------------------------------------------------- */