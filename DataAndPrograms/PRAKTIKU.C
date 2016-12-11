// -----------------------------------------------------------------------
// file            PRAKTIKU.C
// function(s)     Include file fuer Moessbauer lab course
// author          Gerhard Grosse
// version         1.0.00
// created         06.11.95
// revised
// -----------------------------------------------------------------------


// -----------------------------------------------------------------------
// Vital: Include overhead for Lorentzian theories:
// ------------------------------------------------
 # include  "lorentz.c"
// -----------------------------------------------------------------------


// -----------------------------------------------------------------------
   static int  nLrtz, nPara;
// Number of Lorentzians and parameters.
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
   static void Theory (double par[], lrtz_t lor[]);
// Prototype of user routine
// -----------------------------------------------------------------------


// -----------------------------------------------------------------------
 # pragma argsused    // supress an annoying warning ...
   static int  no_lrtz  (int     np,
			 double  p[])
//
// Determine number of Lorentzians.
// -----------------------------------------------------------------------

{
   nPara = np;
   nLrtz = NLINES;
   if (nPara != NPARAMETERS)
      error(0,"Number of parameters does not match theory program");
   return nLrtz;
}


// -----------------------------------------------------------------------
   static void  mk_lrtz  (double  p[],
			  lrtz_t  l[])
//
// Calculate the Lorentzians given the parameter array <p>.
// -----------------------------------------------------------------------

{
   Theory(p,l);
}


// -----------------------------------------------------------------------
   static int  subspec  (int  n,
			 int  subl[])
//
// Return the no. of Lorentzians belonging to the n'th subspectrum and
// write the indices of the Lorentzians to <subl>.
// -----------------------------------------------------------------------

{
  if (n > nLrtz)
     return 0;
  else {
     subl[0] = n-1;
     return 1;
  }
}


// -----------------------------------------------------------------------
// End of file PRAKTIKU.C
// -----------------------------------------------------------------------
