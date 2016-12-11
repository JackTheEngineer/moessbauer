// -----------------------------------------------------------------------
// file            TEMPLATE.C
// function(s)     Theorieprogramm-Schablone fuer Moessbauer-Praktikum
// author          Gerhard Grosse
// version         1.0.00
// created         06.11.95
// revised
// -----------------------------------------------------------------------


// -----------------------------------------------------------------------
// Parameter Definitionen:
// -----------------------

 # define NLINES                 // Anzahl der Lorentz-Linien
 # define NPARAMETERS            // Anzahl der Fit-Parameter

 # define BASELINE   par[0]      // Darf nicht geaendert werden!

// Weitere Parameter ...

// -----------------------------------------------------------------------


// -----------------------------------------------------------------------
// Einbinden der Schnittstelle zu Mos-90:
// --------------------------------------

 # include  "praktiku.c"

// -----------------------------------------------------------------------


// -----------------------------------------------------------------------
// Die eigentliche Theorie-Funktion:
// ---------------------------------

   static void  Theory  (double par[], lrtz_t lor[])

{
   // Linienpositionen ...


   // Linienbreiten ...


   // Linientiefen ...

}


// -----------------------------------------------------------------------
// End of file TEMPLATE.C
// -----------------------------------------------------------------------
