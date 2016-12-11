// -----------------------------------------------------------------------
// file            QUADRUPO.C
// function(s)     Beispiel-Theorieprogramm fuer Moessbauer-Praktikum
// author          Gerhard Grosse
// version         1.0.00
// created         05.11.95
// revised
// -----------------------------------------------------------------------


// -----------------------------------------------------------------------
// Parameter Definitionen:
// -----------------------

 # define NLINES       2         // Anzahl der Lorentz-Linien
 # define NPARAMETERS  5         // Anzahl der Fit-Parameter

 # define BASELINE   par[0]      // Darf nicht geaendert werden!
 # define QUADSPLIT  par[1]      // eQVzz/2: Quadrupolaufspaltung
 # define ISOSHIFT   par[2]      // Isomerieverschiebung
 # define LINEWIDTH  par[3]      // Linienbreite
 # define LINEDEPTH  par[4]      // Linientiefe (max. Absorption)

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
   // Linienposition:

   lor[0].lpos = ISOSHIFT - 0.5 * QUADSPLIT;
   lor[1].lpos = ISOSHIFT + 0.5 * QUADSPLIT;

   // Linienbreiten (beide Linien sind gleich breit):

   lor[0].wdth = LINEWIDTH;
   lor[1].wdth = LINEWIDTH;

   // Linientiefen (beide Linien sind gleich tief):

   lor[0].dpth = LINEDEPTH;
   lor[1].dpth = LINEDEPTH;
}


// -----------------------------------------------------------------------
// End of file QUADRUPO.C
// -----------------------------------------------------------------------
