#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <mkl.h>

/* --------------------- INCLUDE PARTS --------------------- */
#include"headjmart.c"
#include"loadjmart.c"
#include"movemart.c"
#include"markmart.c"
#include"gausmart.c"
#include"heatmart.c"
#include"impact.c"    /* Gregor add */
/* --------------------- INCLUDE PARTS --------------------- */




/* Solve differencial equations by MARKER-IN-CELL + FINITE-DIFFERENCES method */
int main()
{
/* Counters */
int n0,n1,f0;
long int pos0cur0,m1,m2,m3;
double dx,dy,dx1,dx2,dx3;
/**/
/**/
/**/
/* Load configuration from mode.t3c */
loadconf();
/**/
/**/
/**/
/* Output File Cycle */
for (f0=0;f0<fl0num;f0++)
	{
	/* Reload Cur Input File Name, Type */
	for (n1=0;n1<50;n1++) fl1in[n1]=fl0in[f0][n1];
	fl1itp=fl0itp[f0];
	/**/
	/* Load data from input file */
	loader();
	/**/
	/* Postprocessing */
	/* Reload Cur Output File Name, Type */
	for (n1=0;n1<50;n1++) fl1out[n1]=fl0out[f0][n1];
	fl1otp=fl0otp[f0];
	/**/
	/* Save data to output file */
	saver(f0);
	}
/* End Program */
return 0;
}
/* Solve differencial equations by MARKER-IN-CELL + FINITE-DIFFERENCES method */
