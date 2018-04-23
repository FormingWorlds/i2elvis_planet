#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <mkl.h>
/*
#include <omp.h>
*/

/* --------------------- INCLUDE PARTS --------------------- */
#include"headmart.c"
#include"loadmart.c"
#include"movemart.c"
#include"markmart.c"
#include"gausmart.c"
#include"heatmart.c"
#include"impact.c"
#include"core.c"
#include"grain.c"
/* --------------------- INCLUDE PARTS --------------------- */




/* Solve differencial equations by MARKER-IN-CELL + FINITE-DIFFERENCES method */
int main()
{
/* Counters */
int n0,n1,f0,fln3;
long int pos0cur0,m1,m2,m3;
/**/
/**/
/**/
/* Load data from input file */
fln3=loadconf()+1;
/**/
/**/
/**/
/* Load data from input file */
loader();
/**/
/* Read impact history */
impactread();
/**/
/**/
/* Output File Cycle */
for (f0=fln3;f0<fl0num;f0++)
{
/* Reload Cur Output File Name, Type */
for (n1=0;n1<50;n1++) fl1out[n1]=fl0out[f0][n1];
fl1otp=fl0otp[f0];
/**/
/* Reload cyc0max_maxxystep_maxtkstep_maxtmstep */
cyc0max=fl0cyc[f0];
maxxystep=fl0stp[f0][0];
maxtkstep=fl0stp[f0][1];
maxtmstep=fl0stp[f0][2];
/**/
/* General Cycle */
for (n0=0;n0<cyc0max;n0++)
	{
	if (printmod) printf("\n! FILE %s  KRUG ! %d\n",fl1out,n0+1);
	/**/
	/**/
	/**/
	/* Set initial time step */
	timestep=maxtmstep;
	maxvelstep();
	if (printmod) printf("\n !!! INITIAL TIME STEP FOR CYCLE %e YEARS !!! \n",timestep/3.15576e+7);
	/**/
	/**/
	///// START ITERATION WITH TN=0 guess, REFINE GUESS WITH titerate, REFEED TO viterate, cycle 4-6 times
	/**/
	/* vX,vY recalc after Stokes+Continuity equation */
	if(movemod)
		{
		if (printmod) printf("\n EPS, SIG, P, VX, VY CALCULATION...\n");
		viterate(n0);
		if (printmod) printf("EPS, SIG, P, VX, VY  OK!\n");
		}
	/**/
	/**/
	/**/
	/* Time step for marker displacement definition using FD stabilization upgrade [T. Duretz et al., Geochem. Geophys. Geosyst., 12, Q07004, 2011] */
	timestep=maxtmstep;
	maxvelstep();
	if (printmod) printf("\n !!! FINAL TIME STEP FOR CYCLE %e YEARS !!! \n",timestep/3.15576e+7);	
	/**/
	/**/
	/**/
	/* Tk recalc after Heat transport equation */
	if(timedir>0 && tempmod && timestep)
		{
		if (printmod) printf("\n TEMPERATURE CALCULATION...\n");
		titerate(n0);
		if (printmod) printf("TEMPERATURE OK!\n");
		}
	/**/
	///// --> GET BETTER ESTIMATE OF TN
	// END ITERATION
	/**/
	/**/
	/* Move marker */
	if(markmod && timestep)
		{
		/* Time step for markers definition */
		if(timedir<0 && timestep>0) timestep=-timestep;
		if (printmod) printf("\n CURRENT VALID TIME STEP %e YEARS IN CYCLE\n",timestep/3.15576e+7);
		if (printmod) printf("MARKERS MOVE...");
		movemark();
		if (printmod) printf("MARKERS OK!\n");
		}
	/**/
	/**/
	/**/
	/* Reset Vx, Vy, P values */
	if(ratemod)
		{
		if (printmod) printf("VX, VY, P RESET ...");
		/* Vx, Vy Reset Cycle */
		for (m1=0;m1<nodenum;m1++)
			{
			vx[m1]=vy[m1]=0;
			}
		/**/
		/* Pressure in cells Reset Cycle */
		for (m1=0;m1<xnumx1;m1++)
		for (m2=0;m2<ynumy1;m2++)
			{
			/* Pos in pr[] */
			m3=m1*ynumy1+m2;
			/* Recalc P */
			pr[m3]=pinit+((double)(m2)+0.5)*ystpy*GYKOEF*pkf[0];
			}
		if (printmod) printf("VX, VY, P OK!");
		}
	/**/
	/**/
	/*  Impact takes place */
	impact();
	/**/
	/* Handle pebble accretion event(s) since last timestep, Tim (2017-04-11) */
	pebbleaccr();
	/**/
	/* ro[],nu[] Recalc */
	if(gridmod)
		{
		if (printmod) printf("\n RO, NU, CP etc  RECALC AFTER NEW MARKERS POSITIONS...");
		ronurecalc();
		if (printmod) printf("RO, NU, CP etc OK!\n");
		}
	/**/
	/**/
	/**/
	/* Increase Timesum */
	timesum+=timestep;
	/* 1year=3.15576*10^7sek */
	if (printmod) printf("\n %e YEARS IN CYCLE     %e YEARS FROM START\n\n",timestep/3.15576e+7,timesum/3.15576e+7);
	/**/
	/**/
	/**/
	/* Print Results */
/*
	if(n0<cyc0max-1) saver(f0+1,n0);
*/
	/* Print Results */
	}
/**/
/**/
/**/
/* Print Results */
saver(f0+1,n0-1);
/**/
/* Exit if final time is reached */
if(timesum>timeexit) {printf("Reached timeexit (timesum %e yr > timeexit %e yr) >>> QUIT\n",timesum/3.15576e+7,timeexit/3.15576e+7); exit(0);}
/**/
/* Check core evolution (Greg addition) */
core();
/**/
/* Check grain site evolution (Greg addition) */
grain();
/**/
/* Save impact and core history */
impactsave();
/**/
/**/
/* Print Results */
/**/
/* End Output file Names Cycle */
}
/* End Program */
return 0;
}
/* Solve differential equations by MARKER-IN-CELL + FINITE-DIFFERENCES method */
