/* ADD & SOLVE MATRIX BY ECONOMICAL FRONTAL GAUSS METHOD */
int gausmat3(int am, long int mcmax, long int mcmin)
/* wn[] - line koef numbers */
/* wi[] - line koef values */
/* am - mode of action rm=1,2,3 - add, rm=0 - solve */
/* val0[] - matrix contents */
/* lin0[] - line numbers for matrix contents */
/* pos0[] - first pos numbers for line in  val0[], lin0[] */
/* num0[] - pos numbers for line in  val0[], lin0[] */
/* fre0[] - free member for lines */
/* pos0cur - first free position number in val0[], lin0[] */
/* mcmax - current line in num0[] */
/* mcmin - first line in num0[] */
{
/* Counters */
long int m1,m2,m3,nempty;
/* Limits */
long int linbeg,linend;
/* Val Buffer */
double ival;
/**/
/* Space Check */
if (mcmax>=MAXPAR)
	{
	printf("EXIT PROGRAM: Space out in fre0[] %ld",mcmax);
	exit(0);
	}
/**/
/**/
/**/
/* STEP 1: EXEED KOEF ELIMINATION MATRIX BODY FORMATION am>0 */
if (am>0)
	{
/*
printf("\n GAUS %ld %ld   %ld   %ld %e ",am,mcmax,pos0cur,wn[0],wi[0]); getchar();
*/
	/* First, Last line numbers */
	linbeg=wn[1]; linend=wn[wn[0]];
	/**/
	/* Free member reload from buffer */
	fre0[mcmax]=wi[0];
	/**/
	/* Line koef reload from buffer */
	for (m2=1;m2<=wn[0];m2++)
		{
		bufv[wn[m2]]+=wi[m2];
		/**/
		/* First, Last line numbers */
		if(wi[m2])
			{
			if (linbeg>wn[m2]) linbeg=wn[m2];
			if (linend<wn[m2]) linend=wn[m2];
			}
/*
printf("%ld %ld %e",m2,wn[m2],wi[m2]); getchar();
*/
		}
	/**/
	/* Cur line koef recalc in buffer */
	for (m2=linbeg;m2<mcmax;m2++)
		{
		if(bufv[m2])
			{
			ival=bufv[m2];
			bufv[m2]=0;
			/* Check Presence of line */
			if (!num0[m2]) 
				{
				printf("EXIT PROGRAM: Line  %ld absent in matrix when processing Line %ld",m2,mcmax);
				exit(0);
				}
			/* Current Line koef recalc after cur koef and upper line */
			/* 1-st coef of any upper line = 1 */
			for (m3=1;m3<num0[m2];m3++)
				{
				bufv[lin0[pos0[m2]+m3]]-=ival*val0[pos0[m2]+m3];
				}
			/**/
			/* Free member recalc after cur koef upper line av=0,1,2 */
			fre0[mcmax]-=fre0[m2]*ival;
			/* Check last line number */
			linend=MAXV(linend,lin0[pos0[m2]+num0[m2]-1]);
			/**/
			}
		}
	/**/
	/* Cur line save in val0[],lin0[],num0[],pos0[] */
	/* Check Singularity */
	/* Check Presence of line */
	if (!bufv[mcmax])
		{
		printf("EXIT PROGRAM: Matrix is singular at Line %ld",mcmax);
		exit(0);
		}
	pos0[mcmax]=pos0cur;
	num0[mcmax]=0;
	ival=bufv[mcmax];
	for (m2=mcmax;m2<=linend;m2++)
		{
		/**/
		/* Recalc and Save Val>0 */
		if (bufv[m2])
			{
			/* Save Cur Koef */
			lin0[pos0cur]=m2;
			val0[pos0cur]=bufv[m2]/ival;
/*
printf("%ld %ld %ld %ld %e ",mcmax,pos0cur,m2,num0[pos0cur],val0[pos0cur]); getchar();
*/
			pos0cur++;
			/* Check Space */
			if (pos0cur>=MAXMAT) 
				{
				printf("EXIT PROGRAM: Space out in val0[] %ld %ld",mcmax,pos0cur);
				exit(0);
				}
			num0[mcmax]++;
			}
		/* Clear Cur Koef */
		bufv[m2]=0;
		}
	/* Free member recalc */
	fre0[mcmax]/=ival;
	return 0;
	}
/* End STEP 1: EXEED KOEF ELIMINATION MATRIX BODY FORMATION am>0 */
/**/
/**/
/**/
/* STEP 3: SOLUTION CALC CHECK am=0 */
nempty=0;
printf("TOTAL BAND MATRIX SIZE: val0[%ld]\n",pos0cur);
/*
*/
for (m1=mcmax;m1>=mcmin;m1--)
	{
	/* Calc sol0[] */
	if (num0[m1])
		{
		/* Recalc koef after Sol */
		ival=fre0[m1];
		for (m2=0;m2<num0[m1];m2++)
			{
			ival-=val0[pos0[m1]+m2]*sol0[lin0[pos0[m1]+m2]];
			}
		/* Calc Sol */
		sol0[m1]=ival;
		}
	else
		{
		sol0[m1]=0;
/*
printf("%ld %ld %e ",m1,num0[m1],sol0[m1]); getchar();
*/
		nempty++;
		}
/*
printf("%ld %ld %e ",m1,num0[m1],sol0[m1]); getchar();
*/
	}
/* End STEP 3: SOLUTION CALC CHECK am=0 */
/*
printf("%ld %ld %ld %e",nempty,pos0cur,mcmax,(double)(pos0cur+1)/(double)(mcmax+1)); getchar();
*/
return nempty;
}
/* End SOLVE MATRIX BY ECONOMICAL GAUSS METHOD */





/* ADD/SOLVE PARDISO MATRIX */
int gausmat4(int am, long int mcmax, long int mcmin)
/* wn[] - line koef numbers */
/* wi[] - line koef values */
/* am - mode of action rm=1,2,3 - add, rm=0 - solve */
/* ia[] - first coloumn of each row */
/* ja[] - coloumn position of each koef in the respective row from left to right */
/* a[]  - value of each koef */
/* b[] - Right hand side */
/* bufv[] - koefficient values buffer */
/* cur0[] - koefficient use y/n */
{
/* Counters */
long int m1,m2,m3,leftnum,rightnum;
int i;
/* Val Buffer */
double ival;
/**/
/**/
/**/
/* Space Check */
if (mcmax>=MAXPAR)
	{
	printf("EXIT PROGRAM: Space out in fre0[] %ld",mcmax);
	exit(0);
	}
/**/
/**/
/**/
/* Solve Matrix ===================================================== */
if (am==0)
{
/* Solve Matrix with PARDISO */
/* Matrix data. */
/*Calling function making arrays */
int  na=mcmax;
/**/
int mtype = 11; /* Real unsymmetric matrix */
/* RHS and solution vectors. */
int nrhs = 1; /* Number of right hand sides. */
/* Internal solver memory pointer pt, */
/* 32-bit: int pt[64]; 64-bit: long int pt[64] */
/* or void *pt[64] should be OK on both architectures */
void *pt[64];
/* Pardiso control parameters. */
int iparm[64];
int maxfct, mnum, phase, error, msglvl;
/* Auxiliary variables. */
double ddum; /* Double dummy */
int idum; /* Integer dummy. */
/* -------------------------------------------------------------------- */
/* .. Setup Pardiso control parameters. */
/* -------------------------------------------------------------------- */
for (i = 0; i < 64; i++)
	{
	iparm[i] = 0;
	}
iparm[0] = 1; /* No solver default */
iparm[1] = 0; /* 0: minimum degree compaction, 2: Fill-in reordering from METIS */
/* Numbers of processors, value of OMP_NUM_THREADS */
iparm[2] = 1;
iparm[3] = 0; /* No iterative-direct algorithm */
iparm[4] = 0; /* No user fill-in reducing permutation */
iparm[5] = 0; /* Write solution into x */
iparm[6] = 0; /* Not in use */
iparm[7] = 20; /* Max numbers of iterative refinement steps */
iparm[8] = 0; /* Not in use */
iparm[9] = 13; /* Perturb the pivot elements with 1E-13 */
iparm[10] = 1; /* Use nonsymmetric permutation and scaling MPS */
iparm[11] = 0; /* Not in use */
iparm[12] = 0; /* Not in use */
iparm[13] = 0; /* Output: Number of perturbed pivots */
iparm[14] = 0; /* Not in use */
iparm[15] = 0; /* Not in use */
iparm[16] = 0; /* Not in use */
iparm[17] = -1; /* Output: Number of nonzeros in the factor LU */
iparm[18] = -1; /* Output: Mflops for LU factorization */
iparm[19] = 0; /* Output: Numbers of CG Iterations */
maxfct = 1; /* Maximum number of numerical factorizations. */
mnum = 1; /* Which factorization to use. */
msglvl = 1; /* Print statistical information in file */
error = 0; /* Initialize error flag */
/* -------------------------------------------------------------------- */
/* .. Initialize the internal solver memory pointer. This is only */
/* necessary for the FIRST call of the PARDISO solver. */
/* -------------------------------------------------------------------- */
for (i = 0; i < 64; i++)
	{
	pt[i] = 0;
	}
/* -------------------------------------------------------------------- */
/* .. Reordering and Symbolic Factorization. This step also allocates */
/* all memory that is necessary for the factorization. */
/* -------------------------------------------------------------------- */
phase = 11;
PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
	&na, a, ia, ja, &idum, &nrhs,
	iparm, &msglvl, &ddum, &ddum, &error);
if (error != 0) 
	{
	printf("\nERROR during symbolic factorization: %d", error);
	exit(1);
	}	
printf("\nReordering completed ... ");
printf("\nNumber of nonzeros in factors = %d", iparm[17]);
printf("\nNumber of factorization MFLOPS = %d", iparm[18]);
/* -------------------------------------------------------------------- */
/* .. Numerical factorization. */
/* -------------------------------------------------------------------- */
phase = 22;
PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
	&na, a, ia, ja, &idum, &nrhs,
	iparm, &msglvl, &ddum, &ddum, &error);
if (error != 0) 
	{
	printf("\nERROR during numerical factorization: %d", error);
	exit(2);
	}
printf("\nFactorization completed ... ");
/* -------------------------------------------------------------------- */
/* .. Back substitution and iterative refinement. */
/* -------------------------------------------------------------------- */
phase = 33;
iparm[7] = 20; /* Max numbers of iterative refinement steps. */
PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
	&na, a, ia, ja, &idum, &nrhs,
	iparm, &msglvl, b, x, &error);
if (error != 0) 	
	{
	printf("\nERROR during solution: %d", error);
	exit(3);
	}
printf("\nSolve completed ... ");
/* -------------------------------------------------------------------- */
/* .. Termination and release of memory. */
/* -------------------------------------------------------------------- */
if (1==1)
{
phase = -1; /* Release internal memory. */
PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
	&na, &ddum, ia, ja, &idum, &nrhs,
	iparm, &msglvl, &ddum, &ddum, &error);
}
if (printmod) printf("PARDISO SOLVER OK!\n");
return 0;
}
/* End Solve Matrix ===================================================== */
/**/
/**/
/**/
/**/
/**/
/* STEP 0: RELOAD KOEF FROM wi[] to val0[] */
/* Reset coefficients */
for (m2=1;m2<=wn[0];m2++)
	{
	/* wi[]>0 */
	if (wi[m2])
		{
		bufv[wn[m2]]=0.0;
		cur0[wn[m2]]=1;
		}
	}
leftnum=nodenum3+1;
rightnum=-1;
/* Sum up coefficients */
for (m2=1;m2<=wn[0];m2++)
	{
	/* wi[]>0 */
	if (wi[m2])
		{
		/* Check left/right band width */
		m3=wn[m2];
		leftnum=MINV(leftnum,m3);
		rightnum=MAXV(rightnum,m3);
		/* Add coefficient */
		bufv[wn[m2]]+=wi[m2];
		}
	}
/**/
/**/
/**/
/* Free member reload from buffer */
/*
fre0[mcmax]=wi[0];
*/
b[mcmax]=wi[0];
/**/
/* Line koef reload to a[] from buffer */
/*
pos0[mcmax]=pos0cur;
num0[mcmax]=0;
*/
ia[mcmax]=pos0cur+1;
ia[mcmax+1]=pos0cur+1;
for (m2=leftnum;m2<=rightnum;m2++)
	{
/*
	printf("\n GAUS    %ld   %ld   %d %e %ld %ld %e %d %d\n",leftnum,rightnum,wn[0],wi[0],m2,pos0cur,val1[1],lin0[1],num0[mcmax]); 
*/
	/* Marked Line */
	if(cur0[m2]==1)
		{
		/* wi[]>0 */
		if (bufv[m2])
			{
			/* Save Cur Koef */
/*
			lin0[pos0cur]=m2;
			val1[pos0cur]=bufv[m2];
*/
			ja[pos0cur]=m2+1;
			a[pos0cur]=bufv[m2];
			pos0cur++;
			/* Check Space */
			if (pos0cur>=MAXVAL) 
				{
				printf("EXIT PROGRAM: Space out in a[] %ld %ld",mcmax,pos0cur);
				exit(0);
				}
/*
			num0[mcmax]++;
*/
			ia[mcmax+1]++;
			}
		/* Clean marked line */
		bufv[m2]=0.0;
		cur0[m2]=0;
		}
	}
return 0;
}
/* End ADD/SOLVE PARDISO MATRIX */



