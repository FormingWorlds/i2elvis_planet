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
/* --------------------- INCLUDE PARTS --------------------- */



/* Formation of Data file for i2.c program */
int main()
{
/* Counters */
int n1,n2,n3,n4,n5,yn;
long int m1,m2,m3,m4,m5,m6,m7;
long int mm1,mm2,minmx=0,maxmx=0,minmy=0,maxmy=0,mm3;
/* Nonstability for markers, m */
int nonstab;
double pival,wlen,wamp,wymin,wymid,wymax,xnonstab=0,ynonstab=0;
/* Distance val */
double x,y,x1,y1,x2,y2,x3,y3,x4,y4,dx,dy;
/* Initial Rock type for Grid,Rock bodyes Num, Initial Distribution */
double bx[4],by[4];
double cnst,koef,koef1,koef2,ival,ival1,akf,bkf;
/* T K in Edges of Grid */
double t[10];
long int m10,m11,m20,m21,nshift,nshift1,nshift2;
/* Set Pi value */
pival=2.000*asin(1.000);
/**/
/*  Read impact history (Gregor addition) */
impactread();
/**/
/**/
/* Load configuration from mode.t3c */
loadconf();
/**/
/**/
/**/
/* Open File init.t3c */
fl1 = fopen("init.t3c","rt");
printf("Formation of Initial Condition ...\n");
/**/
/**/
/**/
/* Grid Parameters */
ffscanf1();xnumx=atoi(sa);
ffscanf1();ynumy=atoi(sa);
ffscanf1();mnumx=atoi(sa);
if (mnumx<0)
	{
	mnumx=-mnumx;
	ffscanf1();minmx=atoi(sa);
	ffscanf1();maxmx=atoi(sa);
	}
ffscanf1();mnumy=atoi(sa);
if (mnumy<0)
	{
	mnumy=-mnumy;
	ffscanf1();minmy=atoi(sa);
	ffscanf1();maxmy=atoi(sa);
	}
ffscanf1();xsize=atof(sa);
ffscanf1();ysize=atof(sa);
ffscanf1();gamma_eff=atof(sa);
ffscanf1();memory_fe=atof(sa);
ffscanf1();memory_si=atof(sa);
ffscanf1();por_init=atof(sa);
ffscanf1();gr_init=atof(sa);
ffscanf1();znumz=atoi(sa);
ffscanf1();corr2d3d=atoi(sa);
ffscanf1();pinit=atof(sa);
ffscanf1();GXKOEF=atof(sa);
ffscanf1();GYKOEF=atof(sa);
ffscanf1();al2627_init=atof(sa)*1.0e-5;
ffscanf1();fe6056_init=atof(sa)*1.0e-8;
ffscanf1();timesum=atof(sa)*3.15576e+7;
ffscanf1();nonstab=atoi(sa);
/**/
printf("26Al/27Al and 60Fe/56Fe ratios at time 0: %e, %e \n",al2627_init,fe6056_init);
/* Regular Nonstability Read */
if(nonstab<0)
	{
	ffscanf1();wamp=atof(sa);
	ffscanf1();wlen=atof(sa);
	if (nonstab==-2)
		{
		ffscanf1();wymin=atof(sa);
		ffscanf1();wymid=atof(sa);
		ffscanf1();wymax=atof(sa);
		}
	}
/**/
/* Random Nonstability Read */
if(nonstab>0)
	{
	ffscanf1();xnonstab=atof(sa)*xsize/((double)(xnumx-1))/((double)(mnumx));
	ffscanf1();ynonstab=atof(sa)*ysize/((double)(ynumy-1))/((double)(mnumy));
	}
/**/
/**/
/**/
/* First Calc,Check Grid parameters */
gridcheck();
/**/
/**/
/**/
/* Set initial gridlines positions */
gx[0]=0;
for (m1=1;m1<xnumx;m1++) gx[m1]=gx[m1-1]+xstpx;
gy[0]=0;
for (m2=1;m2<ynumy;m2++) gy[m2]=gy[m2-1]+ystpy;
/**/
/**/
/**/
/* Load Markers Types,X,Y from File name */
ffscanf1();
if (sa[0]!='0')
	{
	/* Reload File Name */
	for (n1=0;n1<15;n1++) fl1in[n1]=sa[n1];
	ffscanf1();fl1itp=0; if(sa[0] == 'b') fl1itp=1;
	/* Load data from input file */
	loader();
	}
/* Marker X,Y  set */
else
	{
	/* Marker Num calc */
	marknum=mnumx*mnumy*(xnumx-1+maxmx-minmx)*(ynumy-1+maxmy-minmy);
	printf("Set markers location \n");
	dx=xstpx/(double)(mnumx);
	dy=ystpy/(double)(mnumy);
	mm1=0;
	/* Marker Types add */
	for (m1=minmx;m1<xnumx-1+maxmx;m1++)
	for (m2=minmy;m2<ynumy-1+maxmy;m2++)
		{
		x=(double)(m1)*xstpx+dx/2.0;
		for (m4=0;m4<mnumx;m4++)
			{
			y=(double)(m2)*ystpy+dy/2.0;
			for (m5=0;m5<mnumy;m5++)
				{
				/* X,Y */
				markx[mm1]=(float)(x);
				marky[mm1]=(float)(y);
				/* Add Y */
				y+=dy;
				/* Add mark Num */
				mm1++;
				}
			/* Add X */
			x+=dx;
			}
		}
	}
/**/
/**/
/**/
/* Random Nonstability Set on markers */
if(nonstab>0)
	{
	printf("Set Random nonstability on markers\n");
	for (mm1=0;mm1<marknum;mm1++)
	if(markt[mm1]<100)
		{
		/* Random nonstability Set on X,Y */
		markx[mm1]+=(float)(rand() % (nonstab*2+1) - nonstab)/((float)(nonstab))*(float)(xnonstab);
		marky[mm1]+=(float)(rand() % (nonstab*2+1) - nonstab)/((float)(nonstab))*(float)(ynonstab);
		}
	}
/**/
/**/
/**/
/* Output File name */
ffscanf1();
for (n1=0;n1<15;n1++) fl1out[n1]=sa[n1];
ffscanf1(); if(sa[0] == 'b') fl1otp=1;
/**/
/**/
/**/
/* Set initial Surface Level */
/*
if(erosmod) for (m1=0;m1<xnumx;m1++) ep[m1]=ep[xnumx+m1]=-ystpy;
*/
/**/
/**/
/**/
/* Rock Properties Read */
printf("Read rocks properties \n");
ffscanf1();
while (sa[0]!='~')
	{
	/* Immobile markers Y/N */
	if(sa[0]=='i') 
		{
		n1=atoi(sa+1);
		markim[n1]=1;
		}
	else
		{
		n1=atoi(sa);
		markim[n1]=0;
		}
	ffscanf1();markn0[n1]=atof(sa);
	ffscanf1();markn1[n1]=atof(sa);
	ffscanf1();marks0[n1]=atof(sa);
	ffscanf1();marks1[n1]=atof(sa);
	ffscanf1();marknu[n1]=atof(sa);
	ffscanf1();markdh[n1]=atof(sa);
	ffscanf1();markdv[n1]=atof(sa);
	ffscanf1();markss[n1]=atof(sa);
	ffscanf1();markmm[n1]=atof(sa);
	ffscanf1();markll[n1]=atof(sa);
	ffscanf1();marka0[n1]=atof(sa);
	ffscanf1();marka1[n1]=atof(sa);
	ffscanf1();markb0[n1]=atof(sa);
	ffscanf1();markb1[n1]=atof(sa);
	ffscanf1();marke0[n1]=atof(sa);
	ffscanf1();marke1[n1]=atof(sa);
	ffscanf1();markro[n1]=atof(sa);
	ffscanf1();markbb[n1]=atof(sa);
	ffscanf1();markaa[n1]=atof(sa);
	ffscanf1();markcp[n1]=atof(sa);
	ffscanf1();markkt[n1]=atof(sa);
	ffscanf1();markkf[n1]=atof(sa);
	ffscanf1();markkp[n1]=atof(sa);
	/* Ht in Wt/kg load Y(k)/N() */
	ffscanf1(); if(sa[0]=='k') markht[n1]=atof(sa+1)*markro[n1]; else markht[n1]=atof(sa);
	ffscanf1();
	rocknum=MAXV(rocknum,n1+1);
	}
/**/
/**/
/**/
/* 1-st ro[],nu[], cp[] kt[] for nodes Recalc after markers */
/*
printf("1-st Recalc nodes properties after markers types \n");
ronurecalc();
*/
/**/
/* Second Calc,Check Grid parameters */
gridcheck();
/**/
/**/
/**/
/* Bondary Conditions reset */
for (m1=0;m1<nodenum;m1++)
	{
	/* 0   1   2  */
	/* P,  Vx, Vy */
	m3=m1*3;
	for (m2=0;m2<3;m2++)
		{
		bondm[m3+m2]=0;
		}
	/* T */
	bondm[nodenum3+m1]=0;
	/* Nu */
	mu[m1]=0;
	}
/**/
/* Bondary Conditions Read,Set to Grid */
printf("Set boundary conditions\n");
bondnum=1;
ffscanf1();
while (sa[0]!='~')
	{
	/* Type of bondary geometry: */
	/* 0 - regular box */
	/* 1 - irregular rectangle */
	/* 2 - circular sector */
	yn=0;
	if(sa[0] == 'R') yn=1;
	if(sa[0] == 'C') yn=2;
	/**/
	/* VAR for bondary conditions: Vx,Vy,T,P,Nu  */
	n1=-1;
	if(sa[1] == 'x') n1=0;
	if(sa[1] == 'y') n1=1;
	if(sa[0] == 'P' || sa[1] == 'P') n1=2;
	if(sa[0] == 'T' || sa[1] == 'T') n1=3;
	if(sa[1] == 'u') n1=4;
	if(sa[0] == 'E' || sa[1] == 'E') n1=5;
	if(sa[0] == 'H' || sa[1] == 'H') n1=6;
	if(sa[0] == 'L' || sa[1] == 'L') n1=7;
	if(sa[0] == 'X' || sa[1] == 'X') n1=8;
	if(sa[0] == 'Y' || sa[1] == 'Y') n1=9;
	if(sa[0] == 'M' || sa[1] == 'M') {n1=10; printf("Set new Cell markers location \n");}
	if(n1 == -1) {printf("Unknown Parameter <%s>",sa); exit(0);}
	/**/
	/**/
	/**/
	/* Bond Parameters Read */
	if(yn==0)
		{
		/* Bond Parameters Read for cur regular Box */
		/* VAR___m10___m11___m20___m21___Const___Koef___dm1___dm2 */
		/* Vx    1     x-1   0     0     0       1.0    0     +1 */
		/* Vx    1     x-1   y     y     0       0 */
		ffscanf1();if (sa[0] == 'x') m10=xnumx-1+atoi(sa+1); else m10=atoi(sa);
		ffscanf1();if (sa[0] == 'x') m11=xnumx-1+atoi(sa+1); else m11=atoi(sa);
		ffscanf1();if (sa[0] == 'y') m20=ynumy-1+atoi(sa+1); else m20=atoi(sa);
		ffscanf1();if (sa[0] == 'y') m21=ynumy-1+atoi(sa+1); else m21=atoi(sa);
		}
	else
		{
		/* Read parameters for cur irreg box, circle */
		/* 0  2 */
		/* 1  3 */
		for(n2=0;n2<4;n2++)
			{
			ffscanf1(); if(sa[0]=='m') bx[n2]=atof(sa+1)/xsize; else bx[n2]=atof(sa);
			ffscanf1(); if(sa[0]=='m') by[n2]=atof(sa+1)/ysize; else by[n2]=atof(sa);
			}
		m10=0;
		m11=xnumx-1;
		m20=0;
		m21=ynumy-1;
		}
	/**/
	/**/
	/**/
	/* Read constant */
	ffscanf1();cnst=atof(sa);
	/**/
	/**/
	/**/
	/* Read 1st koef */
	ffscanf1();if (sa[0] == 'c') koef=1.1e+30; else koef=atof(sa);
	/* printf("%d %ld %ld %ld %ld %e %e %e",n1,m10,m11,m20,m21,cnst,koef,ival);getchar(); */
	/* Read 1st shift parameters */
	n2=0;
	if(koef)
		{
		ffscanf1();nshift=atoi(sa)*ynumy;
		/* Vy=Vx boundary condition */
		if (sa[0] == 'x')
			{
			n2=-1;
			nshift=atoi(sa+1)*ynumy;
			}
		/* Vx=Vy boundary condition */
		if (sa[0] == 'y')
			{
			n2=1;
			nshift=atoi(sa+1)*ynumy;
			}
		ffscanf1();nshift+=atoi(sa);
		}
	/**/
	/**/
	/**/
	/* Read 2nd koef */
	ffscanf1(); koef1=atof(sa);
	/* Read 2nd shift parameters */
	n3=0;
	if(koef1)
		{
		ffscanf1();nshift1=atoi(sa)*ynumy;
		/* Vy=Vx boundary condition */
		if (sa[0] == 'x')
			{
			n3=-1;
			nshift1=atoi(sa+1)*ynumy;
			}
		/* Vx=Vy boundary condition */
		if (sa[0] == 'y')
			{
			n3=1;
			nshift1=atoi(sa+1)*ynumy;
			}
		ffscanf1();nshift1+=atoi(sa);
		}
	/**/
	/**/
	/**/
	/* Read 3rd koef */
	ffscanf1(); koef2=atof(sa);
	/* Read 2nd shift parameters */
	n4=0;
	if(koef2)
		{
		ffscanf1();nshift2=atoi(sa)*ynumy;
		/* Vy=Vx boundary condition */
		if (sa[0] == 'x')
			{
			n4=-1;
			nshift2=atoi(sa+1)*ynumy;
			}
		/* Vx=Vy boundary condition */
		if (sa[0] == 'y')
			{
			n4=1;
			nshift2=atoi(sa+1)*ynumy;
			}
		ffscanf1();nshift2+=atoi(sa);
		}
	/**/
	/**/
	/**/
	/* Next VAR for bondary conditions read */
	ffscanf1();
	/**/
	/**/
	/**/
	/* Conditions set for Grid */
	/* m10,m11 - X min,max Num */
	/* m20,m21 - Y min,max Num */
	/* n1 - Variable: Vx(0), Vy(1), T(2) */
	/* nshift - Shift for Symmetry conditions */
	/**/
	for (m1=m10;m1<=m11;m1++)
	for (m2=m20;m2<=m21;m2++)
		{
		/* Node num */
		m4=m1*ynumy+m2;
		/* Node X,Y General coordinats calc */
		x=gx[m1]/xsize;
		y=gy[m2]/ysize;
		/**/
		switch(n1)
			{
			/* Vx Bondary equations */
			case 0:
			/* Y General coordinat recalc */
			y=(gy[m2]+gy[m2+1])/2.0/ysize;
			/*  Vx Parameter number */
			m5=m4*3+1;
			/* Regular box */
			if (yn==0)
				{
				bondm[m5]=bondnum;bondv[bondnum][0]=vx[m4]=cnst;
				bondn[bondnum][0]=bondn[bondnum][1]=bondn[bondnum][2]=0;
				bondv[bondnum][1]=bondv[bondnum][2]=bondv[bondnum][3]=0;
				if (koef)
					{
					bondn[bondnum][0]=(m5)+(nshift)*3+n2+ 1;bondv[bondnum][1]=koef;
					}
				if (koef1)
					{
					bondn[bondnum][1]=(m5)+(nshift1)*3+n3+ 1;bondv[bondnum][2]=koef1;
					}
				if (koef2)
					{
					bondn[bondnum][2]=(m5)+(nshift2)*3+n4+ 1;bondv[bondnum][3]=koef2;
					}
				bondnum++;
				}
			/**/
			/* Irregular rectangle */
			if (yn==1)
				{
				/* Node X General coordinat check */
				x1=bx[0]+(y-by[0])*(bx[1]-bx[0])/(by[1]-by[0]);
				x2=bx[2]+(y-by[2])*(bx[3]-bx[2])/(by[3]-by[2]);
				if(x>=x1 && x<=x2)
					{
					/* Node Y limit box coord calc */
					y1=by[0]+(x-bx[0])*(by[2]-by[0])/(bx[2]-bx[0]);
					y2=by[1]+(x-bx[1])*(by[3]-by[1])/(bx[3]-bx[1]);
					/* Node Y general coordinat check */
					if(y>=y1 && y<=y2)
						{
						bondm[m5]=bondnum;bondv[bondnum][0]=vx[m4]=cnst;
						bondn[bondnum][0]=bondn[bondnum][1]=bondn[bondnum][2]=0;
						bondv[bondnum][1]=bondv[bondnum][2]=bondv[bondnum][3]=0;
						if (koef)
							{
							bondn[bondnum][0]=(m5)+(nshift)*3+n2+ 1;bondv[bondnum][1]=koef;
							}
						if (koef1)
							{
							bondn[bondnum][1]=(m5)+(nshift1)*3+n3+ 1;bondv[bondnum][2]=koef1;
							}
						if (koef2)
							{
							bondn[bondnum][2]=(m5)+(nshift2)*3+n4+ 1;bondv[bondnum][3]=koef2;
							}
						bondnum++;
						}
					}
				}
			/**/
			/* Circular sector */
			if (yn==2)
				{
				/* 0 X1,Y1  2 R1,R2 */
				/* 1 X2,Y2  3 L1,L2 */
				/* Calculate and check distance from first center */
				x1=x-bx[0];
				y1=(y-by[0])*ysize/xsize;
				x4=bx[2]-pow(x1*x1+y1*y1,0.5);
				if(x4>=0)
					{
					/* Calculate and check distance from second center */
					x2=x-bx[1];
					y2=(y-by[1])*ysize/xsize;
					y4=pow(x2*x2+y2*y2,0.5)-by[2]*ysize/xsize;
					if(y4>=0)
						{
						/* Calculate and check angle relatively first or second center */
						x3=x1;y3=y1; if(by[3]<0) {x3=x2; y3=y2;} 
						ival=pow(x3*x3+y3*y3,0.5); if(!ival) ival=1.0;
						ival=x3/ival; ival=asin(ABSV(ival))/pival*180.0;
						if(x3>=0 && y3>0) ival=180.0-ival;
						if(x3<0 && y3>0) ival+=180.0;
						if(x3<0 && y3<=0) ival=360.0-ival;
						if(ival>=ABSV(bx[3]) && ival<=ABSV(by[3]))
							{
							/* Relative coordinats calc */
							x4+=y4;if(!x4) x4=1.0;
							x=y4/x4;
							y=cos(ival/180.0*pival);
							/* Vx recalc after relativ coord */
							bondn[bondnum][0]=bondn[bondnum][1]=bondn[bondnum][2]=0;
							bondv[bondnum][1]=bondv[bondnum][2]=bondv[bondnum][3]=0;
							bondm[m5]=bondnum;bondv[bondnum][0]=vx[m4]=y*(cnst*x+koef*(1.0-x));
							bondnum++;
							}
						}
					}
				}
			break;
			/**/
			/**/
			/**/
			/* Vy Bondary equations */
			case 1:
			/* X General coordinat recalc */
			x=(gx[m1]+gx[m1+1])/2.0/xsize;
			/*  Vx Parameter number */
			m5=m4*3+2;
			/* Regular box */
			if (yn==0)
				{
				bondm[m5]=bondnum;bondv[bondnum][0]=vy[m4]=cnst;
				bondn[bondnum][0]=bondn[bondnum][1]=bondn[bondnum][2]=0;
				bondv[bondnum][1]=bondv[bondnum][2]=bondv[bondnum][3]=0;
				if (koef)
					{
					bondn[bondnum][0]=(m5)+(nshift)*3+n2+ 1;bondv[bondnum][1]=koef;
					}
				if (koef1)
					{
					bondn[bondnum][1]=(m5)+(nshift1)*3+n3+ 1;bondv[bondnum][2]=koef1;
					}
				if (koef2)
					{
					bondn[bondnum][2]=(m5)+(nshift2)*3+n4+ 1;bondv[bondnum][3]=koef2;
					}
				bondnum++;
				}
			/**/
			/* Irregular rectangle */
			if (yn==1)
				{
				/* Node X General coordinat check */
				x1=bx[0]+(y-by[0])*(bx[1]-bx[0])/(by[1]-by[0]);
				x2=bx[2]+(y-by[2])*(bx[3]-bx[2])/(by[3]-by[2]);
				if(x>=x1 && x<=x2)
					{
					/* Node Y limit box coord calc */
					y1=by[0]+(x-bx[0])*(by[2]-by[0])/(bx[2]-bx[0]);
					y2=by[1]+(x-bx[1])*(by[3]-by[1])/(bx[3]-bx[1]);
					/* Node Y general coordinat check */
					if(y>=y1 && y<=y2)
						{
						bondm[m5]=bondnum;bondv[bondnum][0]=vy[m4]=cnst;
						bondn[bondnum][0]=bondn[bondnum][1]=bondn[bondnum][2]=0;
						bondv[bondnum][1]=bondv[bondnum][2]=bondv[bondnum][3]=0;
						if (koef)
							{
							bondn[bondnum][0]=(m5)+(nshift)*3+n2+ 1;bondv[bondnum][1]=koef;
							}
						if (koef1)
							{
							bondn[bondnum][1]=(m5)+(nshift1)*3+n3+ 1;bondv[bondnum][2]=koef1;
							}
						if (koef2)
							{
							bondn[bondnum][2]=(m5)+(nshift2)*3+n4+ 1;bondv[bondnum][3]=koef2;
							}
						bondnum++;
						}
					}
				}
			/**/
			/* Circular sector */
			if (yn==2)
				{
				/* 0 X1,Y1  2 R1,R2 */
				/* 1 X2,Y2  3 L1,L2 */
				/* Calculate and check distance from first center */
				x1=x-bx[0];
				y1=(y-by[0])*ysize/xsize;
				x4=bx[2]-pow(x1*x1+y1*y1,0.5);
				if(x4>=0)
					{
					/* Calculate and check distance from second center */
					x2=x-bx[1];
					y2=(y-by[1])*ysize/xsize;
					y4=pow(x2*x2+y2*y2,0.5)-by[2]*ysize/xsize;
					if(y4>=0)
						{
						/* Calculate and check angle relatively first or second center */
						x3=x1;y3=y1; if(by[3]<0) {x3=x2; y3=y2;} 
						ival=pow(x3*x3+y3*y3,0.5); if(!ival) ival=1.0;
						ival=x3/ival; ival=asin(ABSV(ival))/pival*180.0;
						if(x3>=0 && y3>0) ival=180.0-ival;
						if(x3<0 && y3>0) ival+=180.0;
						if(x3<0 && y3<=0) ival=360.0-ival;
						if(ival>=ABSV(bx[3]) && ival<=ABSV(by[3]))
							{
							/* Relative coordinats calc */
							x4+=y4;if(!x4) x4=1.0;
							x=y4/x4;
							y=sin(ival/180.0*pival);
							/* Vy recalc after relativ coord */
							bondn[bondnum][0]=bondn[bondnum][1]=bondn[bondnum][2]=0;
							bondv[bondnum][1]=bondv[bondnum][2]=bondv[bondnum][3]=0;
							bondm[m5]=bondnum;bondv[bondnum][0]=vy[m4]=y*(cnst*x+koef*(1.0-x));
							bondnum++;
							}
						}
					}
				}
			break;
			/**/
			/**/
			/**/
			/* P Bondary equations */
			case 2:
			/* X,Y General coordinat recalc */
			x=(gx[m1]+gx[m1+1])/2.0/xsize;
			y=(gy[m2]+gy[m2+1])/2.0/ysize;
			/*  P Parameter number */
			m5=m4*3+0;
			/* Regular box */
			if (yn==0)
				{
				bondm[m5]=bondnum;bondv[bondnum][0]=cnst;
				bondn[bondnum][0]=bondn[bondnum][1]=bondn[bondnum][2]=0;
				bondv[bondnum][1]=bondv[bondnum][2]=bondv[bondnum][3]=0;
				if (koef)
					{
					bondn[bondnum][0]=(m5)+(nshift)*3+ 1;bondv[bondnum][1]=koef;
					}
				if (koef1)
					{
					bondn[bondnum][1]=(m5)+(nshift1)*3+1;bondv[bondnum][2]=koef1;
					}
				if (koef2)
					{
					bondn[bondnum][2]=(m5)+(nshift2)*3+1;bondv[bondnum][3]=koef2;
					}
				bondnum++;
				}
			/**/
			/* Irregular rectangle */
			if (yn==1)
				{
				/* Node X General coordinat check */
				x1=bx[0]+(y-by[0])*(bx[1]-bx[0])/(by[1]-by[0]);
				x2=bx[2]+(y-by[2])*(bx[3]-bx[2])/(by[3]-by[2]);
				if(x>=x1 && x<=x2)
					{
					/* Node Y limit box coord calc */
					y1=by[0]+(x-bx[0])*(by[2]-by[0])/(bx[2]-bx[0]);
					y2=by[1]+(x-bx[1])*(by[3]-by[1])/(bx[3]-bx[1]);
					/* Node Y general coordinat check */
					if(y>=y1 && y<=y2)
						{
						bondm[m5]=bondnum;bondv[bondnum][0]=cnst;
						bondn[bondnum][0]=bondn[bondnum][1]=bondn[bondnum][2]=0;
						bondv[bondnum][1]=bondv[bondnum][2]=bondv[bondnum][3]=0;
						if (koef)
							{
							bondn[bondnum][0]=(m5)+(nshift)*3+ 1;bondv[bondnum][1]=koef;
							}
						if (koef1)
							{
							bondn[bondnum][1]=(m5)+(nshift1)*3+1;bondv[bondnum][2]=koef1;
							}
						if (koef2)
							{
							bondn[bondnum][2]=(m5)+(nshift2)*3+1;bondv[bondnum][3]=koef2;
							}
						bondnum++;
						}
					}
				}
			break;
			/**/
			/**/
			/**/
			/* T bondary equations */
			case 3:
			m5=nodenum3+m4;
			/* Regular box */
			if (yn==0)
				{
				bondm[m5]=bondnum;
				bondn[bondnum][0]=bondn[bondnum][1]=bondn[bondnum][2]=0;
				bondv[bondnum][1]=bondv[bondnum][2]=bondv[bondnum][3]=0;
				if(cnst<=0 && !koef && !koef1 && !koef2)
					{
					bondv[bondnum][0]=tk[m4];
					}
				else
					{
					bondv[bondnum][0]=cnst;
					}
				if (koef)
					{
					bondn[bondnum][0]=(m5)+(nshift)+ 1;bondv[bondnum][1]=koef;
					}
				if (koef1)	
					{
					bondn[bondnum][1]=(m5)+(nshift1)+ 1;bondv[bondnum][2]=koef1;
					}
				if (koef2)
					{
					bondn[bondnum][2]=(m5)+(nshift2)+ 1;bondv[bondnum][3]=koef2;
					}
				bondnum++;
				}
			/**/
			/* Irregular rectangle */
			if (yn==1)
				{
				/* Node X General coordinat check */
				x1=bx[0]+(y-by[0])*(bx[1]-bx[0])/(by[1]-by[0]);
				x2=bx[2]+(y-by[2])*(bx[3]-bx[2])/(by[3]-by[2]);
				if(x>=x1 && x<=x2)
					{
					/* Node Y limit box coord calc */
					y1=by[0]+(x-bx[0])*(by[2]-by[0])/(bx[2]-bx[0]);
					y2=by[1]+(x-bx[1])*(by[3]-by[1])/(bx[3]-bx[1]);
					/* Node Y general coordinat check */
					if(y>=y1 && y<=y2)
						{
						bondm[m5]=bondnum;
						bondn[bondnum][0]=bondn[bondnum][1]=bondn[bondnum][2]=0;
						bondv[bondnum][1]=bondv[bondnum][2]=bondv[bondnum][3]=0;
						if(cnst<=0 && !koef && !koef1 && !koef2)
							{
							bondv[bondnum][0]=tk[m4];
							}
						else
							{
							bondv[bondnum][0]=cnst;
							}
						if (koef)
							{
							bondn[bondnum][0]=(m5)+(nshift)+ 1;bondv[bondnum][1]=koef;
							}
						if (koef1)	
							{
							bondn[bondnum][1]=(m5)+(nshift1)+ 1;bondv[bondnum][2]=koef1;
							}
						if (koef2)
							{
							bondn[bondnum][2]=(m5)+(nshift2)+ 1;bondv[bondnum][3]=koef2;
							}
						bondnum++;
						}
					}
				}
			/**/
			break;
			/**/
			/**/
			/**/
			/* Nu Bondary equations */
			case 4:
			/* Regular box */
			if (yn==0)
				{
				mu[m4]=cnst;
				}
			/**/
			/* Irregular rectangle */
			if (yn==1)
				{
				/* Node X General coordinat check */
				x1=bx[0]+(y-by[0])*(bx[1]-bx[0])/(by[1]-by[0]);
				x2=bx[2]+(y-by[2])*(bx[3]-bx[2])/(by[3]-by[2]);
				if(x>=x1 && x<=x2)
					{
					/* Node Y limit box coord calc */
					y1=by[0]+(x-bx[0])*(by[2]-by[0])/(bx[2]-bx[0]);
					y2=by[1]+(x-bx[1])*(by[3]-by[1])/(bx[3]-bx[1]);
					/* Node Y general coordinat check */
					if(y>=y1 && y<=y2)
						{
						mu[m4]=cnst;
						}
					}
				}
			/**/
			/* Circular sector */
			if (yn==2)
				{
				/* 0 X1,Y1  2 R1,R2 */
				/* 1 X2,Y2  3 L1,L2 */
				/* Calculate and check distance from first center */
				x1=x-bx[0];
				y1=(y-by[0])*ysize/xsize;
				x4=bx[2]-pow(x1*x1+y1*y1,0.5);
				if(x4>=0)
					{
					/* Calculate and check distance from second center */
					x2=x-bx[1];
					y2=(y-by[1])*ysize/xsize;
					y4=pow(x2*x2+y2*y2,0.5)-by[2]*ysize/xsize;
					if(y4>=0)
						{
						/* Calculate and check angle relatively first or second center */
						x3=x1;y3=y1; if(by[3]<0) {x3=x2; y3=y2;} 
						ival=pow(x3*x3+y3*y3,0.5); if(!ival) ival=1.0;
						ival=x3/ival; ival=asin(ABSV(ival))/pival*180.0;
						if(x3>=0 && y3>0) ival=180.0-ival;
						if(x3<0 && y3>0) ival+=180.0;
						if(x3<0 && y3<=0) ival=360.0-ival;
						if(ival>=ABSV(bx[3]) && ival<=ABSV(by[3]))
							{
							mu[m4]=cnst;
							}
						}
					}
				}
			break;
			/**/
			/* Erosion surface definition */
			case 5:
			ival=gx[m11]-gx[m10];
			if(!ival) ival=1.0;
			ival1=gx[m1]-gx[m10];
			ep[m1]=ep[xnumx+m1]=cnst+(koef-cnst)*ival1/ival;
/*
printf("E  %e %e %e %e %e %e", ival,ival1,cnst,koef,ep[m1],ep[xnumx+m1]); getchar();
*/
			break;
			/**/
			/* Hydration Upper surface definition */
			case 6:
			ival=gx[m11]-gx[m10];
			if(!ival) ival=1.0;
			ival1=gx[m1]-gx[m10];
			ep[xnumx*2+m1]=cnst+(koef-cnst)*ival1/ival;
			if(ep[xnumx*2+m1]>ep[xnumx*3+m1]) ep[xnumx*2+m1]=ep[xnumx*3+m1];
/*
printf("H  %e %e %e %e %e", ival,ival1,cnst,koef,ep[xnumx*2+m1]); getchar();
*/
			break;
			/**/
			/**/
			/**/
			/* Hydration Lower surface definition */
			case 7:
			ival=gx[m11]-gx[m10];
			if(!ival) ival=1.0;
			ival1=gx[m1]-gx[m10];
			ep[xnumx*3+m1]=cnst+(koef-cnst)*ival1/ival;
/*
printf("L  %e %e %e %e %e", ival,ival1,cnst,koef,ep[xnumx*3+m1]); getchar();
*/
			break;
			/**/
			/**/
			/**/
			/* X coordinates of gridlines definition */
			case 8:
			ival=(double)(m11-m10);
			if(!ival) ival=1.0;
			ival1=(double)(m1-m10);
			gx[m1]=gx[m1-1]+cnst+(koef-cnst)*ival1/ival;
			if(koef1) gx[m1]=gx[m1-1]+exp(log(cnst)+log(koef1/cnst)*ival1/ival);
/*
printf("X  %e %e %e %e %e %e", ival,ival1,cnst,koef,gx[m1-1],gx[m1]); getchar();
*/
			break;
			/**/
			/**/
			/**/
			/* Y coordinates of gridlines definition */
			case 9:
			ival=(double)(m21-m20);
			if(!ival) ival=1.0;
			ival1=(double)(m2-m20);
			gy[m2]=gy[m2-1]+cnst+(koef-cnst)*ival1/ival;
			if(koef1) gy[m2]=gy[m2-1]+exp(log(cnst)+log(koef1/cnst)*ival1/ival);
/*
printf("Y  %e %e %e %e %e %e", ival,ival1,cnst,koef,gy[m2-1],gy[m2]); getchar();
*/
			break;
			/**/
			/**/
			/**/
			/* MARKER grid set to cells */
			case 10:
			dx=(gx[m1+1]-gx[m1])/(double)(nshift);
			dy=(gy[m2+1]-gy[m2])/(double)(nshift1);
			/* Marker Types add */
			x=gx[m1]+dx/2.0;
			for (mm1=0;mm1<nshift;mm1++)
				{
				y=gy[m2]+dy/2.0;
				for (mm2=0;mm2<nshift1;mm2++)
					{
					/* X,Y */
					markx[marknum]=(float)(x);
					marky[marknum]=(float)(y);
					markt[marknum]=0;
					/* Random Nonstability Set on markers */
					if(cnst>0 && koef>0)
						{
						/* Random nonstability Set on X */
						markx[marknum]+=(float)(rand() % ((int)(cnst)*2+1) - (int)(cnst))/((float)(cnst))*(float)(dx*koef);
/*
printf("X  %ld %ld %ld %ld %e %e %e ", m1,m2,mm1,mm2,dx,dy,(float)(rand() % ((int)(cnst)*2+1) - (int)(cnst))/((float)(cnst))*(float)(dx*koef));getchar();
*/
						}
					if(cnst>0 && koef1>0)
						{
						/* Random nonstability Set on X */
						marky[marknum]+=(float)(rand() % ((int)(cnst)*2+1) - (int)(cnst))/((float)(cnst))*(float)(dy*koef1);
/*
printf("Y  %ld %ld %ld %ld %e %e %e ", m1,m2,mm1,mm2,dx,dy,(float)(rand() % ((int)(cnst)*2+1) - (int)(cnst))/((float)(cnst))*(float)(dy*koef1));getchar();
printf("%ld %ld   %ld %ld   %e %e   %e %e ", m1,m2,mm1,mm2,dx,dy,markx[marknum],marky[marknum]);getchar();
*/
						}
					/* Add Y */
					y+=dy;
					/* Add mark Num */
					marknum++;
					}
				/* Add X */
				x+=dx;
				}
			break;
			/**/
			/**/
			/**/
			}
		}
	}
/* Set boundary Conditions for Grid ------------------ */
/**/
/**/
/**/
/* Reset First/Last gridline positions */
gx[0]=0;
gx[xnumx-1]=xsize;
gy[0]=0;
gy[ynumy-1]=ysize;
/**/
/**/
/**/
/* Rock bodies set on markers */
printf("Set rock bodies as marker types \n");
ffscanf1();
while (sa[0]!='~')
	{
	/* Read Type of current marker pattern: */
	/* 0 - simple rectangle */
	/* 1 - simple additional rectangle */
	/* 2 - simple rectangle with changing of markers type */
	/* 3 - simple circular pattern */
	/* 4 - simple circular pattern with changing of markers type */
	yn=atoi(sa);
	/* Immobile (i) and aditional (a)  marker notes, Rock type number or rock type increment */
	ffscanf1();
	n5=0; n1=atoi(sa);
	if(sa[0]=='i' && sa[1]!='a') {n1=atoi(sa+1); n5=1;}
	if(sa[0]=='a' && sa[1]!='i') {n1=atoi(sa+1); n5=2;}
	if((sa[0]=='i' && sa[1]=='a') || (sa[0]=='a' && sa[1]=='i')) {n1=atoi(sa+2); n5=3;}
/*
printf("%d %d %d",yn,n1,n5); getchar();
*/
	/**/
	/* Read coordinates of Cur boxe */
	/* 0  2 */
	/* 1  3 */
	for(n2=0;n2<4;n2++)
		{
		ffscanf1(); if(sa[0]=='m') bx[n2]=atof(sa+1)/xsize; else bx[n2]=atof(sa);
		ffscanf1(); if(sa[0]=='m') by[n2]=atof(sa+1)/ysize; else by[n2]=atof(sa);
		}
	/**/
	/* Rock type read */
	ffscanf1();
	/**/
	/* Add Markers in rectangle */
	if (yn==1)
		{
		printf("Set new marker location \n");
		dx=xstpx/(double)(mnumx);
		dy=ystpy/(double)(mnumy);
		/* Marker Types add */
		for (m1=(int)(bx[0]*(double)(xnumx-1));m1<(int)(bx[3]*(double)(xnumx-1)+1.0);m1++)
		for (m2=(int)(by[0]*(double)(ynumy-1));m2<(int)(by[3]*(double)(ynumy-1)+1.0);m2++)
			{
			x=(double)(m1)*xstpx+dx/2.0;
			for (m4=0;m4<mnumx;m4++)
				{
				y=(double)(m2)*ystpy+dy/2.0;
				for (m5=0;m5<mnumy;m5++)
					{
					/* X,Y */
					markx[marknum]=(float)(x);
					marky[marknum]=(float)(y);
					markt[marknum]=n1;
					/* Random Nonstability Set on markers */
					if(nonstab>0)
						{
						if(markt[marknum]<100)
							{
							/* Random nonstability Set on X,Y */
							markx[marknum]+=(float)(rand() % (nonstab*2+1) - nonstab)/((float)(nonstab))*(float)(xnonstab);
							marky[marknum]+=(float)(rand() % (nonstab*2+1) - nonstab)/((float)(nonstab))*(float)(ynonstab);
							}
						}
					/* Add Y */
					y+=dy;
					/* Add mark Num */
					marknum++;
					}
				/* Add X */
				x+=dx;
				}
			}
		}
	/**/
	/**/
	/* Set cur Body on markers */
	mm3=marknum;
	for (mm1=0;mm1<marknum;mm1++)
	if(yn!=1)
		{
		/* Rectangle yn<3 or circular body yn>=3 */
		if(yn<3)
			{
			/* RECTANGLES -------------------------- */
			/* Rectangles */
			/* Edges */
			/* 0  2 */
			/* 1  3 */
			/**/
			/* Left bond check */
			y=(double)(marky[mm1])/ysize-by[0];
			x=xsize*(bx[0]+y*(bx[1]-bx[0])/(by[1]-by[0]));
			if((double)(markx[mm1])>=x)
				{
				/* Right bond check */
				y=(double)(marky[mm1])/ysize-by[2];
				x=xsize*(bx[2]+y*(bx[3]-bx[2])/(by[3]-by[2]));
				if((double)(markx[mm1])<=x)
					{
					/* Up bond check */
					x=(double)(markx[mm1])/xsize-bx[0];
					y=ysize*(by[0]+x*(by[2]-by[0])/(bx[2]-bx[0]));
					if((double)(marky[mm1])>=y)
						{
						/* Down bond check */
						x=(double)(markx[mm1])/xsize-bx[1];
						y=ysize*(by[1]+x*(by[3]-by[1])/(bx[3]-bx[1]));
						if((double)(marky[mm1])<=y)
							{
							/* Rock type set */
							if(yn==0 && n5==0) {markt[mm1]=n1;}
							if(yn==2 && n5==0) {markt[mm1]+=n1;}
							if(yn==0 && n5==1) {markt[mm1]=n1+100;}
							if(yn==2 && n5==1) {markt[mm1]+=n1; if(markt[mm1]<100) markt[mm1]+=100;}
							if(yn==0 && n5==2) {markt[mm3]=n1; markx[mm3]=markx[mm1]; marky[mm3]=marky[mm1]; mm3++;}
							if(yn==2 && n5==2) {markt[mm3]=markt[mm1]+n1; markx[mm3]=markx[mm1]; marky[mm3]=marky[mm1]; mm3++;}
							if(yn==0 && n5==3) {markt[mm3]=n1+100; markx[mm3]=markx[mm1]; marky[mm3]=marky[mm1]; mm3++;}
							if(yn==2 && n5==3) {markt[mm3]=markt[mm1]+n1; markx[mm3]=markx[mm1]; marky[mm3]=marky[mm1]; if(markt[mm3]<100) markt[mm3]+=100; mm3++;}
							}
						}
					}
				}
			}
		else
			{
			/* CIRCULAR BODIES ----------------------- */
			/* 0 X1,Y1  2 R1,R2 */
			/* 1 X2,Y2  3 L1,L2 */
			/* Calculate and check distance from first center */
			x=(double)(markx[mm1])/xsize;
			y=(double)(marky[mm1])/ysize;
			x1=x-bx[0];
			y1=(y-by[0])*ysize/xsize;
			if((x1*x1+y1*y1)<=(bx[2]*bx[2]))
				{
				/* Calculate and check distance from second center */
				x2=x-bx[1];
				y2=(y-by[1])*ysize/xsize;
				if((x2*x2+y2*y2)>=(by[2]*by[2]*ysize*ysize/xsize/xsize))
					{
					/* Calculate and check angle relatively first or second center */
					x=x1;y=y1; if(by[3]<0) {x=x2; y=y2;} 
					ival=pow(x*x+y*y,0.5); if(!ival) ival=1.0;
					ival=x/ival; ival=asin(ABSV(ival))/pival*180.0;
					if(x>=0 && y>0) ival=180.0-ival;
					if(x<0 && y>0) ival+=180.0;
					if(x<0 && y<=0) ival=360.0-ival;
					if(ival>=ABSV(bx[3]) && ival<=ABSV(by[3]))
						{
						/* Rock type set */
						if(yn==3 && n5==0) {markt[mm1]=n1;}
						if(yn==4 && n5==0) {markt[mm1]+=n1;}
						if(yn==3 && n5==1) {markt[mm1]=n1+100;}
						if(yn==4 && n5==1) {markt[mm1]+=n1; if(markt[mm1]<100) markt[mm1]+=100;}
						if(yn==3 && n5==2) {markt[mm3]=n1; markx[mm3]=markx[mm1]; marky[mm3]=marky[mm1]; mm3++;}
						if(yn==4 && n5==2) {markt[mm3]=markt[mm1]+n1; markx[mm3]=markx[mm1]; marky[mm3]=marky[mm1]; mm3++;}
						if(yn==3 && n5==3) {markt[mm3]=n1+100; markx[mm3]=markx[mm1]; marky[mm3]=marky[mm1]; mm3++;}
						if(yn==4 && n5==3) {markt[mm3]=markt[mm1]+n1; markx[mm3]=markx[mm1]; marky[mm3]=marky[mm1]; if(markt[mm1]<100) markt[mm1]+=100; mm3++;}
						}
					}
				}
			}
		}
	/* Add markers num */
	marknum=mm3;
	}
printf("Total amount of markers = %ld\n",marknum);
/**/
/**/
/**/
/* Sinusoidal Nonstability Set on markers */
if(nonstab<0)
	{
	printf("Set Sin nonstability on markers\n");
	for (mm1=0;mm1<marknum;mm1++)
	if(markt[mm1]<100)
		{
		if(nonstab==-1)
			{
			marky[mm1]+=(float)(wamp*sin(2.0*pival/wlen*((double)(markx[mm1])-xsize*0.5+wlen*0.25)));
			}
		else
			{
			ival=0;
			if((double)(marky[mm1])>wymin && (double)(marky[mm1])<=wymid) ival=wamp*((double)(marky[mm1])-wymin)/(wymid-wymin);
			if((double)(marky[mm1])>wymid && (double)(marky[mm1])<wymax) ival=wamp*((double)(marky[mm1])-wymid)/(wymax-wymid);
			marky[mm1]+=(float)(ival*sin(2.0*pival/wlen*((double)(markx[mm1])-xsize*0.5+wlen*0.25)));
			}
		}
	}
/**/
/**/
/**/
/* Random distribution within cells */
if(1==0 && nonstab>0 && (mnumx>1 || mnumy>1))
	{
	/* Marker cycle */
	m3=mnumx*mnumy;
	for (mm1=0;mm1<marknum;mm1+=m3)
	for (mm2=mm1;mm2<mm1+m3;mm2++)
		{
		/* Random swap of markers within cells */
		m1=(long int) ((double)(rand()%(m3)));
		m4=mm1+m1;
/*
		m1=rand() % (mnumx);
		m2=rand() % (mnumy);
		m4=mm2;
printf("%ld %ld %ld %ld %ld %ld %ld",mnumx,mnumy,m1,m2,mm1,mm2,m4);getchar();
*/
		if(mm2<marknum && m4<marknum)
			{
			ival=(double)markx[mm2]; markx[mm2]=markx[m4]; markx[m4]=(float)ival;
			ival=(double)marky[mm2]; marky[mm2]=marky[m4]; marky[m4]=(float)ival;
			ival=(double)markk[mm2]; markk[mm2]=markk[m4]; markk[m4]=(float)ival;
			n1=markt[mm2]; markt[mm2]=markt[m4]; markt[m4]=n1;
/*
*/
			}
		}
	}
/**/
/**/
/**/
/* Initial distribution of Temperature Read, Set to Grid */
printf("Read distribution of temperature \n");
ffscanf1();
while (sa[0]!='~')
	{
	/* Read Type of current T distribution: */
	/* 0 - simple rectangle */
	/* 1 - simple rectangle with horizontal exponential T distribution */
	/* 2 - simple rectangle with vertical exponential T distribution */
	/* 3 - simple circular distribution */
	/* 4 - Rectangle with cooling age */
	/* 5 - Transition 4->0 */
	/* 6 - Transition 0->4 */
	/* 8 - Impact structure following parametrization of Monteux et al. (2007) */
	yn=atoi(sa);
	/* Read coordinats of Cur boxe */
	/* 0  2 */
	/* 1  3 */
	for(n2=0;n2<4;n2++)
		{
		ffscanf1(); if(sa[0]=='m') bx[n2]=atof(sa+1)/xsize; else bx[n2]=atof(sa);
		ffscanf1(); if(sa[0]=='m') by[n2]=atof(sa+1)/ysize; else by[n2]=atof(sa);
		}
	/**/
	/* Read T K in Edges of Cur boxe */
	ffscanf1(); t[0]=atof(sa);
	ffscanf1(); t[1]=atof(sa);
	ffscanf1(); t[2]=atof(sa);
	ffscanf1(); t[3]=atof(sa);
	if(yn==4 || yn==5 || yn==6) 
		{
		ffscanf1(); t[4]=atof(sa);
		ffscanf1(); t[5]=atof(sa);
		ffscanf1(); t[6]=atof(sa);
		}
	ffscanf1();
	/**/
	/* Initial distribution of Temperature set for cur box */
	for (m1=0;m1<xnumx;m1++)
	for (m2=0;m2<ynumy;m2++)
		{
		/* Node Num */
		m3=m1*ynumy+m2;
		/* Node X,Y General coordinats calc */
		x=gx[m1]/xsize;
		y=gy[m2]/ysize;
		if(yn!=3 && yn!=7)
			{
			/* Node X General coordinat check */
			x1=bx[0]+(y-by[0])*(bx[1]-bx[0])/(by[1]-by[0]);
			x2=bx[2]+(y-by[2])*(bx[3]-bx[2])/(by[3]-by[2]);
			if(x>=x1 && x<=x2)
				{
				/* Node Y limit box coord calc */
				y1=by[0]+(x-bx[0])*(by[2]-by[0])/(bx[2]-bx[0]);
				y2=by[1]+(x-bx[1])*(by[3]-by[1])/(bx[3]-bx[1]);
				/* Node Y general coordinat check */
				if(y>=y1 && y<=y2)
					{
					/* Node XY box coord calc */
					x=(x-x1)/(x2-x1);
					y=(y-y1)/(y2-y1);
/*
printf("%e %e %ld %ld",x,y,m1,m2);getchar();
*/
					/**/
					if(yn==0)
						{
						/* Simple  T recalc after relativ coord */
						tk[m3]=(t[0]+y*(t[1]-t[0]))*(1.0-x)+(t[2]+y*(t[3]-t[2]))*x;
						}
					if(yn==1)
						{
						/* Spec Horis T recalc after relativ coord */
						ival=y*(t[1]-t[0]);
						akf=exp(-markdh[0]/8.314/t[0]);
						bkf=exp(-markdh[0]/8.314/t[2])-akf;
						tk[m3]=ival-markdh[0]/8.314/log(akf+bkf*x);
						}
					if(yn==2)
						{
						/* Spec Vert  T recalc after relativ coord */
						ival=x*(t[2]-t[0]);
						akf=exp(-markdh[0]/8.314/t[0]);
						bkf=exp(-markdh[0]/8.314/t[1])-akf;
						tk[m3]=ival-markdh[0]/8.314/log(akf+bkf*y);
						}
					if(yn==4)
						{
						/* Turcotte & Schubert, 2002: Cooling age T=T1-erfc()(T1-T0) */
						ival=(gy[m2]-y1*ysize)/2.0/pow(t[6]*(t[4]+(t[5]-t[4])*x)*3.15576e+7,0.5);
						tk[m3]=(t[1]+(t[3]-t[1])*x)-(1.0-erf(ival))*((t[1]-t[0])*(1.0-x)+(t[3]-t[2])*x);
/*
printf("%ld %ld %ld  %e %e   %e %e   %e",m1,m2,m3,gx[m1],gy[m2],x,y,tk[m3]);getchar();
*/
						}
					if(yn==5)
						{
						/* Transition 4->0 */
						ival=(gy[m2]-y1*ysize)/2.0/pow(t[6]*(t[4]+(t[5]-t[4])*x)*3.15576e+7,0.5);
						ival=(t[1]+(t[3]-t[1])*x)-(1.0-erf(ival))*((t[1]-t[0])*(1.0-x)+(t[3]-t[2])*x);
						tk[m3]=(t[0]+y*(t[1]-t[0]))*(1.0-x)+(t[2]+y*(t[3]-t[2]))*x;
						tk[m3]=ival*(1.0-x)+tk[m3]*x;
						}
					if(yn==6)
						{
						/* Transition 0->4 */
						ival=(gy[m2]-y1*ysize)/2.0/pow(t[6]*(t[4]+(t[5]-t[4])*x)*3.15576e+7,0.5);
						ival=(t[1]+(t[3]-t[1])*x)-(1.0-erf(ival))*((t[1]-t[0])*(1.0-x)+(t[3]-t[2])*x);
						tk[m3]=(t[0]+y*(t[1]-t[0]))*(1.0-x)+(t[2]+y*(t[3]-t[2]))*x;
						tk[m3]=ival*x+tk[m3]*(1.0-x);
						}
/*
printf("%e %e %ld %ld %ld %e",x,y,m1,m2,m3,tk[m3]);getchar();
*/
					}
				}
			}
		else
			{
			if (yn==3)
				{
				/* CIRCULAR BODIES ----------------------- */
				/* 0 X1,Y1  2 R1,R2 */
				/* 1 X2,Y2  3 L1,L2 */
				/* Calculate and check distance from first center */
				x1=x-bx[0];
				y1=(y-by[0])*ysize/xsize;
				x4=bx[2]-pow(x1*x1+y1*y1,0.5);
				if(x4>=0)
					{
					/* Calculate and check distance from second center */
					x2=x-bx[1];
					y2=(y-by[1])*ysize/xsize;
					y4=pow(x2*x2+y2*y2,0.5)-by[2]*ysize/xsize;
					if(y4>=0)
						{
						/* Calculate and check angle relatively first or second center */
						x3=x1;y3=y1; if(by[3]<0) {x3=x2; y3=y2;} 
						ival=pow(x3*x3+y3*y3,0.5); if(!ival) ival=1.0;
						ival=x3/ival; ival=asin(ABSV(ival))/pival*180.0;
						if(x3>=0 && y3>0) ival=180.0-ival;
						if(x3<0 && y3>0) ival+=180.0;
						if(x3<0 && y3<=0) ival=360.0-ival;
						if(ival>=ABSV(bx[3]) && ival<=ABSV(by[3]))
							{
							/* Relative coordinats calc */
							x4+=y4;if(!x4) x4=1.0;
							x=y4/x4;
							x4=ABSV(by[3])-ABSV(bx[3]); if(!x4) x4=1.0;
							y=(ival-ABSV(bx[3]))/x4;
							/* T recalc after relativ coord */
							tk[m3]=(t[0]*y+t[1]*(1.0-y))*x+(t[2]*y+t[3]*(1.0-y))*(1.0-x);
							}
						}
					}
				}
			if (yn==7)
				{
				/* IMPACT STRUCTURE using parametrization of Monteux et al. (2007)---- */
				/* 0 X1,Y1  2 R1,R2 */
				/* 1 X2,Y2  3 L1,L2 */
				/* Calculate and check distance from first center */
				x1=x-bx[0];                 /* non-dim. x distance from inner centre of isobaric core */
				y1=(y-by[0])*ysize/xsize;   /* non-dim. y distance from inner centre of isobaric core */
				x4=pow(x1*x1+y1*y1,0.5);    /* non-dim. total distance from inner centre of isobaric core*/
                                if(x4>bx[2])                /* point is within the isobaric core */
					{
					/* Calculate and check distance from second center */
					x2=x-bx[1];                 /* non-dim. x distance from outer centre of isobaric core */
					y2=(y-by[1])*ysize/xsize;   /* non-dim. y distance from outer centre of isobaric core */
					y4=pow(x2*x2+y2*y2,0.5);    /* non-dim. total distance from outer centre of isobaric core */
/* Impose decaying thermal anomaly, when outside the isobaric core */
/* Isobaric core itself should be prescribed using a circular structure (yn==3) */
					if(y4>bx[2])
						{
						tk[m3]+=t[0]*pow((bx[2]/y4),4.400);
/*
printf("%ld %ld %e %e %e %e %e %e %e",m1,m2,bx[1],by[1],bx[2],by[2],y4,t[0],t[0]*pow((bx[2]/y4),4.4));getchar();
*/
						}
					}
				}
			}
		}
	}
/**/
/* Set initial porosity and grain size of all materials */
/* added by Gregor (last updated: 24/04/2015) */
for(mm1=0;mm1<marknum;mm1++)
	{
	/* No porosity for sticky air and iron material */
	/* Initial grain size for sticky air and iron material */
	if(markt[mm1]==0 || markt[mm1]>=7)
		{
		markpor[mm1] = 0.000;
		markgr[mm1]  = (float)(gr_init/1.000e+4);
		}
	/**/
	/* Porosity and grain size for solid silicate material */
	if(markt[mm1]==5 || markt[mm1]==6)
		{
		markpor[mm1] = (float)(por_init);
		markgr[mm1]  = (float)(gr_init);
		}
	}
/**/
/**/
/* Set start accretion time of all materials */
/* Fabio/Tim (last update: 28/12/2016) */
for(mm1=0;mm1<marknum;mm1++)
        {
        /* Start time as initial accretion time, NaN for sticky air */
        if(markt[mm1]!=0)
        	{
		markacc[mm1] = (float)(timesum/(3.15576e+7));
        	}
    	}
/**/
/**/
/* Set initial magnetization and magnetization time for all markers */
/* added by Gregor (18/10/2011) */
for(mm1=0;mm1<marknum;mm1++)
	{
	markmg_old[mm1]=0;
	markmg_time[mm1]=-1.000;
	}
/**/
/**/
/* Calc temperature for markers */
for (mm1=0;mm1<marknum;mm1++)
/* Check markers out of grid */
if (markx[mm1]>0 && marky[mm1]>0 && (double)(markx[mm1])<xsize && (double)(marky[mm1])<ysize)
	{
	allintert((double)(markx[mm1]),(double)(marky[mm1]));
	markk[mm1]=(float)(eps[2]);
	if(markt[mm1]!=0)
                {
		marktmax[mm1]=(float)(markk[mm1]); /* Tim: Write initial temperature field for marker max */
		}
	}
/**/
/**/
/* 2-nd ro[],nu[], cp[] kt[] for nodes Recalc after markers */
printf("2-nd Recalc nodes properties after markers types \n");
ronurecalc();
/**/
/**/
/**/
/* Close Configuration file */
fclose(fl1);
/**/
/* Third  Calc,Check Grid parameters */
gridcheck();
/**/
/* Call Impact subroutine */
impact();
/**/
/* Check core evolution (Greg addition) */
core();
/**/
/* Save information in output file */
printf("Saving information to output file \n");
printmod=1;
saver(0,0);
/**/
/*  Save impact and core history (Gregor addition) */
impactsave();
/**/
/**/
printf("OK!\n");
return 0;
}
/* Formation of Data file for i2mart.c program */
