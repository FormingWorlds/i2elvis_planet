/* LEFT+Right Side or Err for X-Stokes Equation */
/* Stokes equation initial form */
/* dSIGxx/dX + dSIGxy/dY - dP/dX = - RO*Gx */
double xstokserr(long int m1, long int m2, int ynerr)
/* m1,m2 - node X,Y number */
/* ynerr - Err Calc Y(1)/N(0) */
{
/* Counters */
int n1,n2;
long int v[4];
/* Err Buf */
double leftx,rightx,nueff;
/* Distances */
double xkf=(gx[m1+1]-gx[m1-1])/2.0,ykf=gy[m2+1]-gy[m2];
double GXKOEFX=GXKOEF;
/**/
/**/
/**/
/* Staggered Nodes num */
/*                      [0]                [2]   */
/*                  Ro0,Sxy0,Nu0                 */
/*                                               */
/*         Pr1,Sxx1    <Vx0>  Pr3,Sxx3           */
/*                                               */
/*                      [1]                [3]   */
/*                  RO1,Sxy1,Nu1                 */
/*                                               */
v[0]=m1*ynumy+m2;v[1]=v[0]+1;
v[2]=v[0]+ynumy;v[3]=v[2]+1;
/**/
/**/
/**/
/* Spherical Gravity acceleration from gravity potential */
if(0==0)
	{
	GXKOEFX=-(gp[v[3]]-gp[v[1]])/xkf;
/*
printf("GX %ld %ld  %e %e  %e    %e %e ",m1,m2,gp[v[3]],gp[v[1]],GXKOEFX,gp[v[3]],GXKOEF);getchar();
*/
	}
/**/
/**/
/**/
/* RIGHT parts of X-Stokes */
/* dSIGik/dXk - dP/dXi = - Gi*RO */
rightx  = -GXKOEFX*(ro[v[0]]+ro[v[1]])/2.0;
/**/
/**/
/**/
/* Return val for LSQ err ----------------------------*/
if(ynerr==1)
	{
	/* LEFT part of X-Stokes */
	/* dSIGxx/dX + dSIGxy/dY - dP/dX = - RO*Gx */
	/**/
	/* dSIGxx/dX */
	leftx =(sxx[v[3]]-sxx[v[1]])/xkf;
	/* dSIGxy/dY */
	leftx+=(sxy[v[1]]-sxy[v[0]])/ykf;
	/* -dP/dX */
	leftx-=(pr[v[3]]-pr[v[1]])/xkf;
	/**/
	/* FD stabilization mechanism included [T. Duretz et al., Geochem. Geophys. Geosyst., 12, Q07004, 2011] */
	/* -gx/dt/2*(vx*dRHO/dx+vy*dRHO/dy) */
	leftx-=timestep/2.0*GXKOEFX*(vx[v[0]]*(ro[v[2]]+ro[v[3]]-ro[v[0]-ynumy]-ro[v[1]-ynumy])/(gx[m1+1]-gx[m1-1])/2.0+((vy[v[0]]+vy[v[1]])*(gx[m1]-gx[m1-1])+(vy[v[0]-ynumy]+vy[v[1]-ynumy])*(gx[m1+1]-gx[m1]))/(gx[m1+1]-gx[m1-1])/2.0*(ro[v[1]]-ro[v[0]])/(gy[m2+1]-gy[m2]));
	/**/
	/**/
	/* Effective NU calc */
	nueff=(nu[v[0]]+nu[v[1]])/2.0;
	/**/
	/* X-STOKES Error */
	leftx=(leftx-rightx)/nueff;
	/**/
	/* Min,Max Value of Nu Save */
	errbuf[9]=nueff;
	/* Min,Max Value of P Save */
	errbuf[10]=MAXV(pr[v[1]],pr[v[3]]);
	errbuf[11]=MINV(pr[v[1]],pr[v[3]]);
	/**/
	return leftx;
	}
/**/
/**/
/**/
/* Set Initial Num of lines -------------------------------------- */
wn[0]=2;
/* Save Right part Save for X-Stokes ---------------------*/
wi[0]=rightx;
/**/
/* Add Coefficients for left parts of X-Stokes ----------------*/
/* Staggered Nodes num */
/*                      [0]                [2]   */
/*                  Ro0,Sxy0,Nu0                 */
/*                                               */
/*         Pr1,Sxx1    <Vx0>  Pr3,Sxx3           */
/*                                               */
/*                      [1]                [3]   */
/*                  RO1,Sxy1,Nu1                 */
/*                                               */
/*  0(P) 1(Vx)  2(Vy)  */
/* -dP/dX */
wn[1]=v[1]*3+0;
wi[1]=+1.0/xkf;
wn[2]=v[3]*3+0;
wi[2]=-1.0/xkf;
/* dSIGxx/dX */
sxxcalc(m1,m2+1,-1.0/xkf);
sxxcalc(m1+1,m2+1,1.0/xkf);
/* dSIGxy/dY */
sxycalc(m1,m2,-1.0/ykf);
sxycalc(m1,m2+1,1.0/ykf);
/**/
/**/
/* FD stabilization mechanism included [T. Duretz et al., Geochem. Geophys. Geosyst., 12, Q07004, 2011] */
/* -gx/dt/2*(vx*dRHO/dx+vy*dRHO/dy) */
wn[wn[0]+1]=v[0]*3+1;
wi[wn[0]+1]=-timestep/2.0*GXKOEFX*(ro[v[2]]+ro[v[3]]-ro[v[0]-ynumy]-ro[v[1]-ynumy])/(gx[m1+1]-gx[m1-1])/2.0;
wn[wn[0]+2]=v[0]*3+2;
wi[wn[0]+2]=-timestep/2.0*GXKOEFX*(gx[m1]-gx[m1-1])/(gx[m1+1]-gx[m1-1])/2.0*(ro[v[1]]-ro[v[0]])/(gy[m2+1]-gy[m2]);
wn[wn[0]+3]=v[1]*3+2;
wi[wn[0]+3]=-timestep/2.0*GXKOEFX*(gx[m1]-gx[m1-1])/(gx[m1+1]-gx[m1-1])/2.0*(ro[v[1]]-ro[v[0]])/(gy[m2+1]-gy[m2]);
wn[wn[0]+4]=(v[0]-ynumy)*3+2;
wi[wn[0]+4]=-timestep/2.0*GXKOEFX*(gx[m1+1]-gx[m1])/(gx[m1+1]-gx[m1-1])/2.0*(ro[v[1]]-ro[v[0]])/(gy[m2+1]-gy[m2]);
wn[wn[0]+5]=(v[1]-ynumy)*3+2;
wi[wn[0]+5]=-timestep/2.0*GXKOEFX*(gx[m1+1]-gx[m1])/(gx[m1+1]-gx[m1-1])/2.0*(ro[v[1]]-ro[v[0]])/(gy[m2+1]-gy[m2]);
wn[0]+=5;
/**/
/**/
/**/
/*
for(n1=0;n1<=vn[0];n1++)printf("Vx %e %ld \n",wi[n1],wn[n1]);getchar();
*/
return 0;
}
/* Left+Right Side or Err for X-Stokes Equation */




/* LEFT+Right Side or Err for Y-Stokes Equation */
/* Stokes equation initial form */
/* dSIGyy/dY + dSIGxy/dX - dP/dY = - RO*Gy */
double ystokserr(long int m1, long int m2, int ynerr)
/* m1,m2 - node X,Y number */
/* ynerr - Err Calc Y(1)/N(0) */
{
/* Counters */
int n1;
long int v[4];
/* Err Buf */
double lefty,righty,nueff;
/* Distances */
double xkf=gx[m1+1]-gx[m1],ykf=(gy[m2+1]-gy[m2-1])/2.0;
double GYKOEFY=GYKOEF;
/**/
/**/
/**/
/* Staggered Nodes num */
/*                                               */
/*                Pr2,Syy2                       */
/*                                               */
/*         [0]                  [2]              */
/*    Ro0,Sxy0,Nu0   <Vy0>   Ro2,Sxy2,Nu2        */
/*                                               */
/*                Pr3,Syy3                       */
/*                                               */
/*         [1]                  [3]              */
/*                                               */
v[0]=m1*ynumy+m2;v[1]=v[0]+1;
v[2]=v[0]+ynumy;v[3]=v[2]+1;
/**/
/**/
/**/
/* Spherical Gravity */
if(0==0)
	{
	GYKOEFY=-(gp[v[3]]-gp[v[2]])/ykf;
/*
printf("GY %ld %ld  %e %e  %e    %e %e ",m1,m2,gp[v[3]],gp[v[1]],GYKOEFY,gp[v[3]],GXKOEF);getchar();
*/
	}
/**/
/**/
/**/
/* RIGHT part of Y-Stokes */
/* dSIGik/dXk - dP/dXi = - Gi*RO */
righty  = -GYKOEFY*(ro[v[0]]+ro[v[2]])/2.0;
/**/
/**/
/**/
/**/
/* Return val for LSQ err ----------------------------*/
if(ynerr==1)
	{
	/* LEFT part of Y-Stokes */
	/* dSIGyy/dY + dSIGxy/dX - dP/dY = - RO*Gy */
	/**/
	/* dSIGyy/dY */
	lefty =(syy[v[3]]-syy[v[2]])/ykf;
	/* dSIGxy/dX */
	lefty+=(sxy[v[2]]-sxy[v[0]])/xkf;
	/* -dP/dY */
	lefty-=(pr[v[3]]-pr[v[2]])/ykf;
	/**/
	/**/
	/* FD stabilization mechanism included [T. Duretz et al., Geochem. Geophys. Geosyst., 12, Q07004, 2011] */
	/* -gy/dt/2*(vy*dRHO/dy+vx*dRHO/dx) */
	lefty-=timestep/2.0*GYKOEFY*(vy[v[0]]*(ro[v[1]]+ro[v[3]]-ro[v[0]-1]-ro[v[2]-1])/(gy[m2+1]-gy[m2-1])/2.0+((vx[v[0]]+vx[v[2]])*(gy[m2]-gy[m2-1])+(vx[v[0]-1]+vx[v[2]-1])*(gy[m2+1]-gy[m2]))/(gy[m2+1]-gy[m2-1])/2.0*(ro[v[2]]-ro[v[0]])/(gx[m1+1]-gx[m1]));
	/**/
	/**/
	/* Effective NU calc */
	nueff=(nu[v[0]]+nu[v[2]])/2.0;
	/**/
	/* Y-STOKES Error */
	lefty=(lefty-righty)/nueff;
	/**/
	/* Min,Max Value of Nu Save */
	errbuf[9]=nueff;
	/* Min,Max Value of P Save */
	errbuf[10]=MAXV(pr[v[2]],pr[v[3]]);
	errbuf[11]=MINV(pr[v[2]],pr[v[3]]);
	/**/
	return lefty;
	}
/**/
/**/
/**/
/* Set Initial Num of lines -------------------------------------- */
wn[0]=2;
/* Save Right parts Save for Y-Stokes ---------------------*/
wi[0]=righty;
/**/
/* Add Coefficients for left parts of Y-Stokes ----------------*/
/* Staggered Nodes num */
/*                                               */
/*                Pr2,Syy2                       */
/*                                               */
/*         [0]                  [2]              */
/*    Ro0,Sxy0,Nu0   <Vy0>   Ro2,Sxy2,Nu2        */
/*                                               */
/*                Pr3,Syy3                       */
/*                                               */
/*         [1]                  [3]              */
/*                                               */
/*  0(P) 1(Vx)  2(Vy)  */
/* -dP/dY */
wn[1]=v[2]*3+0;
wi[1]=+1.0/ykf;
wn[2]=v[3]*3+0;
wi[2]=-1.0/ykf;
/* dSIGyy/dY */
syycalc(m1+1,m2,-1.0/ykf);
syycalc(m1+1,m2+1,1.0/ykf);
/* dSIGxy/dX */
sxycalc(m1,m2,-1.0/xkf);
sxycalc(m1+1,m2,1.0/xkf);
/**/
/**/
/* FD stabilization mechanism included [T. Duretz et al., Geochem. Geophys. Geosyst., 12, Q07004, 2011] */
/* -gy/dt/2*(vy*dRHO/dy+vx*dRHO/dx) */
wn[wn[0]+1]=v[0]*3+2;
wi[wn[0]+1]=-timestep/2.0*GYKOEFY*(ro[v[1]]+ro[v[3]]-ro[v[0]-1]-ro[v[2]-1])/(gy[m2+1]-gy[m2-1])/2.0;
wn[wn[0]+2]=v[0]*3+1;
wi[wn[0]+2]=-timestep/2.0*GYKOEFY*(gy[m2]-gy[m2-1])/(gy[m2+1]-gy[m2-1])/2.0*(ro[v[2]]-ro[v[0]])/(gx[m1+1]-gx[m1]);
wn[wn[0]+3]=v[2]*3+1;
wi[wn[0]+3]=-timestep/2.0*GYKOEFY*(gy[m2]-gy[m2-1])/(gy[m2+1]-gy[m2-1])/2.0*(ro[v[2]]-ro[v[0]])/(gx[m1+1]-gx[m1]);
wn[wn[0]+4]=(v[0]-1)*3+1;
wi[wn[0]+4]=-timestep/2.0*GYKOEFY*(gy[m2+1]-gy[m2])/(gy[m2+1]-gy[m2-1])/2.0*(ro[v[2]]-ro[v[0]])/(gx[m1+1]-gx[m1]);
wn[wn[0]+5]=(v[2]-1)*3+1;
wi[wn[0]+5]=-timestep/2.0*GYKOEFY*(gy[m2+1]-gy[m2])/(gy[m2+1]-gy[m2-1])/2.0*(ro[v[2]]-ro[v[0]])/(gx[m1+1]-gx[m1]);
wn[0]+=5;
/**/
/**/
/**/
/*
for(n1=0;n1<=vn[0];n1++)printf("Vy %e %ld \n",wi[n1],wn[n1]);getchar();
*/
return 0;
}
/* Left+Right Side or Err for Y-Stokes Equation */




/* Left side or Err for Continuity Equation  */
/* dVx/dX+dVy/dY = 0 */
double conterr(long int m1, long int m2, int ynerr)
/* m1,m2 - node X,Y number */
/* ynerr - Err Calc Y(1)/N(0) */
{
/* Counter */
long int v[4],m1min,m1max,m2min,m2max,m3;
/* Val Buffer */
double divv=0,xi;
int n1,nx,ny;
/**/
/* Staggered Nodes num */
/*   [0]       Vy0      [2] */
/*                          */
/*   Vx0        <P3>    Vx2 */
/*            Exx3,Eyy3     */
/*                          */
/*   [1]       Vy1      [3] */
v[0]=(m1-1)*ynumy+(m2-1);v[1]=v[0]+1;
v[2]=v[0]+ynumy;v[3]=v[2]+1;
/**/
/**/
/**/
/* Calc, Check Fd limits for dVx/dX */
m1min=m1-1-stoksfd; if(m1min<0) m1min=0;
m1max=m1+stoksfd; if(m1max>xnumx-1) m1max=xnumx-1;
/**/
/* Current X position */
xi=(gx[m1-1]+gx[m1])/2.0;
/**/
/* Load distances to xn[] */
for (m3=m1min;m3<=m1max;m3++)
	{
	xn[m3-m1min]=gx[m3];
	}
/**/
/* Calc maximal position in xn[] */
nx=(int)(m1max-m1min);
/**/
/* Calc Vx coefficients for EPSxx */
fdweight(nx,1,xi);
/* Reload coefficients to cn[] */
for (m3=0;m3<=nx;m3++)
	{
	cn[m3][2]=cn[m3][1];
	}
/**/
/**/
/**/
/* Calc, Check Fd limits dVy/dY */
m2min=m2-1-stoksfd; if(m2min<0) m2min=0;
m2max=m2+stoksfd; if(m2max>ynumy-1) m2max=ynumy-1;
/**/
/* Current Y position */
xi=(gy[m2-1]+gy[m2])/2.0;
/**/
/* Load distances to xn[] */
for (m3=m2min;m3<=m2max;m3++)
	{
	xn[m3-m2min]=gy[m3];
	}
/**/
/* Calc maximal position in xn[] */
ny=(int)(m2max-m2min);
/**/
/* Calc Vy coefficients for EPSyy */
fdweight(ny,1,xi);
/**/
/* Return dVx/dX+dVy/dY err ----------------------------*/
if(ynerr==1)
	{
	/* dVx/dX */
	/* Add Vx with koefficients */
	for (m3=m1min;m3<=m1max;m3++)
		{
		divv+=vx[m3*ynumy+(m2-1)]*cn[m3-m1min][2];
		}
	/**/
	/* dVy/dY */
	/* Add Vy with koefficients */
	for (m3=m2min;m3<=m2max;m3++)
		{
		divv+=vy[(m1-1)*ynumy+m3]*cn[m3-m2min][1];
		}
	/**/
	return divv;
	}
/**/
/**/
/**/
/* Add continuity equation */
/* Set Initial Num of lines -------------------------------------- */
wn[0]=0;
/**/
/* Save Right part for Contin ---------------------*/
wi[0]=0;
/**/
/**/
/**/
/* Add Coefficients for left parts of dVx/dX ----------------*/
/*  0(P) 1(Vx)  2(Vy)  */
/* Sxx=2Nu*Exx, Exx=dVx/dX */
/* Add Vx with koefficients */
for (m3=m1min;m3<=m1max;m3++)
	{
	wn[wn[0]+1+m3-m1min]=(m3*ynumy+(m2-1))*3+1;
	wi[wn[0]+1+m3-m1min]=cn[m3-m1min][2];
	}
/**/
/* Add total Num of lines */
wn[0]+=nx+1;
/**/
/**/
/**/
/* Add Coefficients for left parts of Syy ----------------*/
/*  0(P) 1(Vx)  2(Vy)  */
/* Sxx=2Nu*Exx, Exx=dVx/dX */
/* Add Vy with koefficients */
for (m3=m2min;m3<=m2max;m3++)
	{
	wn[wn[0]+1+m3-m2min]=((m1-1)*ynumy+m3)*3+2;
	wi[wn[0]+1+m3-m2min]=cn[m3-m2min][1];
	}
/**/
/* Add total Num of lines */
wn[0]+=ny+1;
/**/
/**/
/**/
/* Check Boundary conditions around Cell */
divv=1.0;
if (!bondm[v[0]*3+1]) divv=0;
if (!bondm[v[0]*3+2]) divv=0;
if (!bondm[v[1]*3+2]) divv=0;
if (!bondm[v[2]*3+1]) divv=0;
/**/ 
/*
for(n1=0;n1<3;n1++)printf("Cont %e %d \n",wi[n1],wn[n1]);getchar();
*/
return divv;
}
/* Left side or Err for Continuity Equation */



/* Left side or Err for vX Boundary Condition Equation */ 
/* V=CONST+KOEF*Vn */
double xbonderr(long int mcmax, int ynerr)
/* mcmax - numer of cur Vx in sol[] */
/* ynerr - Err Calc Y(1)/N(0) */
{
/* Val Buffer */
double leftx=0;
int n1;
/**/
/**/
/**/
/* Error Calc */
if (ynerr)
	{
	/* Add Const */
	leftx=x[mcmax]-bondv[bondm[mcmax]][0];
	/* Add Koef */
	for (n1=0;n1<3;n1++)
		{
		if(bondn[bondm[mcmax]][n1]) 
			{
			leftx-=bondv[bondm[mcmax]][n1+1]*x[bondn[bondm[mcmax]][n1]-1];
			}
		}
	/**/
	return leftx;
	}
/**/
/**/
/**/
/* Add X CONST */
wn[0]=1;
wi[0]=bondv[bondm[mcmax]][0];
wn[1]=mcmax;
wi[1]=1.0;
/* Add X PAR1,PAR2,PAR3 */
for (n1=0;n1<3;n1++)
	{
	if(bondn[bondm[mcmax]][n1]) 
		{
		wn[0]+=1;
		wn[wn[0]]=bondn[bondm[mcmax]][n1]-1;
		wi[wn[0]]=-bondv[bondm[mcmax]][n1+1];
		}
	}
/**/
/*
for(n1=0;n1<3;n1++)printf("%e %ld \n",wi[n1],wn[n1]);getchar();
*/
return 0;
}
/* Left side or Err for vX Boundary Condition Equation */ 





/* Left side or Err for vY Boundary Condition Equation */ 
/* V=CONST+KOEF*Vn */
double ybonderr(long int mcmax, int ynerr)
/* mcmax - numer of cur Vx in sol[] */
/* ynerr - Err Calc Y(1)/N(0) */
{
/* Val Buffer */
double lefty=0;
int n1;
/**/
/**/
/**/
/* Error Calc */
if (ynerr)
	{
	/* Add Const */
	lefty=x[mcmax]-bondv[bondm[mcmax]][0];
	/* Add Koef */
	for (n1=0;n1<3;n1++)
		{
		if(bondn[bondm[mcmax]][n1]) 
			{
			lefty-=bondv[bondm[mcmax]][n1+1]*x[bondn[bondm[mcmax]][n1]-1];
			}
		}
	/**/
	return lefty;
	}
/**/
/**/
/**/
/* Add Y CONST */
wn[0]=1;
wi[0]=bondv[bondm[mcmax]][0];
wn[1]=mcmax;
wi[1]=1.0;
/* Add Y PAR1,PAR2,PAR3 */
for (n1=0;n1<3;n1++)
	{
	if(bondn[bondm[mcmax]][n1]) 
		{
		wn[0]+=1;
		wn[wn[0]]=bondn[bondm[mcmax]][n1]-1;
		wi[wn[0]]=-bondv[bondm[mcmax]][n1+1];
		}
	}
/*
for(n1=0;n1<3;n1++)printf("%e %d \n",wi[n1],wn[n1]);getchar();
*/
return 0;
}
/* Left side or Err for vY Boundary Condition Equation */ 




/* Left side or Err for P Boundary Condition Equation */ 
/* P=CONST+KOEF*Pn */
double pbonderr(long int mcmax, int ynerr)
/* mcmax - numer of cur P  in sol[] */
/* ynerr - Err Calc Y(1)/N(0) */
{
/* Val Buffer */
double leftp;
int n1;
/**/
/**/
/**/
/* Error Calc */
if (ynerr)
	{
	/* Add Const */
	leftp=x[mcmax]-bondv[bondm[mcmax]][0];
	/* Add Koef */
	for (n1=0;n1<3;n1++)
		{
		if(bondn[bondm[mcmax]][n1]) 
			{
			leftp-=bondv[bondm[mcmax]][n1+1]*x[bondn[bondm[mcmax]][n1]-1];
			}
		}
	/**/
	return leftp*leftp;
	}
/**/
/**/
/**/
/* Add P CONST */
wn[0]=1;
wi[0]=bondv[bondm[mcmax]][0];
wn[1]=mcmax;
wi[1]=1.0;
/* Add P PAR1,PAR2,PAR3 */
for (n1=0;n1<3;n1++)
	{
	if(bondn[bondm[mcmax]][n1]) 
		{
		wn[0]+=1;
		wn[wn[0]]=bondn[bondm[mcmax]][n1]-1;
		wi[wn[0]]=-bondv[bondm[mcmax]][n1+1];
		}
	}
/*
for(n1=0;n1<3;n1++)printf("%e %d \n",wi[n1],wn[n1]);getchar();
*/
return 0;
}
/* Left side or Err for P Boundary Condition Equation */ 





/* Solve XY-Stokes+Continuity equations by vX,vY,P mixed arbitrary order Finite Diff method */
void viterate(int m0)
{
/* Counters */
long int m1,m2,m3,mcmax,mcmax0,mcmax1,dm1,dm2,ccc=0,n1;
double vxmin,vxmax,vymin,vymax,minvx=1e+30,maxvx=-1e+30,minvy=1e+30,maxvy=-1e+30,mindx,pmpa,minnu=1e+50,maxnu=-1e+50,minpr=1e+50,maxpr=-1e+50;
int printyn=printmod;
/**/
/* Val buffer */
double ival,ival1,mukoef;
/* Err koef */
double bondsum=0,bondnum=0;
double stoksum=0,stoknum=0;
double contsum=0,contnum=0;
double gamma=6.67384e-11,dx02,dx24,dy12,dy23;
double pival,rosum,ronum,dx,dy;
pival=2.0*asin(1.0);
/**/
/**/
/**/
/* Gravity field compute */
if(0==0)
{
/* Compute average density of the planet */
rosum=ronum=0;
/* Node  Cycle */
for (m1=1;m1<xnumx;m1++)
for (m2=1;m2<ynumy;m2++)
	{
	/* Pos GP in sol0[] */
	mcmax=m1*ynumy+m2;
	/**/
	/* Calc Distace */
	dx=(gx[m1-1]+gx[m1]-xsize)/2.0;
	dy=(gy[m2-1]+gy[m2]-ysize)/2.0;
	ival=pow(dx*dx+dy*dy,0.5);
	/* Compose boundary condition equations */
	if(ival<GYKOEF)
		{
		ival1=(gx[m1]-gx[m1-1])*(gy[m2]-gy[m2-1]);
		rosum+=ival1*(ro[mcmax-1-ynumy]+ro[mcmax-ynumy]+ro[mcmax-1]+ro[mcmax])/4.0;
		ronum+=ival1;
		}
	}
/* Compute mass and gravity of the planet */
if(ronum>0)
	{
	rosum/=ronum;
	GXKOEF=4.0/3.0*pival*GYKOEF*rosum*gamma;
/*
printf("%e %e",rosum,GXKOEF);getchar();
printf("%e %e %e",rosum,GYKOEF,GXKOEF);getchar();
printf("%e %e %e",rosum,rosum*4.0/3.0*pival*GYKOEF*GYKOEF*GYKOEF,GXKOEF);getchar();
*/
	}
/**/
/* Clear New Solution */
/* P, Vx,Vy */
pos0cur=0;
/*
mukoef=1e-1;
*/
mukoef=2.0/xstpx/xstpx+2/ystpy/ystpy;
/**/
/**/
/* Add Matrix by gravity potential equation */
/* Node  Cycle */
for (m1=0;m1<xnumx;m1++)
for (m2=0;m2<ynumy;m2++)
	{
	/* Pos GP in sol0[] */
	mcmax=m1*ynumy+m2;
	/**/
	if(m1>0 && m2>0)
		{
		/* Calc Distace */
		dx=(gx[m1-1]+gx[m1]-xsize)/2.0;
		dy=(gy[m2-1]+gy[m2]-ysize)/2.0;
		ival=pow(dx*dx+dy*dy,0.5);
		/* Compose boundary condition equations */
/*
		if(m1==1 || m1==xnumx-1 || m2==1 || m2==ynumy-1 || ival>GYKOEF)
*/
		if(m1==1 || m1==xnumx-1 || m2==1 || m2==ynumy-1)
			{
			wn[0]=1;
			wn[1]=mcmax;
			/* Gravity potential for constant concentric acceleration */
			wi[0]=-GXKOEF*GYKOEF*GYKOEF/ival*mukoef;
/*
			wi[0]=0;
*/
			wi[1]=1.0*mukoef;
		        gausmat4(1,mcmax,0);
			}
		/* Gravity potential per unit mass equation */
		else
			{
			wn[0]=8;
			wi[0]=4.0/1.5*pival*gamma*(ro[mcmax-1-ynumy]+ro[mcmax-ynumy]+ro[mcmax-1]+ro[mcmax])/4.0;
			/* Cross */
			/*       m1-1   m1    m1+1 */
			/* m2-1         G1         */
			/* m2     G0    G2     G4  */
			/* m2+1         G3         */
			/* d(dG/dX)/dX = ((G4-G2)/dX24-(G2-G0)/dX02)/(dX02/2+dX24/2)       */
			dx02=(gx[m1]-gx[m1-2])/2.0;
			dx24=(gx[m1+1]-gx[m1-1])/2.0;
			wn[1]=mcmax-ynumy;
			wi[1]= 2.0/dx02/(dx02+dx24);
			wn[2]=mcmax;
			wi[2]=-2.0/dx02/(dx02+dx24);
			wn[3]=mcmax+ynumy;
			wi[3]= 2.0/dx24/(dx02+dx24);
			wn[4]=mcmax;
			wi[4]=-2.0/dx24/(dx02+dx24);
			/* d(dG/dY)/dY = ((G3-G2)/dY23-(G2-G1)/dY12)/(dY12/2+dY23/2)       */
			dy12=(gy[m2]-gy[m2-2])/2.0;
			dy23=(gy[m2+1]-gy[m2-1])/2.0;
			wn[5]=mcmax-1;
			wi[5]= 2.0/dy12/(dy12+dy23);
			wn[6]=mcmax;
			wi[6]=-2.0/dy12/(dy12+dy23);
			wn[7]=mcmax+1;
			wi[7]= 2.0/dy23/(dy12+dy23);
			wn[8]=mcmax;
			wi[8]=-2.0/dy23/(dy12+dy23);
			if(ival>GYKOEF)
				{
				wn[0]=9;
				wi[0]-=GXKOEF*GYKOEF*GYKOEF/ival*mukoef*1e+18;
				wn[9]=mcmax;
				wi[9]=1.0*mukoef*1e+18;
				}
		        gausmat4(1,mcmax,0);
			}
		}
	/* Add ghost parameters */
	else
		{
		wn[0]=1;
		wi[0]=0;
		wn[1]=mcmax;
		wi[1]=1.0*mukoef;
	        gausmat4(1,mcmax,0);
		}
	/**/
	}
/* Solve Matrix ==================================================== */
if (printmod) printf("Number of positions in global matrix = %ld  Number of unknown = %ld \n",pos0cur,nodenum);
gausmat4(0,nodenum,0);
/*
gausmat3(0,nodenum-1,0);
gausmat4(0,nodenum,0);
*/
/* End Solve Matrix ==================================================== */
/**/
/* Reload GRAVITY POTENTIAL */
for (m1=0;m1<nodenum;m1++) 
	{
	gp[m1]=x[m1];
	}
/**/
/* New function included by Taras Gerya, compute magnitude of gravity acceleration */
/* Node Cycle */
for (m1=1;m1<xnumx-1;m1++)
for (m2=1;m2<ynumy-1;m2++)
	{
	/* Pos GP in sol0[] */
	mcmax=m1*ynumy+m2;
		{
		/**/
		/* gx-component */
		ival=(gp[mcmax+ynumy]+gp[mcmax+ynumy+1]-gp[mcmax]-gp[mcmax+1])/(gx[m1+1]-gx[m1-1]);
		/**/
		/* gy-component */
		ival1=(gp[mcmax+1]+gp[mcmax+ynumy+1]-gp[mcmax]-gp[mcmax+1])/(gy[m2+1]-gy[m2-1]);
		/**/
		/* Calculate magnitude of gravity acceleration for each node [m/s^2] */
		ga[mcmax]=pow(ival*ival+ival1*ival1,0.500);
		}
	}
/* End upgrade */
}
/**/
/**/
/**/
/* Clear New Solution */
/* P, Vx,Vy */
pos0cur=0;
mukoef=nubeg*strmin;
/**/
/* Add Matrix by vX-vY-Stokes, Continuity, Boundary, EPS, SIG, Equations */
/* Node  Cycle */
for (m1=0;m1<xnumx;m1++)
for (m2=0;m2<ynumy;m2++)
	{
	/* Pos P,Vx,Vy in sol0[] */
	mcmax=(m1*ynumy+m2)*3;
/*
printf("P-Vx-Vy Cycle %ld %ld    %ld",m1,m2,mcmax); getchar();
*/
	/**/
	/**/
	/**/
	/* Add Continuity equation for Cells ------------------------------------------------ */
	if(m1 && m2)
		{
/*
printf("Pr %ld %ld    %ld",m1,m2,bondm[mcmax+0]); getchar();
*/
		if(!bondm[mcmax+0]) 
			{
			/* Check Pressure in Cell calculation Y/N */ 
	                conterr(m1,m2,0);
			/* Rescale coefficients */
			for (n1=0;n1<=wn[0];n1++) wi[n1]*=mukoef;
	        	gausmat4(1,mcmax+0,0);
			}
		else
			{
			/* Add P-Boundary */
			pbonderr(mcmax+0,0);
			/* Rescale coefficients */
			for (n1=0;n1<=wn[0];n1++) wi[n1]*=mukoef;
	        	gausmat4(1,mcmax+0,0);
			}
		}
	/* Add ghost parameters */
	else
		{
		wn[0]=1;
		wi[0]=0;
		wn[1]=mcmax+0;
		wi[1]=mukoef;
	        gausmat4(1,mcmax+0,0);
		}
	/**/
	/**/
	/**/
	/* Add vX-Equations --------------------------------------------------- */
	if(m2<ynumy-1)
		{
/*
printf("Vx %ld %ld    %ld",m1,m2,bondm[mcmax+1]); getchar();
*/
		if(!bondm[mcmax+1]) 
			{
			/* Add vX-Stokes */
	                xstokserr(m1,m2,0);
			/**/
        		/* Add matrix */
			/* vX */
			gausmat4(1,mcmax+1,0);
/*
printf("Vx %ld %ld    %ld",m1,m2,bondm[mcmax+1]); getchar();
*/
			}
		else
			{
			/* Continuity Equation Vx boundary condition */
			if(bondv[bondm[mcmax+1]][1]>1e+30)
				{
				/* m1 m2 increment definition */
				m3=(bondn[bondm[mcmax+1]][0]-1-(mcmax+1))/3;
				dm1=dm2=0;
				if(m3>=ynumy) dm1=1;
				if(m3==1 || m3>ynumy) dm2=1;
/*
printf("Vx(%ld) P(%ld) boundary conditions for X=%ld Y=%ld Node",mcmax+1,bondn[bondm[mcmax+1]][0]-2,m1,m2); getchar();
*/
				if(m3<0 || m3>ynumy+1 || !(m1+dm1) || !(m2+dm2) || !bondm[bondn[bondm[mcmax+1]][0]-2])
					{
					printf("EXIT PROGRAM Inconsistent Vx(%ld) P(%ld) boundary conditions for X=%ld Y=%ld Node",mcmax+1,bondn[bondm[mcmax+1]][0]-2,m1,m2);
					exit(0);
					}
		                conterr(m1+dm1,m2+dm2,0);
				}
			else
				{
				/* Add vX Simple Boundary */
				xbonderr(mcmax+1,0);
				}
			/**/
			/* Vx boundary condition Add */
			/* Rescale coefficients */
			for (n1=0;n1<=wn[0];n1++) wi[n1]*=mukoef;
	        	gausmat4(1,mcmax+1,0);
			}
		}
	/* Add ghost parameters */
	else
		{
		wn[0]=1;
		wi[0]=0;
		wn[1]=mcmax+1;
		wi[1]=mukoef;
	        gausmat4(1,mcmax+1,0);
		}
	/**/
	/**/
	/**/
	/* Add vY-Equations --------------------------------------------------- */
	if(m1<xnumx-1)
		{
/*
printf("Vy %ld %ld    %ld",m1,m2,bondm[mcmax+2]); getchar();
*/
		if(!bondm[mcmax+2]) 
			{
			/* Add vX-Stokes */
        	        ystokserr(m1,m2,0);
			/**/
        		/* Add matrix */
			/* vY */
			gausmat4(1,mcmax+2,0);
/*
printf("Vy2 %ld %ld    %ld",m1,m2,bondm[mcmax+8]); getchar();
*/
			}
		else
			{
			/* Continuity Equation Vy boundary condition */
			if(bondv[bondm[mcmax+2]][1]>1e+30)
				{
				/* m1 m2 increment definition */
				m3=(bondn[bondm[mcmax+2]][0]-1-(mcmax+2))/3;
				dm1=dm2=0;
				if(m3>=ynumy) dm1=1;
				if(m3==1 || m3>ynumy) dm2=1;
				if(m3<0 || m3>ynumy+1 || !(m1+dm1) || !(m2+dm2) || !bondm[bondn[bondm[mcmax+2]][0]-3])
/*
printf("Vy(%ld) P(%ld) boundary conditions for X=%ld Y=%ld Node",mcmax+2,bondn[bondm[mcmax+2]][0]-3,m1,m2); getchar();
*/
					{
					printf("EXIT PROGRAM Inconsistent Vy(%ld) P(%ld) boundary conditions for X=%ld Y=%ld Node",mcmax+2,bondn[bondm[mcmax+2]][0]-3,m1,m2);
					exit(0);
					}
	        	        conterr(m1+dm1,m2+dm2,0);
				}
			/* Simple Vy boundary condition */
			else
				{
				/* Add vY Simple Boundary */
				ybonderr(mcmax+2,0);
				}
			/**/
			/* Vy boundary condition Add */
			/* Rescale coefficients */
			for (n1=0;n1<=wn[0];n1++) wi[n1]*=mukoef;
	        	gausmat4(1,mcmax+2,0);
			}
		}
	/* Add ghost parameters */
	else
		{
		wn[0]=1;
		wi[0]=0;
		wn[1]=mcmax+2;
		wi[1]=mukoef;
	        gausmat4(1,mcmax+2,0);
		}
	/**/
	/**/
	/**/
	}
/* End  Add Matrix By vX-vY-Stokes, Continuity Equations */
/**/
/**/
/**/
/* Solve Matrix ================================================ */
if (printmod) printf("Number of positions in global matrix = %ld  Number of unknown = %ld \n",pos0cur,nodenum3);
gausmat4(0,nodenum3,0);
/* Solve Matrix ================================================ */
/**/
/**/
/**/
/* Reload P, Vx, Vy Results */
/* Node  Cycle */
for (m1=0;m1<xnumx;m1++)
	{
	for (m2=0;m2<ynumy;m2++)
		{
		/* Pos P,Vx,Vy in sol0[] */
		mcmax0=(m1*ynumy+m2)*3;
		/* Pos in vx[], vy[], pr[], etc. */
		mcmax1=m1*ynumy+m2;
		/**/	
		/**/	
		/**/	
		/* Reload/Recalc P */	
		if(m1 && m2) 
			{
			/* Reload P */	
			pr[mcmax1]=x[mcmax0+0];
			}
		/**/	
		/**/	
		/**/	
		/* Reload Vx */	
		if(m2<ynumy-1)
			{ 
			vx[mcmax1]=x[mcmax0+1];
			}
		/**/	
		/**/	
		/**/	
		/* Reload Vy */	
		if(m1<xnumx-1)
			{
			vy[mcmax1]=x[mcmax0+2];
			}
/*	
printf("\n %ld %ld %ld %e %e",m1,m2,mcmax1,vx[mcmax1],vy[mcmax1]); getchar();
*/	
		}
	}
/* End Reload P, Vx, Vy Results */
/**/	
/**/	
/**/	
/* Recalc EPS, SIG Results */
/* Node  Cycle */
for (m1=1;m1<xnumx;m1++)
for (m2=1;m2<ynumy;m2++)
	{
	/* Pos in vx[], vy[], pr[], etc. */
	mcmax1=m1*ynumy+m2;
	/**/	
	/**/	
	/**/	
	/* Sxx,Exx */	
	sxx[mcmax1]=sxxcalc(m1,m2,0); exx[mcmax1]=eps[0];
	/**/	
	/* Syy,Eyy */	
	syy[mcmax1]=syycalc(m1,m2,0); eyy[mcmax1]=eps[0];
	/**/	
	/* Sxy,Exy */	
	if(m1<xnumx-1 && m2<ynumy-1)
		{
		sxy[mcmax1]=sxycalc(m1,m2,0); exy[mcmax1]=eps[0];
		}
	}
/* End Recalc EPS, SIG Results */
/**/	
/**/	
/**/	
/* Vx,Vy max-min definition */
for (m1=0;m1<xnumx;m1++)
for (m2=0;m2<ynumy;m2++)
	{
	/* Pos of Vx in sol0[] */
	m3=m1*ynumy+m2;
	/**/
	/* Min,Max Vx definition */
	if(m2<ynumy-1)
		{
		minvx=MINV(minvx,vx[m3]);
		maxvx=MAXV(maxvx,vx[m3]);
		}
	/* Min,Max Vy definition */
	if(m1<xnumx-1)
		{
		minvy=MINV(minvy,vy[m3]);
		maxvy=MAXV(maxvy,vy[m3]);
		}
	}
/* Max Vx,Vy Diff in Grid Calc */
vxmin=minvx;
vymin=minvy;
vxmax=maxvx;
vymax=maxvy;
maxvx-=minvx;
maxvy-=minvy;
maxvx=MAXV(maxvx,maxvy);
mindx=MINV(xstpx,ystpy);
/**/
/**/
/**/
/* Check Error */
/* Node  Cycle */
for (m1=0;m1<xnumx;m1++)
for (m2=0;m2<ynumy;m2++)
	{
	/* Pos of P,Vx,Vy in sol0[] */
	mcmax0=(m1*ynumy+m2)*3;
	/**/
	/**/
	/**/
	/* Check Continuity equation for Cells =========================== */
	if(m1 && m2) 
		{
		ival=conterr(m1,m2,1);
/*	
printf("\n %ld %ld   %e %e %e     %e ",m1,m2,ival,maxvx,mindx,ival/(maxvx/mindx)); getchar();
*/	
		contsum+=ival*ival;
       	        contnum+=1.0;
		/* Print Results */
		ival/=maxvx/mindx;
		if (printmod && ABSV(ival)>DIVVMIN)
			{
			printf("\n Large Continuity err at X=%ld Y=%ld:   Err=%e",m1,m2,ival);
			}
		}
	/**/
	/**/
	/**/
	/* Check vX-Equations for nodes =========================== */
	if(m2<ynumy-1)
		{
		if(!bondm[mcmax0+1]) 
			{
			/* Add vX-Stokes */
	                ival=xstokserr(m1,m2,1);
        	        stoksum+=ival*ival;
                	stoknum+=1.0;
	                /* Min,Max Nu value Calc */
        	        minnu=MINV(minnu,errbuf[9]);
	                maxnu=MAXV(maxnu,errbuf[9]);
        	        /* Min,Max Pr value Calc */
                	maxpr=MAXV(maxpr,errbuf[10]);
	                minpr=MINV(minpr,errbuf[11]);
			/* Print Results */
			ival/=maxvx/mindx/mindx;
			if (printmod && ABSV(ival)>STOKSMIN)
				{
				printf("\n Large X stokes err at X=%ld Y=%ld:   Err=%e",m1,m2,ival);
				}
	                }
		else
			{
			if(bondv[bondm[mcmax0+1]][1]<1e+30)
				{
				/* Add vX-Boundary */
				ival=xbonderr(mcmax0+1,1);
				bondsum+=ival*ival;
	                	bondnum+=1.0;
				}
			}
		}
	/**/
	/**/
	/**/
	/* Check vY-Equations for nodes =========================== */
	if(m1<xnumx-1)
		{
		if(!bondm[mcmax0+2]) 
			{
			/* Add vX-vY-Stokes */
	                ival=ystokserr(m1,m2,1);
        	        stoksum+=ival*ival;
                	stoknum+=1.0;
	                /* Min,Max Nu value Calc */
        	        minnu=MINV(minnu,errbuf[9]);
	                maxnu=MAXV(maxnu,errbuf[9]);
        	        /* Min,Max Pr value Calc */
	                maxpr=MAXV(maxpr,errbuf[10]);
        	        minpr=MINV(minpr,errbuf[11]);
			/* Print Results */
			ival/=maxvx/mindx/mindx;
			if (printmod && ABSV(ival)>STOKSMIN)
				{
				printf("\n Large Y stokes err at X=%ld Y=%ld:   Err=%e",m1,m2,ival);
				}
                	}
		else
			{
			if(bondv[bondm[mcmax0+2]][1]<1e+30)
				{
				/* Add vX-vY-Boundary */
				ival=ybonderr(mcmax0+2,1);
				bondsum+=ival*ival;
                		bondnum+=1.0;
				}
			}
		}
/*	
printf("\n %ld %ld   %e %e %e ",m1,m2,minnu,maxnu,errbuf[9]); getchar();
*/	
	}
stoksum=pow(stoksum/stoknum,0.5)/(maxvx/mindx/mindx);
contsum=pow(contsum/contnum,0.5)/(maxvx/mindx);
bondsum=pow(bondsum/bondnum,0.5)/maxvx;
/* End Check Error */
/**/
/**/
/**/
/* Print Results */
if (printmod)
	{
	 printf("\n KRUG %2d \n",m0+1);
	 printf("X-VELOCITY: min = %e max = %e \n",vxmin,vxmax);
	 printf("Y-VELOCITY: min = %e max = %e \n",vymin,vymax);
	 printf("PRESSURE: min = %e max = %e \n",minpr,maxpr);
	 printf("VISKOS: min = %e max = %e \n",minnu,maxnu);
	 printf("STOKES: num = %e err = %e \n",stoknum,stoksum);
	 printf("CONT : num = %e err = %e \n",contnum,contsum);
	 printf("BOUND V: num = %e err = %e \n",bondnum,bondsum);
	 }
/**/
/* Printf EPS, SIG Results */
/*
for (m1=100;m1<xnumx;m1++)
for (m2=0;m2<ynumy;m2++)
	{
	allintere(gx[m1],gy[m2]);
	ival=pow(0.5*(eps[6]*eps[6]+eps[8]*eps[8])+eps[4]*eps[4],0.5);
printf("%ld %ld %e %e   %e   %e %e %e   %e",m1,m2,gx[m1],gy[m2],eps[10],eps[6],eps[8],eps[4],ival);getchar();
	}
*/	
/**/	
/**/
/**/
/**/
}
/* End Solve XY-Stokes+Continuity equations by vX,vY,P mixed arbitrary order Finite Diff method */





/* Max Vx,Vy in nodes serch time step recalc */
void maxvelstep()
{
double maxvx=0,maxvy=0,ival;
long int m1,m2,m3;
/**/
/**/
/**/
/* Node  Cycle */
for (m1=0;m1<xnumx;m1++)
for (m2=0;m2<ynumy;m2++)
	{
	/* Pos of in Vx,Vy */
	m3=m1*ynumy+m2;
	/**/
	if(m2<ynumy-1)
		{
		ival=ABSV(vx[m3]);
		maxvx=MAXV(maxvx,ival);
		}
	if(m1<xnumx-1)
		{
		ival=ABSV(vy[m3]);
		maxvy=MAXV(maxvy,ival);
		}
	}
if (maxvx)
	{
	maxvx=(maxxystep*xstpx)/maxvx;
	if(printmod) printf("\n !!! MAX VALID TIME STEP FOR Vx-MARKER %e YEAR !!!\n",maxvx/3.15576e+7);
	timestep=MINV(maxvx,timestep);
	}
if (maxvy)
	{
	maxvy=(maxxystep*ystpy)/maxvy;
	if(printmod) printf("\n !!! MAX VALID TIME STEP FOR Vy-MARKER %e YEAR !!!\n",maxvy/3.15576e+7);
	timestep=MINV(maxvy,timestep);
	}
}
/* Max Vx,Vy in nodes serch time step recalc */




/* Weight of FD calculation for after Fornberg (1996) */
void fdweight(int n, int m, double xi)
/* n - maximal index 0-n */
/* m - required derivative order 0-m */
/* xi - derivation point coordinate */
{
/* Counters */
int i,j,k,mn;
double c1,c2,c3,c4,c5,kk;
/**/
/**/
/**/
c1=1.0;
c4=xn[0]-xi;
for(k=0;k<=m;k++)
	{
	for(j=0;j<=n;j++)
		{
		cn[j][k]=0;
		}
	}
/**/
cn[0][0]=1.0;
for(i=1;i<=n;i++)
	{
	mn=i;if(mn>m) mn=m;
	c2=1.0;
	c5=c4;
	c4=xn[i]-xi;
	for(j=0;j<i;j++)
		{
		c3=xn[i]-xn[j];
		c2*=c3;
		for(k=mn;k>0;k--)
			{
			kk=(double)(k);
			cn[i][k]=c1*(kk*cn[i-1][k-1]-c5*cn[i-1][k])/c2;
			}
		cn[i][0]=-c1*c5*cn[i-1][0]/c2;
		for(k=mn;k>0;k--)
			{
			kk=(double)(k);
			cn[j][k]=(c4*cn[j][k]-kk*cn[j][k-1])/c3;
			}
		cn[j][0]=c4*cn[j][0]/c3;
		}
	c1=c2;
	}
/*
for(i=0;i<=n;i++)printf("FD %d %d %e %d %e %e\n",n,m,xi,i,xn[i],cn[i][m]);getchar();
*/
}
/* Weight of FD calculation after Fornberg (1996) */



/* Left side or Value for Sxx  Equation */ 
/* Sxx=2Nu*Exx, Exx=dVx/dX */
double sxxcalc(long int m1, long int m2, double ynval)
/* m1,m2 - node X,Y number */
/* ynval - Val Sxx Calc Y(0)/N(koefficient) */
{
/* Exx horisontal position */
double xi=(gx[m1-1]+gx[m1])/2.0,leftsxx=0,nueff;
long int m1min,m1max,m3,v[4];
int n1,n;
/**/
/**/
/**/
/* Staggered Nodes num */
/*   [0]                [2] */
/*   Nu0                Nu2 */
/*                          */
/*   Vx0    Sxx3,Exx3   Vx2 */
/*                          */
/*   [1]                [3] */
/*   Nu1                Nu3 */
/*                          */
v[0]=(m1-1)*ynumy+(m2-1);v[1]=v[0]+1;
v[2]=v[0]+ynumy;v[3]=v[2]+1;
/**/
/**/
/**/
/* Effective viscosity calc */
nueff=nd[v[3]];
/**/
/**/
/**/
/* Calc, Check Fd limits */
m1min=m1-1-stoksfd; if(m1min<0) m1min=0;
m1max=m1+stoksfd; if(m1max>xnumx-1) m1max=xnumx-1;
/**/
/* Load distances to xn[] */
for (m3=m1min;m3<=m1max;m3++)
	{
	xn[m3-m1min]=gx[m3];
	}
/**/
/* Calc maximal position in xn[] */
n=(int)(m1max-m1min);
/**/
/* Calc Vx coefficients for EPSxx */
fdweight(n,1,xi);
/**/
/* Return Sxx,Exx val ----------------------------*/
if(ynval==0)
	{
	/* Exx=dVx/dX */
	/* Add Vx with koefficients */
	for (m3=m1min;m3<=m1max;m3++)
		{
		leftsxx+=vx[m3*ynumy+(m2-1)]*cn[m3-m1min][1];
		}
	/**/
	/* Save Exx */
	eps[0]=leftsxx;
	/**/
	/* Calc Sxx=2Nu*Exx */
	leftsxx*=2.0*nueff;
	/**/
	return leftsxx;
	}
/**/
/**/
/**/
/* Add Coefficients for left parts of Sxx ----------------*/
/*  0(P) 1(Vx)  2(Vy)  */
/* Sxx=2Nu*Exx, Exx=dVx/dX */
/* Add Vx with koefficients */
for (m3=m1min;m3<=m1max;m3++)
	{
	wn[wn[0]+1+m3-m1min]=(m3*ynumy+(m2-1))*3+1;
	wi[wn[0]+1+m3-m1min]=ynval*2.0*nueff*cn[m3-m1min][1];
	}
/**/
/* Add total Num of lines */
wn[0]+=n+1;
/**/
/**/
/**/
/*
for(n1=0;n1<=wn[0];n1++)printf("Exx %e %ld \n",wi[n1],wn[n1]);getchar();
*/
return 0;
}
/* Left side or Value for Sxx  Equation */ 




/* Left side or Value for Syy  Equation */ 
/* Syy=2Nu*Eyy, Eyy=dVy/dY */
double syycalc(long int m1, long int m2, double ynval)
/* m1,m2 - node X,Y number */
/* ynval - Val Syy Calc Y(0)/N(koefficient) */
{
/* Eyy vertical position */
double xi=(gy[m2-1]+gy[m2])/2.0,leftsyy=0,nueff;
long int m2min,m2max,m3,v[4];
int n1,n;
/**/
/**/
/**/
/* Staggered Nodes num */
/*   [0]       Vy0      [2] */
/*   Nu0                Nu2 */
/*                          */
/*          Sxx3,Exx3       */
/*                          */
/*   [1]       Vy1      [3] */
/*   Nu1                Nu3 */
/*                          */
v[0]=(m1-1)*ynumy+(m2-1);v[1]=v[0]+1;
v[2]=v[0]+ynumy;v[3]=v[2]+1;
/**/
/**/
/**/
/* Effective viscosity calc */
nueff=nd[v[3]];
/**/
/**/
/**/
/* Calc, Check Fd limits */
m2min=m2-1-stoksfd; if(m2min<0) m2min=0;
m2max=m2+stoksfd; if(m2max>ynumy-1) m2max=ynumy-1;
/**/
/* Load distances to xn[] */
for (m3=m2min;m3<=m2max;m3++)
	{
	xn[m3-m2min]=gy[m3];
	}
/**/
/* Calc maximal position in xn[] */
n=(int)(m2max-m2min);
/**/
/* Calc Vy coefficients for EPSyy */
fdweight(n,1,xi);
/**/
/* Return Syy,Eyy val ----------------------------*/
if(ynval==0)
	{
	/* Eyy=dVy/dY */
	/* Add Vy with koefficients */
	for (m3=m2min;m3<=m2max;m3++)
		{
		leftsyy+=vy[(m1-1)*ynumy+m3]*cn[m3-m2min][1];
		}
	/**/
	/* Save Eyy */
	eps[0]=leftsyy;
	/**/
	/* Calc Syy=2Nu*Eyy */
	leftsyy*=2.0*nueff;
	/**/
	return leftsyy;
	}
/**/
/**/
/**/
/* Add Coefficients for left parts of Syy ----------------*/
/*  0(P) 1(Vx)  2(Vy)  */
/* Sxx=2Nu*Exx, Exx=dVx/dX */
/* Add Vy with koefficients */
for (m3=m2min;m3<=m2max;m3++)
	{
	wn[wn[0]+1+m3-m2min]=((m1-1)*ynumy+m3)*3+2;
	wi[wn[0]+1+m3-m2min]=ynval*2.0*nueff*cn[m3-m2min][1];
	}
/**/
/* Add total Num of lines */
wn[0]+=n+1;
/**/
/**/
/**/
/*
for(n1=0;n1<=wn[0];n1++)printf("Eyy %e %ld\n",wi[n1],wn[n1]);getchar();
*/
return 0;
}
/* Left side or Value for Syy  Equation */ 





/* Left side or Value for Sxy  Equation */ 
/* Sxy=2Nu*Exy, Exy=1/2(dVx/dY+dVy/dX) */
double sxycalc(long int m1, long int m2, double ynval)
/* m1,m2 - node X,Y number */
/* ynval - Val Syy Calc Y(0)/N(koefficient) */
{
/* Exy position */
double xi,leftsxy=0,nueff;
long int m1min,m1max,m2min,m2max,m3;
int n1,nx,ny;
/**/
/**/
/**/
/* Staggered Nodes num */
/*  [0]                [2]            */
/*                                    */
/*                     Vx2            */
/*                                    */
/*  [1]     Vy1        [3]       Vy3  */
/*                   Exy3,Nu3         */
/*                                    */
/*                     Vx3            */
/*                                    */
/* Effective viscosity calc */
nueff=nu[m1*ynumy+m2];
/**/
/**/
/**/
/* dVx/dY */
xi=gy[m2];
/* Calc, Check Fd limits */
m2min=m2-1-stoksfd;
if(m2min<0) m2min=0;
m2max=m2+stoksfd;
if(m2max>ynumy-2) m2max=ynumy-2;
/**/
/* Load distances to xn[] */
for (m3=m2min;m3<=m2max;m3++)
	{
	xn[m3-m2min]=(gy[m3]+gy[m3+1])/2.0;
	}
/**/
/* Calc maximal position in xn[] */
nx=(int)(m2max-m2min);
/**/
/* Calc Vx coefficients for EPSxy */
fdweight(nx,1,xi);
/* Reload coefficients to cn[] */
for (m3=0;m3<=nx;m3++)
	{
	cn[m3][2]=cn[m3][1];
	}
/**/
/**/
/**/
/* dVy/dX */
xi=gx[m1];
/* Calc, Check Fd limits */
m1min=m1-1-stoksfd;
if(m1min<0) m1min=0;
m1max=m1+stoksfd;
if(m1max>xnumx-2) m1max=xnumx-2;
/**/
/* Load distances to xn[] */
for (m3=m1min;m3<=m1max;m3++)
	{
	xn[m3-m1min]=(gx[m3]+gx[m3+1])/2.0;
	}
/**/
/* Calc maximal position in xn[] */
ny=(int)(m1max-m1min);
/**/
/* Calc Vy coefficients for EPSxy */
fdweight(ny,1,xi);
/**/
/* Return Sxy,Exy val ----------------------------*/
if(ynval==0)
	{
	/* Exy=1/2(dVx/dY+dVy/dX)=0 */
	/* 1/2dVx/dY */
	/* Add Vx with koefficients */
	for (m3=m2min;m3<=m2max;m3++)
		{
		leftsxy+=vx[m1*ynumy+m3]*cn[m3-m2min][2]/2.0;
		}
	/**/
	/* 1/2dVy/dX */
	/* Add Vy with koefficients */
	for (m3=m1min;m3<=m1max;m3++)
		{
		leftsxy+=vy[m3*ynumy+m2]*cn[m3-m1min][1]/2.0;
		}
	/**/
	/* Save Eyy */
	eps[0]=leftsxy;
	/**/
	/* Calc Sxy=2Nu*Exy */
	leftsxy*=2.0*nueff;
	/**/
	return leftsxy;
	}
/**/
/**/
/**/
/* Add Coefficients for left parts of Syy ----------------*/
/*  0(P) 1(Vx)  2(Vy)  */
/* Sxy=2Nu*Exy, Exy=1/2(dVx/dY+dVy/dX) */
/**/
/* 1/2dVx/dY */
/* Add Vx with koefficients */
for (m3=m2min;m3<=m2max;m3++)
	{
	wn[wn[0]+1+m3-m2min]=(m1*ynumy+m3)*3+1;
	wi[wn[0]+1+m3-m2min]=ynval*2.0*nueff*cn[m3-m2min][2]/2.0;
	}
/**/
/* Add total Num of lines */
wn[0]+=nx+1;
/**/
/* 1/2dVy/dX */
/* Add Vy with koefficients */
for (m3=m1min;m3<=m1max;m3++)
	{
	wn[wn[0]+1+m3-m1min]=(m3*ynumy+m2)*3+2;
	wi[wn[0]+1+m3-m1min]=ynval*2.0*nueff*cn[m3-m1min][1]/2.0;
	}
/**/
/* Add total Num of lines */
wn[0]+=ny+1;
/**/
/**/
/**/
/*
for(n1=0;n1<=wn[0];n1++)printf("Exy %e %ld \n",wi[n1],wn[n1]);getchar();
*/
return 0;
}
/* Left side or Value for Sxy  Equation */ 


