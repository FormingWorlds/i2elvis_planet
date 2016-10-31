/* Move markers by using Runge-Kutta method */
void movemark()
{
/* Vx, Vy buffer */
double vx0,vx1,vx2,vx3,vx4,vy0,vy1,vy2,vy3,vy4,xnew,ynew,xold=0,yold=0;
/**/
long int mm1,marknum1;
/* Erosion-Sedimentation Y/N */
int n1;
int mm2;
/* Nonstability for immobile markers */
double xnonstab=0.50,ynonstab=0.60,mnu,mpb,mtk;
int nonstab=100;
/**/
/**/
/* Save number of markers */
marknum1=marknum;
/**/
/**/
/**/
/* Move markers */
for (mm1=0;mm1<marknum;mm1++)
{
/* Marker type */
mm2=(int)markt[mm1]; if (mm2>=100) mm2-=100;
if( ((markx[mm1]>=0 && marky[mm1]>=0 && (double)(markx[mm1])<=xsize && (double)(marky[mm1])<=ysize) || outgrid!=1) && !markim[mm2] )
	{
	/* Random displacement for immobile markers */
	if(markt[mm1]>=100)
		{
		/* X,Y save */
		xold=(double)(markx[mm1]);
		yold=(double)(marky[mm1]);
		/* Random nonstability Set on X,Y */
		if(nonstab>0)
			{
			markx[mm1]+=(float)(rand() % (nonstab*2+1) - nonstab)/((float)(nonstab))*(float)(xnonstab*mardx);
			marky[mm1]+=(float)(rand() % (nonstab*2+1) - nonstab)/((float)(nonstab))*(float)(ynonstab*mardy);
			}
		}
	/**/
	/**/
	/**/
	/* Calc/Save strain rate and stress invariants, viscosity */
	allintere((double)(markx[mm1]),(double)(marky[mm1]));
	markeii=pow(0.5*(eps[6]*eps[6]+eps[8]*eps[8])+eps[4]*eps[4],0.5);
	marksii=pow(0.5*(eps[7]*eps[7]+eps[9]*eps[9])+eps[5]*eps[5],0.5);
	/* Correct strain rate for numerical diffusion */
	markrii=1.0;
	if(0==0 && markv[mm1]>0) markrii*=pow(0.5*marksii/markeii/markv[mm1],0.87);
	mpb=eps[10]*1e-5;
	mtk=(double)(markk[mm1]);
	/* Save pressure */
	markp[mm1]=eps[10];
	/* Compute viscosity */
	if(markt[mm1]<50 && mtk>0)
		{
		if(markt[mm1]<20)
			{
			mnu=viscalc(mtk,mpb,(double)(markx[mm1]),(double)(marky[mm1]),mm1,mm2);
			}
		else
			{
			mnu=viscalc(mtk,mpb,(double)(markx[mm1]),(double)(marky[mm1]),mm1,mm2-20);
			}
		markv[mm1]=mnu;
		}
	/**/
	/**/
	/**/
	/* Motion Calc ///////////////////////////////// */
	/**/
	if(markmod==1)
		{
		/* Vx, Vy simple calc */
		allinterv((double)(markx[mm1]),(double)(marky[mm1]));
		vx0=eps[11]; vy0=eps[12];
		}
	else
		{
		/* Vx, Vy 4 Runge-Kutta coef calc */
		allinterv((double)(markx[mm1]),(double)(marky[mm1]));
		vx1=eps[11]; vy1=eps[12];
		/**/
		allinterv((double)(markx[mm1])+vx1*timestep/2.0,(double)(marky[mm1])+vy1*timestep/2.0);
		vx2=eps[11]; vy2=eps[12];
		/**/
		allinterv((double)(markx[mm1])+vx2*timestep/2.0,(double)(marky[mm1])+vy2*timestep/2.0);
		vx3=eps[11]; vy3=eps[12];
		/**/
		allinterv((double)(markx[mm1])+vx3*timestep,(double)(marky[mm1])+vy3*timestep);
		vx4=eps[11]; vy4=eps[12];
		/**/
		/* Vx,Vy, EpsXX, EpsYY, EpsXY calc after Runge-Kutta */
		vx0=(vx1+2.0*vx2+2.0*vx3+vx4)/6.0;
		vy0=(vy1+2.0*vy2+2.0*vy3+vy4)/6.0;
		}
	/**/
	/* Orthogonal motion only */
	if (outgrid==2)
		{
		if(markx[mm1]<0 || (double)(markx[mm1])>xsize) vy0=0;
		if(marky[mm1]<0 || (double)(marky[mm1])>ysize) vx0=0;
		}
	/**/
	/**/
	/**/
	/* Normal/Immobile markers */
	if(markt[mm1]<100)
		{
		/* Normal markers */
		/* X,Y calc after Runge-Kutta */
		markx[mm1]+=(float)(timestep*vx0);
		marky[mm1]+=(float)(timestep*vy0);
		/**/
		}
	else
		{
		/* Immobile markers */
		/* X,Y calc after Runge-Kutta */
		xnew=(double)(markx[mm1])+timestep*vx0;
		ynew=(double)(marky[mm1])+timestep*vy0;
		/**/
		/* X,Y reset for immobile marker */
		markx[mm1]=(float)(xold);
		marky[mm1]=(float)(yold);
/**/
		/* Check new position, add marker */
		if(xnew>=0 && ynew>=0 && xnew<=xsize && ynew<=ysize)
			{
			/* Type save */
			markt[marknum1]=markt[mm1]-100;
			/* X,Y calc after Runge-Kutta */
			markx[marknum1]=(float)(xnew);
			marky[marknum1]=(float)(ynew);
			/* Temperature, Viscosity */
			markk[marknum1]=0;
			markv[marknum1]=0;
			/* Add additional markers counter */
			marknum1++;
			}
		/**/
		}
	/**/
	/**/
	/**/
	/* Out of grid marker reset */
	if(markx[mm1]<0 || marky[mm1]<0 || (double)(markx[mm1])>xsize || (double)(marky[mm1])>ysize)
		{
		markk[mm1]=0;
		}
	/* Motion Calc ///////////////////////////////// */
	}
}
/**/
/**/
/**/
/* Mark num */
if(marknum1>MAXMRK) {printf("Space   out in markx[]"); exit(0);}
/**/
/**/
/**/
/* Reset additional markers */
mm1=0;
while(marknum1>marknum && mm1<marknum)
	{
	/* Reload marker */
	if((markx[mm1]<0 || marky[mm1]<0 || (double)(markx[mm1])>xsize || (double)(marky[mm1])>ysize) && markt[mm1]<100)
		{
		/* Decrease additional markers counter */
		marknum1--;
		/* Type save */
		markt[mm1]=markt[marknum1];
		/* Temperature, Viscosity */
		markk[mm1]=0;
		markv[mm1]=0;
		/* X,Y reload */
		markx[mm1]=markx[marknum1];
		marky[mm1]=marky[marknum1];
		}
	/* Increase markers counter */
	mm1++;
	}
printf("\n Number of markers: OLD = %ld NEW = %ld \n",marknum,marknum1-1);
/* Set new marker number */
marknum=marknum1;
}
/* End Move markers by using Runge-Kutta method */




/* ro[],nu[] recalc after marker positions */
void ronurecalc()
{
/* Counters */
long int m1,m2,m3,m1min,m1max,m2min,m2max;
int mm2,yn=3,mm3,n1,n2;
int mag;                                    /* Greg: local marker magnetization parameter */
int n_count,no_iterate,no_iterate_tot=100;  /* Greg: Iteration parameters */
long int mm1;
double swt,swt1,swt2,celdx,celdy,ddx,ddy;
double mag_time,mga,mnu,mro,mcp,mkt,mht,mbb,mwa,dmwa,mro0;
double anu,aro,acp,akt,aht,abb,markwlev;
double wnu,wro,wcp,wkt,wht,wbb,dywa;
/* TD Database variables */
 double H0,H1,H2,H3,W0,W1,W2,W3,R0,R1,R2,R3,n,e,hydryl;
/*Newly added parameters */
double a_fact,c_fact,cont_r,decay_heat_al,decay_heat_fe,decay_heat_k,decay_heat_u,decay_heat_uu,decay_heat_th;
double eta_fluid,eta_fluid_fe,funct,funct_drv,funct_error,grain_corr,grain_guess,grain_old,grain_stress,heat_flux;
double init_al,init_fe,init_k,init_ka,init_u,init_uu,init_th;
double mht_al,mht_fe,mht_k,mht_ka,mht_u,mht_uu,mht_th,p_ref,por_orig,por_change,por_ref1,por_ref2;
double timesumA,xmelt,xmelt_fe;
/* As fit was done by Schwenn & Goetze (1978) using old fashioned non-SI units, apply these units for calculation */
double E_act=85000.000;        /* Apparent activation energy for sintering proceses [cal/mol] [Schwenn & Goetze, Tectonophysics, 48, 41-60 (1978)] */
double fact_a=4.000e-5;        /* Factor for porosity change calculation due to sintering [cm^3/(bar^(3/2)*s)] [Schwenn & Goetze, Tectonophysics, 48, 41-60 (1978)] */
double grain_orig=1.000e-6;    /* Assumed grain size [m] [Henke et al., Astronomy & Astrophysics, 537, A45 (2012)] */
double pivalue=3.141592654;    /* Define Pi */
double r_gas=1.9872;           /* Define universal gas constant [cal/(K mol)], for Schwenn & Goetze (1978) calculation */
/**/
double m_mars=6.4185e+23;      /* Present-day mass of Mars: 6.4185e23 kg [Lodders & Fegley, Planetary Scientist's Companion (1998)] */
/**/
/* RO, NU equations var */
double mpb=1.0,mtk=300.0,numax=0,numin=0,radmax;
/**/
/**/
printf("\n Number of nodes = %ld  Number of markers = %ld \n",nodenum,marknum);
/**/
/**/
/**/
/* Hydration front progress */
printf("timestep %e timesum %e densimod %d \n",timestep,timesum,densimod);
/**/
/**/
/**/
/* ADD MARKERS TO THE v-CELLS ========================== */
/* Clear ro[],nu[] wt */
for (m1=0;m1<nodenum;m1++)
	{
	ro0[m1]=0;
	et0[m1]=0;
	nu0[m1]=0;
	nd0[m1]=0;
	cp0[m1]=0;
	kt0[m1]=0;
	ht0[m1]=0;
	tk0[m1]=0;
	po0[m1]=0;
	sol0[m1]=0;
	sol1[m1]=0;
	sol0[nodenum+m1]=0;
	sol1[nodenum+m1]=0;
	}
/**/
/* Maximum planetary radius */
/*radmax=0;*/
/*for (mm1=0;mm1<marknum;mm1+=gridmod)
       {
       celdx=markx[mm1]-xsize/2.0;
       celdy=marky[mm1]-ysize/2.0;
       celdy=pow(celdx*celdx+celdy*celdy,0.5);
       if(markt[mm1]>1)
               {
               radmax=MAXV(radmax,celdy);
               }
       }
radmax+=(xsize/2.0-radmax)*2.0/3.0; */
/**/
/* Use a sticky air layer with constant outer radius, by Greg (last update: 26/09/2011) */
radmax=xsize/2.0*0.95;
/**/
/* Add ro[] nu[] etc. using selected markers */
for (mm1=0;mm1<marknum;mm1+=gridmod)
{
/* Marker type */
mm2=(int)markt[mm1]; if (mm2>=100) mm2-=100;
/* Check markers out of grid */
/**/
if(markx[mm1]>0 && marky[mm1]>0 && (double)(markx[mm1])<xsize && (double)(marky[mm1])<ysize && markk[mm1]>0 && markt[mm1]<50)
	{
	/* Remove Plumes */
	/* P, T parameters calc */
	allintere((double)(markx[mm1]),(double)(marky[mm1]));
	mpb=eps[10]*1.0000e-5;
	mtk=(double)(markk[mm1]);
/**/
	/* Newly added by Greg */
	mga=eps[49];      /* gravity acceleration [m/s^2] on the marker level */
        xmelt=eps[21];    /* silicate melt fraction [0-1] on the marker level */
        xmelt_fe=eps[50]; /* iron melt fraction [0 or 1] on the marker level */
/**/
	/* Temperature reset for sticky air to simulate impact-induced greenhouse atmosphere or space */
	/* by Greg (last modified: 05/03/2011) */
	mwa=0;
	/**/
	if((M_init+M_acc)<=(0.100*m_mars))
		{
		/* No impact-induced atmosphere, T is set to equilibrium temperature of fast rotating body at current Mars distance */
                /* Using assumptions: */
		/* Present day mean surface temperature [K] [Lodders & Fegley, The Planetary Scientist's Companion (1998)] */
		if(mm2<2) mtk=markk[mm1]=290.000;
		}
	if((M_init+M_acc)>(0.100*m_mars))
		{
		/* Impact-induced atmosphere when M>0.1*M_mars [Tyburczy et al., EPSL, 80, 201-207 (1986) (Fig. 3)] */
		/* Estimation of temperature of impact atmosphere [Abe, Earth Moon Planets, 108, 9-14 (2011) (Fig. 1)] */
        	/* if(mm2<2) mtk=markk[mm1]=2000.000; */
		if(mm2<2) mtk=markk[mm1]=290.000;
		}
	/**/
	/* Marker Properties */
	/* Viscosity calc */
	mnu=markv[mm1];
	if(mnu<=0)
		{
		markeii=pow(0.5*(eps[6]*eps[6]+eps[8]*eps[8])+eps[4]*eps[4],0.5);
		marksii=pow(0.5*(eps[7]*eps[7]+eps[9]*eps[9])+eps[5]*eps[5],0.5);
		markrii=1.0;
		if(markt[mm1]<20)
			{
			mnu=viscalc(mtk,mpb,(double)(markx[mm1]),(double)(marky[mm1]),mm1,mm2);
			}
		else
			{
			mnu=viscalc(mtk,mpb,(double)(markx[mm1]),(double)(marky[mm1]),mm1,mm2-20);
			}
		}
	celdx=markx[mm1]-xsize/2.0;
	celdy=marky[mm1]-ysize/2.0;
	celdy=pow(celdx*celdx+celdy*celdy,0.5);
	/**/
	/* Viscosity of outer stabilizing material is limited to 1e22 Pa s to ensure good pressure solution */
        if(celdy>radmax) mnu=1.000e+22;
	/**/
	/**/
	/* Compute porosity of solid silicates taking into account for cold pressing, based on [Henke et al., Astronomy & Astrophysics, 537, A45 (2012)] */
	/* added by Greg (last modified: 24/02/2014) */
	/**/
	/* No porosity for sticky air material, iron material and molten silicates */
	if(mm2==0 || mm2>=7)
		{
		markpor[mm1] = 0.000;
		}
	/**/
	/* Effect of cold pressing and sintering on solid silicate porosity [non-dim.] */
	if((mm2==5 || mm2==6) && (markpor[mm1]>0.000 && markpor[mm1]<=0.850))
		{
		/* Reference pressure [bar] */
		p_ref = 0.130;
		/**/
		/* Effect of cold pressing by self-gravity at low temperatures */
		if(markk[mm1]<700.000)
			{
			/* Assume here that porosity can only decrease, based on Henke et al. (2012) eq. (15) */
			markpor[mm1] = (float)(MINV((double)(markpor[mm1]),(0.420+0.460*pow((pow((mpb/p_ref),1.720)+1.000),(-1.000)))));
			/**/
			/* Introduce artifical upper cut-off for porosity to make sure Henke et al. eq. (15) can be used at all times */
			markpor[mm1] = (float)(MINV((double)(markpor[mm1]),0.850));
			}
		/**/
		/* Effect of hot sintering at higher temperatures */
		if(markk[mm1]>=700.000)
			{
			/* Reset parameters before each calculation to zero */
			c_fact       = 0.000;
			cont_r       = 0.000;
			funct        = 0.000;
			funct_drv    = 0.000;
			funct_error  = 1.000;    /* Set on purpose to large value */
			grain_corr   = 0.000;
			grain_guess  = 0.000;
			grain_old    = 0.000;
			grain_stress = 0.000;
			n_count      = 1;
			por_change   = 0.000;
			por_orig     = 0.000;
			/**/
			/* Estimate original porosity before onset of hot sintering, based on Henke et al. (2012) eq. (15) */
			por_orig     = MINV(por_init,(0.420+0.460*pow((pow((mpb/p_ref),1.720)+1.000),(-1.000))));
			/**/
			/* Introduce artificial upper-cut-off for porosity to make sure Henke et al. eq. (15) can be used at all times */
			por_orig     = MINV(por_orig,0.850);
			/**/
			/* Compute radius of contact area of overlaping spherical grains [m] */
			/**/
			/* Initial guess for corrected grain size [m] */
			grain_guess  = grain_orig;
			/**/
			/* Repeat calculation, if necessary */
			if(n_count>1) recalculate: printf("Repeat Newton-Raphson calculation \n");
			/**/
			/* For this purpose use Newton-Raphson method to compute corrected grain size [m] */
			for(no_iterate=1;no_iterate<no_iterate_tot;no_iterate=no_iterate+1)
				{
				/* Save old corrected grain size [m] */
				grain_old   = grain_guess;
				/**/
				/* Function to solve, assuming 8 as number of grain contacts, based on Henke et al. (2012) eq. (34)-(37) */
				funct       = pow(grain_guess,3.000)+pow(grain_orig,3.000)-3.000*pow(grain_guess,2.000)*grain_orig*pow(((1.000-por_orig)/(1.000-(double)(markpor[mm1]))),(1.000/3.000))+pow(grain_orig,3.000)*(1.000-por_orig)/(1.000-(double)(markpor[mm1]));
				/**/
				/* Derivative of function to solve */
				funct_drv   = 3.000*pow(grain_guess,2.000)-6.000*grain_guess*grain_orig*pow(((1.000-por_orig)/(1.000-(double)(markpor[mm1]))),(1.000/3.000));
				/**/
				/* Compute new value for corrected grain size [m] */
				grain_guess = grain_old-(funct/funct_drv);
				}
			/**/
			/* Check whether the solution is varying significantly */
			funct_error = ABSV(grain_guess-grain_old);
			/**/
			/* In case result is unstable or unphysical, repeat with larger guess value of grain size */
			if((funct_error>5.000e-2) || (grain_guess<=grain_orig))
				{
				n_count     = n_count+1;
				grain_guess = 1.250*(double)(n_count)*grain_orig;
			 	goto recalculate;
				}
			/**/
			/* Save corrected grain size [m], when larger than original grain size */
			if(grain_guess>grain_orig) grain_corr = grain_guess;
			/**/
			/* Compute radius of each of the 8 individual contact areas [m] */
			/* based on Henke et al. (2012), eq. (37) */
			cont_r       = sqrt(pow(grain_corr,2.000)-pow(grain_orig,2.000)*pow(((1.000-por_orig)/(1.000-(double)(markpor[mm1]))),(2.000/3.000)));
			/**/
			/* Compute average cross-section of the cell occupied by one grain unit [m^2] */
			/* based on Henke et al. (2012) eq. (41) */
			c_fact       = 2.000*sqrt(3.000)*(pow(grain_corr,2.000)-pow(cont_r,2.000));
			/**/
			/* Compute effective stress on the contact faces of two grains [bar] */
			/* based on Henke et al. (2012) eq. (40) */
			grain_stress = c_fact*mpb/(pivalue*pow(cont_r,2.000));
			/**/
			/* Compute change of porosity [non-dim.] (based on Henke et al. (2012) eq. (39)+(44)) */
			if(markpor[mm1]>0.100)
				{ 
				/* When using Schwenn & Goetze (1978) model, grain size unit is cm, stress is in bar and temperature in Celsius */
				por_change = fact_a*(1.000-(double)(markpor[mm1]))*pow(grain_stress,1.500)/pow((100.000*grain_corr),3.000)*exp(-E_act/(r_gas*mtk))*timestep;
				}
			/**/
			/* Based on Henke et al. (2012) eq. (43) */
			if(markpor[mm1]<=0.100)
				{
				/* When using Schwenn & Goetze (1978) model, grain unit is cm, stress is in bar and temperature in Celsius */
				por_change = 10.000*(double)(markpor[mm1])*fact_a*(1.000-(double)(markpor[mm1]))*pow(grain_stress,1.500)/pow((100.000*grain_corr),3.000)*exp(-E_act/(r_gas*mtk))*timestep;
				}
			/**/
			/* Assume here that porosity can only decrease */
			if(por_change>0.000) markpor[mm1] = (float)((double)(markpor[mm1])-por_change);
			}
		}
	/**/
	/* Make sure porosity stays within the range 0 to 100% */
	if(markpor[mm1]<0.000) markpor[mm1] = 0.000;
	if(markpor[mm1]>1.000) markpor[mm1] = 1.000;
	/**/
	/**/
	mro0=mro=dencalc(mtk,mpb,(double)(markpor[mm1]),(double)(markx[mm1]),(double)(marky[mm1]),mm2);
	mbb=eps[20];
	mcp=markcp[mm2];
	mkt=(markkt[mm2]+markkf[mm2]/(mtk+77.0))*exp(markkp[mm2]*mpb);
	/**/
	/* Thermodynamic database use for ro, Cp */
	if (densimod==2)
	if(mm2==5 || mm2==6 || mm2==11 || mm2==12 || mm2==25 || mm2==26 || mm2==31 || mm2==32)
		{
		/* Compute TD variables */
		eps[47]=mkt;
		eps[48]=mnu;
		tdbasecalc(mtk,mpb,mm2,mm1);
		mro=eps[41];
		mcp=eps[43];
		mbb=eps[44];
		mkt=eps[47];
		mnu=eps[48];
		/**/
		}
	/**/
	/**/
	/* Modify thermal conductivity to better simulate the cooling of a magma ocean/iron diapirs */
	/* by Greg (last modified: 30/10/2013) */
	/**/
	/* Only valid for molten silicates, not used for very first timestep */
        if((mm2==25 || mm2==26) && (timesum>core_form_time))
		{
		/* Assume crystal mush between 40 and 60 % of Si melt */
		if(xmelt>=0.4000 && xmelt<=0.6000)  /* Use rheological boundary value for silicate melts as suggested by [Solomatov, Treatise on Geophysics Vol. 9, 91-119 (2007)] */
			{
			/* Assume immediate drop of viscosity */
			eta_fluid = 100.000;
			/* Simplified assumption that viscosity of crystal mush decreases exponentially with melt fraction [width of transition zone based on Lejeune and Richet, JGR, 100, 1995] */
			/* eta_fluid = markn0[mm2]*exp((-172.690)*(xmelt-0.400)); */
			/* Simplified assumption that viscosity of crystal mush decreases linearly with melt fraction [width of transition zone based on Lejeune and Richet, JGR, 100, 1995] */
			/* eta_fluid = mnu-((mnu-100.0000)/0.2000)*(xmelt-0.4000); */
			/**/
			/* Apply the soft turbulence model [Solomatov, Treatise on Geophysics Vol. 9, 91-119 (2007)] to compute expected heat flux [W/m^2] */
			heat_flux = 0.089*pow(markkt[mm2],(2.000/3.000))*pow(ABSV(markk[mm1]-markk[0]),(4.000/3.000))*pow((mro*mro*markbb[mm2]*mga*mcp/eta_fluid),(1.000/3.000));
			/**/
			/* Keep mkt at standard value [W/(m*K)], if heat flux is very small */
			if(heat_flux<=0.0001)
				{
				mkt = markkt[mm2];
				}
			/**/
			/* Calculate new thermal conductivity [W/(m*K)] taking numerically used viscosity into account */
			else if(heat_flux>0.0001)
				{
				mkt = pow((heat_flux/0.089),(3.000/2.000))/(pow((markk[mm1]-markk[0]),2.000)*mro)*pow((markbb[mm2]*mga*mcp/mnu),(-1.000/2.000));
				}
			/* End of critical melt fraction loop 1 */
			}
		/**/
		/* Fully molten silicates above 60 % of Si melt */
		if(xmelt>0.6000)
			{
			/* First calculate the correct heat flux for largely molten silicates using eta_fluid = 100 Pa s [Reese et al., JGR, 115, E05004 (2010)] */
                        eta_fluid = 100.0000;
			/**/
			/* Apply the soft turbulence model [Solomatov, Treatise on Geophysics Vol. 9, 91-119 (2007)] to compute expected heat flux [W/m^2] */
			heat_flux = 0.089*pow(markkt[mm2],(2.000/3.000))*pow(ABSV(markk[mm1]-markk[0]),(4.000/3.000))*pow((mro*mro*markbb[mm2]*mga*mcp/eta_fluid),(1.000/3.000));
			/**/
			/* Keep mkt at standard value [W/(m*K)], if heat flux is very small */
			if(heat_flux<=0.0001)
				{
				mkt = markkt[mm2];
				}
			/**/
			/* Calculate new thermal conductivity [W/(m*K)] taking numerically used viscosity into account */
			else if(heat_flux>0.0001)
				{
				mkt = pow((heat_flux/0.089),(3.000/2.000))/(pow((markk[mm1]-markk[0]),2.000)*mro)*pow((markbb[mm2]*mga*mcp/mnu),(-1.000/2.000));
				}
			/* End of critical melt fraction loop 2 */
			}
		/* End of molten silicates loop */
		}
	/**/
	/**/
	/* Only valid for iron, not used for very first timestep */
        if((mm2==7 || mm2==8 || mm2==9 || mm2==10 || mm2==17 || mm2==18 || mm2==19) && (timesum>core_form_time))
		{
		if(xmelt_fe==0.0000)  /* Iron is solid */
			{
			/**/
			/* First calculate the correct heat flux for solid iron using eta_fluid_fe = 1e12 Pa s [Yunker & Van Orman, EPSL, 254, 203-213 (2007)] */
                        eta_fluid_fe = 1.0000e12;
			/**/
			/* Apply the soft turbulence model [Solomatov, Treatise on Geophysics Vol. 9, 91-119 (2007)] to compute expected heat flux [W/m^2] */
			heat_flux = 0.089*pow(markkt[mm2],(2.000/3.000))*pow(ABSV(markk[mm1]-markk[0]),(4.000/3.000))*pow((mro*mro*markbb[mm2]*mga*mcp/eta_fluid_fe),(1.000/3.000));
			/**/
			/* Keep mkt at standard value [W/(m*K)], if heat flux is very small */
			if(heat_flux<=0.0001)
				{
				mkt = markkt[mm2];
				}
			/**/
			/* Calculate new thermal conductivity [W/(m*K)] taking numerically used viscosity into account */
			else if(heat_flux>0.0001)
				{
				mkt = pow((heat_flux/0.089),(3.0000/2.0000))/(pow((markk[mm1]-markk[0]),2.000)*mro)*pow((markbb[mm2]*mga*mcp/mnu),(-1.000/2.000));
				}
			}
		/**/
		if(xmelt_fe==1.0000) /* Iron is molten */
			{
			/**/
			/* First calculate the correct heat flux for molten iron using eta_fluid_fe = 1e-2 Pa s [Rubie et al., Treatise on Geophysics Vol. 9, 51-90 (2007)] */
                        eta_fluid_fe = 0.010;
			/**/
			/* Apply the soft turbulence model [Solomatov, Treatise on Geophysics Vol. 9, 91-119, (2007)] to compute expected heat flux [W/m^2] */
			heat_flux = 0.089*pow(markkt[mm2],(2.000/3.000))*pow(ABSV(markk[mm1]-markk[0]),(4.000/3.000))*pow((mro*mro*markbb[mm2]*mga*mcp/eta_fluid_fe),(1.000/3.000));
			/**/
			/* Keep mkt at standard value [W/(m*K)], if heat flux is very small */
			if(heat_flux<=0.0001)
				{
				mkt = markkt[mm2];
				}
			/**/
			/* Calculate new thermal conductivity [W/(m*K)] taking numerically used viscosity into account */
			else if(heat_flux>0.0001)
				{
				mkt = pow((heat_flux/0.089),(3.000/2.000))/(pow((markk[mm1]-markk[0]),2.000)*mro)*pow((markbb[mm2]*mga*mcp/mnu),(-1.000/2.000));
				}
			}
		}
	/**/
	/* Consider effect of material porosity on thermal conductivity of solid silicates in small planetesimals */
	/* by Greg, based on [Henke et al., Astronomy & Astrophysics, 537, A45 (2012)] (last modified: 30/06/2015) */
	if(mm2==5 || mm2==6)
		{
		/* Define non-dim. parameters */
		por_ref1 = 0.080;
		por_ref2 = 0.167;
		a_fact   = -1.200;
		/**/
		/* based on Henke et al. (2012) eq. (28) */
		if(markpor[mm1]<0.200)
			{
			mkt = markkt[mm2]*exp(-(double)(markpor[mm1])/por_ref1);
			}
		/**/
		/* based on Henke et al. (2012) eq. (29) */
		if(markpor[mm1]>=0.200 && markpor[mm1]<=0.400)
			{
			mkt = pow((pow((markkt[mm2]*exp(-(double)(markpor[mm1])/por_ref1)),4.000)+pow((markkt[mm2]*exp(a_fact-((double)(markpor[mm1])/por_ref2))),4.000)),0.250);
			}
		/**/
		/* based on Henke et al. (2012) eq. (30) */
		if(markpor[mm1]>0.400)
			{
			mkt = markkt[mm2]*exp(a_fact-((double)(markpor[mm1])/por_ref2));
			}
		}
	/**/
	/* Introduce cut-off values for thermal conductivity [W/(m*K)] (due to numerical reasons) */
	if(mkt<=0.0010)
		{
		mkt = 0.0010;  /* Otherwise zeros get into the global matrix */
		}
	if(mkt>=k_cutoff)
		{
		mkt = k_cutoff;
		}
	/**/
	/**/
	/* Check whether magnetic minerals in solid silicates are below their Curie temperatures [K] */
	/* by Greg (last modified: 18/10/2011) */
	/**/
	/* Define the marker magnetization parameter */
	mag = 3;  /* 0: no magnetization, 1: hematite magnetization, 2: magnetite and hematite magnetization, 3: no magnetization possible */
	/**/
	/* Dynamo is inactive */
	if(dynamo==0)
		{
		/* Solid silicates */
		if(mm2==5 || mm2==6)
			{
			/* BOTH magnetite and hematite were already previously magnetized [Lowrie, Fundamentals of Geophysics, p. 244 (2003)] */
			if(mtk<=851.16 && markmg_old[mm1]==2)
				{
				mag = 2;
				}
			else if((mtk>851.16 && mtk<=948.16) && markmg_old[mm1]==2)
				{
				mag = 1;
				}
			/**/
			/* ONLY hematite was already previously magnetized [Lowrie, Fundamentals of Geophysics, p. 240 (2003)] */
			else if(mtk<=851.16 && markmg_old[mm1]==1)
				{
				mag = 1;
				}
			else if((mtk>851.16 && mtk<=948.16) && markmg_old[mm1]==1)
				{
				mag = 1;
				}
			/**/
			/* Silicates demagnetize completely, thus magnetization time is also reset */
			else if(mtk>948.16)
				{
				mag = 0;
				mag_time=-1.000;
				}
			}
		/* Other material than solid silicates can never magnetize or demagnetize, magnetization time is reset */
		else if(mm2<5 || mm2>6)
			{
			mag = 3;
			mag_time=-1.000;
			}
		}
	/**/
	/* Dynamo is active */
	else if(dynamo==1)
		{
		/**/
		/* Solid silicates */
		if(mm2==5 || mm2==6)
			{
			/* T below Curie T of both magnetite and hematite [Lowrie, Fundamentals of Geophysics, p. 244 (2003)] */
			if(mtk<=851.16)
				{
				mag = 2;
				/**/
				/* If mineral magnetization is new, save magnetization time [Ma] */
				if(markmg_old[mm1]<=1 || markmg_old[mm1]==3)
					{
					mag_time = timesum/(3600.000*24.000*365.250*1.000e+6);
					}
				}
			/* T below Curie T of hematite [Lowrie, Fundamentals of Geophysics, p. 240 (2003)] */
			else if(mtk>851.16 && mtk<=948.16)
				{
				mag = 1;
				/**/
				/* If mineral magnetization is new, save magnetization time [Ma] */
				if(markmg_old[mm1]==0 || markmg_old[mm1]==3)
					{
					mag_time = timesum/(3600.000*24.000*365.250*1.000e+6);
					}
				}
			/* T above Curie T of both minerals, no magnetization, magnetization time is reset */
			else if(mtk>948.16)
				{
				mag = 0;
				mag_time = -1.000;
				}
			}
		/* Other material than solid silicates magnetize or demagnetize, magnetization time is reset */
		else if(mm2<5 || mm2>6)
			{
			mag = 3;
			mag_time = -1.000;
			}
		}
	/**/
	/* Save magnetization data for next timestep */
	markmg_old[mm1]  = mag;
	markmg_time[mm1] = mag_time;
	/**/
	/**/
	/**/
	/* No Radiogenic heating in sticky air material */
	if(mm2==0)
		{
        	mht = 0.0000e-10;
		}
	/**/
	/* Calculate time-dependent radiogenic heating, values taken from [Barr & Canup, Icarus, 198, 163-177 (2008)] */
	/* (last updated: 31/10/2016) by Tim */
	/* 26Al estimates updated with information from [Castillo-Rogez et al. (2009)] and [Moskovitz & Gaidos (2011)] */
	/**/
	timesumA = timesum/(365.250*24.000*3600.000);                  /* Time in years */
	/**/
	if(mm2==5 || mm2==6 || mm2==25 || mm2==26)  /* Internal heating by lithophile radiogenic isotopes in molten/solid silicate phase */
	  	{
		/* Decay constants [1/a] */
	  	decay_heat_al = 9.6673e-7;                            		/*  26Al, assuming t1/2 = 0.717 Myr  */
	  	decay_heat_k  = 4.9867e-10;                            		/*  40K   */
          	decay_heat_u  = 9.8458e-10;                            		/*  235U  */
          	decay_heat_uu = 1.5541e-10;                            		/*  238U  */
          	decay_heat_th = 4.9511e-11;                            		/*  232Th */
		/**/
		/* Initial heating rate [W/kg] */
          	init_al       = 1.89716e+23*al2627_init*4.999e-13/(3.264e+13);  /*  26Al, Moskovitz & Gaidos (2009), Eq. (5)  */
		/* Assume potassium depletion in target body: 300 ppm K in silicate material [Waenke & Dreibus, Phil. Trans. A, 349, 285-293 (1994)] */
          	init_k        = 0.40689*1.430e-11;                     		/*  40K   */
          	init_u        = 2.990e-12;                             		/*  235U  */
          	init_uu       = 1.600e-12;                             		/*  238U  */
          	init_th       = 1.000e-12;                             		/*  232Th */
		/**/
		/* Individual heating terms [W/kg] */
	  	mht_al        = init_al*exp(-timesumA*decay_heat_al);
	  	mht_k         = init_k*exp(-timesumA*decay_heat_k);
	  	mht_u         = init_u*exp(-timesumA*decay_heat_u);
	  	mht_uu        = init_uu*exp(-timesumA*decay_heat_uu);
	  	mht_th        = init_th*exp(-timesumA*decay_heat_th);
		/**/
		/* Total heating in silicate phase [W/m^3] */
	  	mht           = (mht_al+mht_k+mht_u+mht_uu+mht_th)*mro;
	  	}
	/**/
	if(mm2==7 || mm2==8 || mm2==9 || mm2==10 || mm2==17 || mm2==18 || mm2==19)  /* Internal heating by siderophile radiogenic isotopes in iron phase */
	  	{
		/* Decay constants [1/a] */
		decay_heat_fe = 2.6660e-7;              	                /*  60Fe, assuming t1/2 = 2.60 Myr from [Wallner et al., Phys. Rev. Lett., (2015)]  */
	  	decay_heat_k  = 4.9867e-10;   		                        /*  40K   */
		/**/
		/* Initial heating rate [W/kg] */
		init_fe       = 1.971e+24*fe6056_init*4.3451e-13/(1.183e+14);   /*  60Fe, Moskovitz & Gaidos (2009), Eq. (5)  */
		/* Assume potassium depletion in target body: 300 ppm K in iron material */
          	init_ka       = 0.090625*1.430e-11;                    		/*  40K   */
		/**/
		/* Individual heating terms [W/kg] */
	  	mht_fe        = init_fe*exp(-timesumA*decay_heat_fe);
	  	mht_ka        = init_ka*exp(-timesumA*decay_heat_k);
		/**/
		/* Total heating in iron phase [W/m^3] */
	  	mht           = (mht_fe+mht_ka)*mro;
	  	}
	/**/
        if(mm2==15)                     /* Internal heating in chondritic protocore material */
        	{
		/**/
		/* Decay constants [1/a] */
	  	decay_heat_al = 9.6673e-7;                                      /*  26Al, assuming t1/2 = 0.717 Myr  */
          	decay_heat_u  = 9.8458e-10;                            		/*  235U  */
          	decay_heat_uu = 1.5541e-10;                            		/*  238U  */
          	decay_heat_th = 4.9511e-11;                            		/*  232Th */
		decay_heat_fe = 2.6660e-7;                                      /*  60Fe, assuming t1/2 = 2.60 Myr from [Wallner et al., Phys. Rev. Lett., (2015)]  */
	  	decay_heat_k  = 4.9867e-10;                            		/*  40K   */
		/**/
		/* Initial heating rate [W/kg] */
          	init_al       = 1.89716e+23*al2627_init*4.999e-13/(3.264e+13);  /*  26Al, Moskovitz & Gaidos (2009), Eq. (5)  */ 
          	init_u        = 2.990e-12;                             		/*  235U  */
          	init_uu       = 1.600e-12;                            	 	/*  238U  */
          	init_th       = 1.000e-12;                             		/*  232Th */
          	init_fe       = 1.971e+24*fe6056_init*4.3451e-13/(1.183e+14);   /*  60Fe, Moskovitz & Gaidos (2009), Eq. (5)  */
		/* Assume NO potassium depletion in chondritic material */
          	init_k        = 1.430e-11;                             		/*  40K   */
		/**/
		/* Individual heating terms [W/kg] */
	  	mht_al        = init_al*exp(-timesumA*decay_heat_al);
	  	mht_u         = init_u*exp(-timesumA*decay_heat_u);
	  	mht_uu        = init_uu*exp(-timesumA*decay_heat_uu);
	  	mht_th        = init_th*exp(-timesumA*decay_heat_th);
	  	mht_fe        = init_fe*exp(-timesumA*decay_heat_fe);
	  	mht_k         = init_k*exp(-timesumA*decay_heat_k);
		/**/
		/* Total heating in chondritic phase [W/m^3] */
	  	mht           = (mht_al+mht_fe+mht_k+mht_u+mht_uu+mht_th)*mro;
	  	}
	/**/
	/**/
	/* Interpolation from markers to nodes ====================================*/
	/* Marker weight calculation using dimension of current Cell */
	celdx=gx[wn[0]+1]-gx[wn[0]];
	celdy=gy[wn[1]+1]-gy[wn[1]];
	/**/
	swt1=1.0;
	/**/
	/* Horizontal, vertical limits for interpolation calc */
	m1min=wn[0]-intermod1; if(m1min<0) m1min=0;
	m1max=wn[0]+1+intermod1; if(m1max>xnumx-1) m1max=xnumx-1;
	/**/
	m2min=wn[1]-intermod1; if(m2min<0) m2min=0;
	m2max=wn[1]+1+intermod1; if(m2max>ynumy-1) m2max=ynumy-1;
	/**/
	/* Interpolation weights calc after Fornberg (1996) */
	nodewt(m1min,m1max,m2min,m2max,(double)(markx[mm1]),(double)(marky[mm1]),0,0);
	/**/
	/* Interpolate ro,Nu etc to nodes using interpolation coefficients */
	/* Add Normal viscosity nd */
	ddx=1.0-ABSV(cn[wn[0]-m1min][1]-0.5);
	ddy=1.0-ABSV(cn[wn[1]-m2min][0]-0.5);
	swt2=swt1*ddx*ddy;
	m3=(wn[0]+1)*ynumy+wn[1]+1;
	/**/
	/* New interpolation scheme */
        if(ddx>0.5 && ddy>0.5)
                {
                swt2=swt1*4.0*(ddx-0.5)*(ddy-0.5);
                if(viscmod==0) nd0[m3]+=mnu*swt2;
                if(viscmod==1) nd0[m3]+=log(mnu)*swt2;
                if(viscmod==2) nd0[m3]+=1.0/mnu*swt2;
                sol1[nodenum+m3]+=swt2;
                }
	/**/
	/* Interpolate ro,Nu etc to nodes using interpolation coefficients */
	/* Reload horizontal coefficients to cn[] */
	for (m1=m1min;m1<=m1max;m1++)
	for (m2=m2min;m2<=m2max;m2++)
		{
		/* Current node num, wt */
		m3=m1*ynumy+m2;
		swt=swt1*cn[m1-m1min][1]*cn[m2-m2min][0];
		ddx=cn[m1-m1min][1];
		ddy=cn[m2-m2min][0];
		swt2=swt1*ddx*ddy;
		/**/
		/* Add physical properties: ro,nu, etc. */
                /* Add Physical Properties: ro,nu, etc. */
                if(ddx>0.5 && ddy>0.5)
                        {
                        swt2=swt1*4.0*(ddx-0.5)*(ddy-0.5);
                        if(viscmod==0) nu0[m3]+=mnu*swt2;
                        if(viscmod==1) nu0[m3]+=log(mnu)*swt2;
                        if(viscmod==2) nu0[m3]+=1.0/mnu*swt2;
                        sol0[nodenum+m3]+=swt2;
                        }
		/* Add other physical properties: ro, cp etc. */
		ro0[m3]+=mro*swt;
		et0[m3]+=mbb*swt;
		cp0[m3]+=mcp*swt;
		kt0[m3]+=mkt*swt;
		ht0[m3]+=mht*swt;
		sol0[m3]+=swt;
		/**/
		/* Add T */
		if(!markim[mm2])
			{
			tk0[m3]+=mtk*swt;
			sol1[m3]+=swt;
			}
		}
	/* End Interpolation from markers to nodes ====================================*/
	}
}
/**/
/**/
/* Recalc ro[] nu[] */
for (m3=0;m3<nodenum;m3++)
{
/* Shear viscosity */
if(sol0[nodenum+m3])
	{
	/* Viscosity recalc check */
	if(mu[m3])
		{
		nu0[m3]=mu[m3];
		}
	else
		{
		nu0[m3]/=sol0[nodenum+m3];
		if(viscmod==1) nu0[m3]=exp(nu0[m3]);
		if(viscmod==2) nu0[m3]=1.0/nu0[m3];
		}
	/* Min,Max NU limitation */
	/**/
	if(nu0[m3]<nubeg) nu0[m3]=nubeg;
	if(nu0[m3]>nuend) nu0[m3]=nuend;
	/* Min,Max NU definition for nu contrast limit */
	if(numin==0 || nu0[m3]<numin) numin=nu0[m3];
	if(numax==0 || nu0[m3]>numax) numax=nu0[m3];
	nu[m3]=nu0[m3];
	/**/
	/* Reset weight */
	sol0[nodenum+m3]=0;
	}
/* Normal viscosity */
if(sol1[nodenum+m3])
	{
	/* Viscosity recalc check */
	if(mu[m3])
		{
		nd0[m3]=mu[m3];
		}
	else
		{
		nd0[m3]/=sol1[nodenum+m3];
		if(viscmod==1) nd0[m3]=exp(nd0[m3]);
		if(viscmod==2) nd0[m3]=1.0/nd0[m3];
		}
	/* Min,Max NU limitation */
	/**/
	if(nd0[m3]<nubeg) nd0[m3]=nubeg;
	if(nd0[m3]>nuend) nd0[m3]=nuend;
	/* Min,Max NU definition for nu contrast limit */
	if(numin==0 || nd0[m3]<numin) numin=nd0[m3];
	if(numax==0 || nd0[m3]>numax) numax=nd0[m3];
	nd[m3]=nd0[m3];
	/* Old pressure from markers */
	po[m3]=po0[m3]/sol1[nodenum+m3];
	/**/
	/* Reset weight */
	sol1[nodenum+m3]=0;
	}
/* Other material properties */
if(sol0[m3])
	{
	/* Material constants recalc */
	ro[m3]=ro0[m3]/sol0[m3];
	et[m3]=et0[m3]/sol0[m3];
	cp[m3]=cp0[m3]/sol0[m3];
	kt[m3]=kt0[m3]/sol0[m3];
	ht[m3]=ht0[m3]/sol0[m3];
	/**/
	/* Advective addition for T K in nodes recalc */
	if (sol1[m3])
		{
		tk[m3]=tk0[m3]/sol1[m3];
		/**/
		sol1[m3]=0;
		}
	/**/
	/* Reset weight */
	sol0[m3]=0;
	}
}
printf("Min, Max viscosity %e %e \n",numin,numax);
/**/
/**/
/* Reset advective temperature */
for (m3=0;m3<nodenum;m3++) tk3[m3]=0;
/**/
/* Set Upper/Lower limits for nu[] after given contrast */
if(nucontr>=1.0) numax=numin*nucontr; else numin=numax*nucontr;
for (m3=0;m3<nodenum;m3++)
	{
	if(nu[m3]<numin) nu[m3]=numin;
	if(nu[m3]>numax) nu[m3]=numax;
	if(nd[m3]<numin) nd[m3]=numin;
	if(nd[m3]>numax) nd[m3]=numax;
	}
/**/
/* Set Boundary conditions for T */
if (printmod) printf("\n AVERAGE TEMPERATURE CORRECTION FOR BOUNDARY CONDITIONS ...\n");
tkrecalc();
if (printmod) printf("AVERAGE TEMPERATURE OK!\n");
/**/
/**/
/**/
/* ADD MARKERS TO THE v-CELLS ========================== */
}
/* End ro[],nu[] recalc after marker positions */
/**/
/**/
/* Calc ro for given P,T */
double dencalc(double mtk, double mpb, double markporos, double x, double y, int mm2)
/* mtk     - T [K] */
/* mpb     - P [bar] */
/* markpor - porosity of solid silicates [non-dim.] */
/* x,y     - XY location of point for Vx,Vy calc */
/* mm2     - Rock number */
{
/* Val buffer */
double ival;
eps[20]=0;
/**/
/* Constant density */
if (densimod==0) return markro[mm2];
/**/
/* Ro=ro0*(1-bro*(TK-298.15))*(1+aro*(Pkbar-0.001)) */
/* Adiabatic term: al=bro/(1-bro*(Tk-298.15)) */
/**/
/* Added material porosity term (by Greg, last modified: 10/02/2014) */
/* addition based on Henke et al. (2012) eq. (3) */
ival=markro[mm2]*(1.0-markbb[mm2]*(mtk-298.15))*(1.0+markaa[mm2]*(mpb-1.0)*1e-3)*(1.000-markporos);
eps[20]=markbb[mm2]/(1.0-markbb[mm2]*(mtk-298.15));
return  ival;
}
/* End Calc ro for given P,T */
/**/
/**/
/* Number of nearest left vertical line find */
long int m1serch(double x)
/* x - X coordinate */
{
/* Variables */
long int m1,m10=0,m11=xnumx-1;
/**/
/* Search cycle */
do
	{
	m1=(m10+m11)/2;
	if (gx[m1]>x) m11=m1; else m10=m1;
	}
while((m11-m10)>1);
if(m10>xnumx-2) m10=xnumx-2;
/**/
return m10;
}
/* Number of nearest left vertical line found */
/**/
/**/
/* Number of nearest upper horizontal line found */
long int m2serch(double y)
/* y - Y coordinate */
{
/* Variables */
long int m2,m20=0,m21=ynumy-1;
/**/
/* Search cycle */
do
	{
	m2=(m20+m21)/2;
	if (gy[m2]>y) m21=m2; else m20=m2;
	}
while((m21-m20)>1);
if(m20>ynumy-2) m20=ynumy-2;
/**/
return m20;
}
/* Number of nearest upper horizontal line found */
/**/
/**/
/* Calculation of Vx,Vy,EPSxx,EPSyy,EPSxy,SIGxx,SIGyy,SIGxy,T,T0,T1,T2,P for current location by linear interpolation */
/* Staggered Nodes num */
/*   [0]                [3]                [6]   */
/*  T0,xy0    Vy0     T3,xy3     Vy3             */
/*                                               */
/*   Vx0    P4,xx4,yy4  Vx3    P7,xx7,yy7        */
/*                                               */
/*   [1]                [4]                [7]   */
/*  T,xy1     Vy1     T4,xy4     Vy4             */
/*                                               */
/*   Vx1    P5,xx5,yy5  Vx4    P8,xx8,yy8        */
/*                                               */
/*   [2]                [5]                [8]   */
/*                                               */
/*                                               */
void allinter(double x, double y)
/* x,y - XY location of point for Vx,Vy calc */
{
/* Counters */
long int m1,m2,m3,m10,m20,m1min,m1max,m2min,m2max;
/* en-NormalisedDistance */
double ival;
/**/
/**/
/* Check X,Y */
if(x<0) x=0; else if(x>xsize) x=xsize;
if(y<0) y=0; else if(y>ysize) y=ysize;
/**/
/**/
/**/
/* Upper Left Node X,Y Num */
wn[0]=m10=m1serch(x);
wn[1]=m20=m2serch(y);
/**/
/**/
/**/
/* Buffer clear */
for (m1=0;m1<=15;m1++) eps[m1]=0;
/**/
/**/
/**/
/* T interpolation ------------------------ */
/* Horizontal,Vertical limits for interpolation calc */
m1min=m10-intermod; if(m1min<0) m1min=0;
m1max=m10+1+intermod; if(m1max>xnumx-1) m1max=xnumx-1;
/**/
m2min=m20-intermod; if(m2min<0) m2min=0;
m2max=m20+1+intermod; if(m2max>ynumy-1) m2max=ynumy-1;
/**/
/* Interpolation weights calc after Fornberg (1996) */
nodewt(m1min,m1max,m2min,m2max,x,y,0,0);
/**/
/* T interpolate after interpolation weights */
for (m1=m1min;m1<=m1max;m1++)
for (m2=m2min;m2<=m2max;m2++)
	{
	/* Current node num, wt */
	m3=m1*ynumy+m2;
	ival=cn[m1-m1min][1]*cn[m2-m2min][0];
	eps[0]+=ival*tk0[m3];
	eps[1]+=ival*tk1[m3];
	eps[2]+=ival*tk[m3];
	eps[3]+=ival*tk2[m3];
	}
/**/
/* Wt for nodes save for T */
wn[2]=m1min; wn[3]=m1max;
for (m1=m1min;m1<=m1max;m1++)
	{
	cn[m1-m1min][5]=cn[m1-m1min][1];
	}
wn[4]=m2min; wn[5]=m2max;
for (m2=m2min;m2<=m2max;m2++)
	{
	cn[m2-m2min][4]=cn[m2-m2min][0];
	}
/* End T interpolation ------------------------ */
/**/
/**/
/**/
/* SIGxy,EPSxy, SIGxy*EPSxy interpolation ------------------------ */
/* Horizontal,Vertical limits for interpolation calc */
m1min=m10; if(m1min<1) m1min=1; if(m1min>xnumx-3) m1min=xnumx-3;
m1max=m1min+1+intermod; if(m1max>xnumx-2) m1max=xnumx-2;
m1min=m1min-intermod; if(m1min<1) m1min=1;
/**/
m2min=m20; if(m2min<1) m2min=1; if(m2min>ynumy-3) m2min=ynumy-3;
m2max=m2min+1+intermod; if(m2max>ynumy-2) m2max=ynumy-2;
m2min=m2min-intermod; if(m2min<1) m2min=1;
/**/
/* Interpolation weights calc after Fornberg (1996) */
nodewt(m1min,m1max,m2min,m2max,x,y,0,0);
/**/
/* SIGxy,EPSxy Interpolate after interpolation weights */
for (m1=m1min;m1<=m1max;m1++)
for (m2=m2min;m2<=m2max;m2++)
	{
	/* Current node num, wt */
	m3=m1*ynumy+m2;
	ival=cn[m1-m1min][1]*cn[m2-m2min][0];
	eps[4]+=ival*exy[m3];
	eps[5]+=ival*sxy[m3];
	eps[13]+=ival*exy[m3]*sxy[m3];
	}
/* End SIGxy,EPSxy interpolation ------------------------ */
/**/
/**/
/**/
/* SIGxx,EPSxx,SIGyy,EPSyy,P, SIGxx*EPSxx, SIGyy*EPSyy interpolation ------------------------ */
/* Horizontal, vertical limits for interpolation calc */
m1min=m10; if(x>(gx[m10]+gx[m10+1])/2.0) m1min+=1;
if(m1min<1) m1min=1; if(m1min>xnumx-2) m1min=xnumx-2;
m1max=m1min+1+intermod; if(m1max>xnumx-1) m1max=xnumx-1;
m1min=m1min-intermod; if(m1min<1) m1min=1;
/**/
m2min=m20; if(y>(gy[m20]+gy[m20+1])/2.0) m2min+=1;
if(m2min<1) m2min=1; if(m2min>ynumy-2) m2min=ynumy-2;
m2max=m2min+1+intermod; if(m2max>ynumy-1) m2max=ynumy-1;
m2min=m2min-intermod; if(m2min<1) m2min=1;
/**/
/* Interpolation weights calc after Fornberg (1996) */
nodewt(m1min,m1max,m2min,m2max,x,y,-1,-1);
/**/
/* SIGxx,EPSxx,SIGyy,EPSyy,P Interpolate after interpolation weights */
for (m1=m1min;m1<=m1max;m1++)
for (m2=m2min;m2<=m2max;m2++)
	{
	/* Current node num, wt */
	m3=m1*ynumy+m2;
	ival=cn[m1-m1min][1]*cn[m2-m2min][0];
	eps[6]+=ival*exx[m3];
	eps[7]+=ival*sxx[m3];
	eps[8]+=ival*eyy[m3];
	eps[9]+=ival*syy[m3];
	eps[10]+=ival*pr[m3];
	eps[14]+=ival*exx[m3]*sxx[m3];
	eps[15]+=ival*eyy[m3]*syy[m3];
	}
/**/
/* End SIGxx,EPSxx,SIGyy,EPSyy,P interpolation ------------------------ */
/**/
/**/
/**/
/* Vx interpolation ------------------------ */
/* Horizontal,Vertical limits for interpolation calc */
m1min=m10-intermod; if(m1min<0) m1min=0;
m1max=m10+1+intermod; if(m1max>xnumx-1) m1max=xnumx-1;
/**/
m2min=m20; if(y<(gy[m20]+gy[m20+1])/2.0) m2min-=1;
if(m2min<0) m2min=0; if(m2min>ynumy-3) m2min=ynumy-3;
m2max=m2min+1+intermod; if(m2max>ynumy-2) m2max=ynumy-2;
m2min=m2min-intermod; if(m2min<0) m2min=0;
/**/
/* Interpolation weights calc after Fornberg (1996) */
nodewt(m1min,m1max,m2min,m2max,x,y,0,+1);
/**/
/* Vx Interpolate after interpolation weights */
for (m1=m1min;m1<=m1max;m1++)
for (m2=m2min;m2<=m2max;m2++)
	{
	/* Current node num, wt */
	m3=m1*ynumy+m2;
	ival=cn[m1-m1min][1]*cn[m2-m2min][0];
	eps[11]+=ival*vx[m3];
	}
/* End Vx interpolation ------------------------ */
/**/
/**/
/**/
/* Vy interpolation ------------------------ */
/* Horizontal,Vertical limits for interpolation calc */
m1min=m10; if(x<(gx[m10]+gx[m10+1])/2.0) m1min-=1;
if(m1min<0) m1min=0; if(m1min>xnumx-3) m1min=xnumx-3;
m1max=m1min+1+intermod; if(m1max>xnumx-2) m1max=xnumx-2;
m1min=m1min-intermod; if(m1min<0) m1min=0;
/**/
m2min=m20-intermod; if(m2min<0) m2min=0;
m2max=m20+1+intermod; if(m2max>ynumy-1) m2max=ynumy-1;
/**/
/* Interpolation weights calc after Fornberg (1996) */
nodewt(m1min,m1max,m2min,m2max,x,y,+1,0);
/**/
/* Vx Interpolate after interpolation weights */
for (m1=m1min;m1<=m1max;m1++)
for (m2=m2min;m2<=m2max;m2++)
	{
	/* Current node num, wt */
	m3=m1*ynumy+m2;
	ival=cn[m1-m1min][1]*cn[m2-m2min][0];
	eps[12]+=ival*vy[m3];
	}
/* End Vy interpolation ------------------------ */
/**/
/**/
/**/
}
/* Calculation of Vx,Vy,EPSxx,EPSyy,EPSxy,SIGxx,SIGyy,SIGxy,T,T0,T1,T2,P for current location by linear interpolation */
/**/
/**/
/* Weights for horisontal and vertical nodes calculation for marker interpolation */
void nodewt(long int m1min, long int m1max, long int m2min, long int m2max, double x, double y, int ynx, int yny)
/* m1min,m1max, m2min,m2max - node X,Y number limits */
/* x,y - coordinates of current point */
/* ynx, yny - Type of shifts: No(0), Back(-1), Forw(1) */
{
/* Eyy vertical position */
long int m3;
int nx,ny;
/**/
/**/
/**/
/* Weigths in horizontal directions */
/* Load distances to xn[] */
if(ynx<0)
	{
	for (m3=m1min;m3<=m1max;m3++)
		{
		xn[m3-m1min]=(gx[m3]+gx[m3-1])/2.0;
		}
	}
if(ynx==0)
	{
	for (m3=m1min;m3<=m1max;m3++)
		{
		xn[m3-m1min]=gx[m3];
		}
	}
if(ynx>0)
	{
	for (m3=m1min;m3<=m1max;m3++)
		{
		xn[m3-m1min]=(gx[m3]+gx[m3+1])/2.0;
		}
	}
/**/
/* Calc maximum position in xn[] */
nx=(int)(m1max-m1min);
/**/
/* Calc coefficients for horizontal direction */
fdweight(nx,0,x);
/**/
/* Reload horizontal coefficients to cn[] */
for (m3=0;m3<=nx;m3++)
	{
	cn[m3][1]=cn[m3][0];
	}
/**/
/**/
/**/
/* Weigths in vertical directions */
/* Load distances to xn[] */
if(yny<0)
	{
	for (m3=m2min;m3<=m2max;m3++)
		{
		xn[m3-m2min]=(gy[m3]+gy[m3-1])/2.0;
		}
	}
if(yny==0)
	{
	for (m3=m2min;m3<=m2max;m3++)
		{
		xn[m3-m2min]=gy[m3];
		}
	}
if(yny>0)
	{
	for (m3=m2min;m3<=m2max;m3++)
		{
		xn[m3-m2min]=(gy[m3]+gy[m3+1])/2.0;
		}
	}
/**/
/* Calc maximum position in xn[] */
ny=(int)(m2max-m2min);
/**/
/* Calc coefficients for horizontal direction */
fdweight(ny,0,y);
/**/
/**/
/**/
}
/* Weights for horizontal and vertical nodes calculation for marker interpolation */
/**/
/**/
/* Calculation of Vx,Vy by interpolation */
void allinterv(double x, double y)
/* x,y - XY location of point for Vx,Vy calc */
{
/* Counters */
long int m1,m2,m3,m10,m20,m1min,m1max,m2min,m2max;
/* en - normalised distance */
double ival;
/**/
/**/
/* Check X,Y */
if(x<0) x=0; else if(x>xsize) x=xsize;
if(y<0) y=0; else if(y>ysize) y=ysize;
/**/
/**/
/**/
/* Up Left Node X,Y Num */
wn[0]=m10=m1serch(x);
wn[1]=m20=m2serch(y);
/**/
/**/
/**/
/* Vx interpolation ------------------------ */
/* Buffer clear */
eps[11]=0;
/* Horizontal, vertical limits for interpolation calc */
m1min=m10-intermod; if(m1min<0) m1min=0;
m1max=m10+1+intermod; if(m1max>xnumx-1) m1max=xnumx-1;
/**/
m2min=m20; if(y<(gy[m20]+gy[m20+1])/2.0) m2min-=1;
if(m2min<0) m2min=0; if(m2min>ynumy-3) m2min=ynumy-3;
m2max=m2min+1+intermod; if(m2max>ynumy-2) m2max=ynumy-2;
m2min=m2min-intermod; if(m2min<0) m2min=0;
/**/
/* Interpolation weights calc after Fornberg (1996) */
nodewt(m1min,m1max,m2min,m2max,x,y,0,+1);
/**/
/* Vx interpolate after interpolation weights */
for (m1=m1min;m1<=m1max;m1++)
for (m2=m2min;m2<=m2max;m2++)
	{
	/* Current node num, wt */
	m3=m1*ynumy+m2;
	ival=cn[m1-m1min][1]*cn[m2-m2min][0];
	eps[11]+=ival*vx[m3];
	}
/* End Vx interpolation ------------------------ */
/**/
/**/
/**/
/* Vy interpolation ------------------------ */
/* Buffer clear */
eps[12]=0;
/* Horizontal, vertical limits for interpolation calc */
m1min=m10; if(x<(gx[m10]+gx[m10+1])/2.0) m1min-=1;
if(m1min<0) m1min=0; if(m1min>xnumx-3) m1min=xnumx-3;
m1max=m1min+1+intermod; if(m1max>xnumx-2) m1max=xnumx-2;
m1min=m1min-intermod; if(m1min<0) m1min=0;
/**/
m2min=m20-intermod; if(m2min<0) m2min=0;
m2max=m20+1+intermod; if(m2max>ynumy-1) m2max=ynumy-1;
/**/
/* Interpolation weights calc after Fornberg (1996) */
nodewt(m1min,m1max,m2min,m2max,x,y,+1,0);
/**/
/* Vy Interpolate after interpolation weights */
for (m1=m1min;m1<=m1max;m1++)
for (m2=m2min;m2<=m2max;m2++)
	{
	/* Current node num, wt */
	m3=m1*ynumy+m2;
	ival=cn[m1-m1min][1]*cn[m2-m2min][0];
	eps[12]+=ival*vy[m3];
	}
/* End Vy interpolation ------------------------ */
/**/
}
/* Calculation of Vx,Vy by interpolation */
/**/
/**/
/* Calculation of EPSxx,EPSyy,EPSxy, P by interpolation */
void allintere(double x, double y)
/* x,y - XY location of point for Vx,Vy calc */
{
/* Counters */
long int m1,m2,m3,m10,m20,m1min,m1max,m2min,m2max;
/* en - normalised distance */
double ival;
/**/
/**/
/* Check X,Y */
if(x<0) x=0; else if(x>xsize) x=xsize;
if(y<0) y=0; else if(y>ysize) y=ysize;
/**/
/**/
/**/
/* Up Left Node X,Y Num */
wn[0]=m10=m1serch(x);
wn[1]=m20=m2serch(y);
/**/
/**/
/**/
/* EPSxy interpolation ------------------------ */
/* Buffer clear */
eps[4]=0;
eps[5]=0;
/* Horizontal, vertical limits for interpolation calc */
m1min=m10; if(m1min<1) m1min=1; if(m1min>xnumx-3) m1min=xnumx-3;
m1max=m1min+1+intermod; if(m1max>xnumx-2) m1max=xnumx-2;
m1min=m1min-intermod; if(m1min<1) m1min=1;
/**/
m2min=m20; if(m2min<1) m2min=1; if(m2min>ynumy-3) m2min=ynumy-3;
m2max=m2min+1+intermod; if(m2max>ynumy-2) m2max=ynumy-2;
m2min=m2min-intermod; if(m2min<1) m2min=1;
/**/
/* Interpolation weights calc after Fornberg (1996) */
nodewt(m1min,m1max,m2min,m2max,x,y,0,0);
/**/
/* EPSxy interpolate after interpolation weights */
for (m1=m1min;m1<=m1max;m1++)
for (m2=m2min;m2<=m2max;m2++)
	{
	/* Current node num, wt */
	m3=m1*ynumy+m2;
	ival=cn[m1-m1min][1]*cn[m2-m2min][0];
	eps[4]+=ival*exy[m3];
	eps[5]+=ival*sxy[m3];
	}
/* End SIGxy,EPSxy interpolation ------------------------ */
/**/
/**/
/**/
/* EPSxx,EPSyy, P interpolation ------------------------ */
/* Buffer clear */
eps[6]=eps[8]=eps[10]=0;
eps[7]=eps[9]=eps[49]=0;
/* Horizontal, vertical limits for interpolation calc */
m1min=m10; if(x>(gx[m10]+gx[m10+1])/2.0) m1min+=1;
if(m1min<1) m1min=1; if(m1min>xnumx-2) m1min=xnumx-2;
m1max=m1min+1+intermod; if(m1max>xnumx-1) m1max=xnumx-1;
m1min=m1min-intermod; if(m1min<1) m1min=1;
/**/
m2min=m20; if(y>(gy[m20]+gy[m20+1])/2.0) m2min+=1;
if(m2min<1) m2min=1; if(m2min>ynumy-2) m2min=ynumy-2;
m2max=m2min+1+intermod; if(m2max>ynumy-1) m2max=ynumy-1;
m2min=m2min-intermod; if(m2min<1) m2min=1;
/**/
/* Interpolation weights calc after Fornberg (1996) */
nodewt(m1min,m1max,m2min,m2max,x,y,-1,-1);
/**/
/* EPSxx,EPSyy,P Interpolate after interpolation weights */
for (m1=m1min;m1<=m1max;m1++)
for (m2=m2min;m2<=m2max;m2++)
	{
	/* Current node num, wt */
	m3=m1*ynumy+m2;
	ival=cn[m1-m1min][1]*cn[m2-m2min][0];
	eps[6]+=ival*exx[m3];
	eps[7]+=ival*sxx[m3];
	eps[8]+=ival*eyy[m3];
	eps[9]+=ival*syy[m3];
	eps[10]+=ival*pr[m3];
	eps[49]+=ival*ga[m3];
	}
/* End EPSxx,EPSyy,P interpolation ------------------------ */
/**/
}
/* Calculation of EPSxx,EPSyy,EPSxy, P by Interpolation */
/**/
/**/
/* Calculation of Vx,Vy, EPSxx*SIGxx,EPSyy*SIGyy,EPSxy*SIGxy by Interpolation */
void allinters(double x, double y)
/* x,y - XY location of point for Vx,Vy calc */
{
/* Counters */
long int m1,m2,m3,m10,m20,m1min,m1max,m2min,m2max;
/* en - normalised distance */
double ival;
/**/
/**/
/* Check X,Y */
if(x<0) x=0; else if(x>xsize) x=xsize;
if(y<0) y=0; else if(y>ysize) y=ysize;
/**/
/**/
/**/
/* Up Left Node X,Y Num */
wn[0]=m10=m1serch(x);
wn[1]=m20=m2serch(y);
/**/
/**/
/**/
/* SIGxy*EPSxy interpolation ------------------------ */
/* Buffer clear */
eps[13]=0;
/* Horizontal, vertical limits for interpolation calc */
m1min=m10; if(m1min<1) m1min=1; if(m1min>xnumx-3) m1min=xnumx-3;
m1max=m1min+1+intermod; if(m1max>xnumx-2) m1max=xnumx-2;
m1min=m1min-intermod; if(m1min<1) m1min=1;
/**/
m2min=m20; if(m2min<1) m2min=1; if(m2min>ynumy-3) m2min=ynumy-3;
m2max=m2min+1+intermod; if(m2max>ynumy-2) m2max=ynumy-2;
m2min=m2min-intermod; if(m2min<1) m2min=1;
/**/
/* Interpolation weights calc after Fornberg (1996) */
nodewt(m1min,m1max,m2min,m2max,x,y,0,0);
/**/
/* SIGxy,EPSxy Interpolate after interpolation weights */
for (m1=m1min;m1<=m1max;m1++)
for (m2=m2min;m2<=m2max;m2++)
	{
	/* Current node num, wt */
	m3=m1*ynumy+m2;
	ival=cn[m1-m1min][1]*cn[m2-m2min][0];
	eps[13]+=ival*exy[m3]*sxy[m3];
	}
/* End SIGxy*EPSxy interpolation ------------------------ */
/**/
/**/
/**/
/* SIGxx*EPSxx, SIGyy*EPSyy interpolation ------------------------ */
/* Buffer clear */
eps[14]=eps[15]=0;
/* Horizontal, vertical limits for interpolation calc */
m1min=m10; if(x>(gx[m10]+gx[m10+1])/2.0) m1min+=1;
if(m1min<1) m1min=1; if(m1min>xnumx-2) m1min=xnumx-2;
m1max=m1min+1+intermod; if(m1max>xnumx-1) m1max=xnumx-1;
m1min=m1min-intermod; if(m1min<1) m1min=1;
/**/
m2min=m20; if(y>(gy[m20]+gy[m20+1])/2.0) m2min+=1;
if(m2min<1) m2min=1; if(m2min>ynumy-2) m2min=ynumy-2;
m2max=m2min+1+intermod; if(m2max>ynumy-1) m2max=ynumy-1;
m2min=m2min-intermod; if(m2min<1) m2min=1;
/**/
/* Interpolation weights calc after Fornberg (1996) */
nodewt(m1min,m1max,m2min,m2max,x,y,-1,-1);
/**/
/* SIGxx,EPSxx,SIGyy,EPSyy,P interpolate after interpolation weights */
for (m1=m1min;m1<=m1max;m1++)
for (m2=m2min;m2<=m2max;m2++)
	{
	/* Current node num, wt */
	m3=m1*ynumy+m2;
	ival=cn[m1-m1min][1]*cn[m2-m2min][0];
	eps[14]+=ival*exx[m3]*sxx[m3];
	eps[15]+=ival*eyy[m3]*syy[m3];
	}
/* End SIGxx*EPSxx,SIGyy*EPSyy interpolation ------------------------ */
/**/
/**/
/**/
/* Vx interpolation ------------------------ */
/* Buffer clear */
eps[11]=0;
/* Horizontal, vertical limits for interpolation calc */
m1min=m10-intermod; if(m1min<0) m1min=0;
m1max=m10+1+intermod; if(m1max>xnumx-1) m1max=xnumx-1;
/**/
m2min=m20; if(y<(gy[m20]+gy[m20+1])/2.0) m2min-=1;
if(m2min<0) m2min=0; if(m2min>ynumy-3) m2min=ynumy-3;
m2max=m2min+1+intermod; if(m2max>ynumy-2) m2max=ynumy-2;
m2min=m2min-intermod; if(m2min<0) m2min=0;
/**/
/* Interpolation weights calc after Fornberg (1996) */
nodewt(m1min,m1max,m2min,m2max,x,y,0,+1);
/**/
/* Vx interpolate after interpolation weights */
for (m1=m1min;m1<=m1max;m1++)
for (m2=m2min;m2<=m2max;m2++)
	{
	/* Current node num, wt */
	m3=m1*ynumy+m2;
	ival=cn[m1-m1min][1]*cn[m2-m2min][0];
	eps[11]+=ival*vx[m3];
	}
/* End Vx interpolation ------------------------ */
/**/
/**/
/**/
/* Vy interpolation ------------------------ */
/* Buffer clear */
eps[12]=0;
/* Horizontal,Vertical limits for interpolation calc */
m1min=m10; if(x<(gx[m10]+gx[m10+1])/2.0) m1min-=1;
if(m1min<0) m1min=0; if(m1min>xnumx-3) m1min=xnumx-3;
m1max=m1min+1+intermod; if(m1max>xnumx-2) m1max=xnumx-2;
m1min=m1min-intermod; if(m1min<0) m1min=0;
/**/
m2min=m20-intermod; if(m2min<0) m2min=0;
m2max=m20+1+intermod; if(m2max>ynumy-1) m2max=ynumy-1;
/**/
/* Interpolation weights calc after Fornberg (1996) */
nodewt(m1min,m1max,m2min,m2max,x,y,+1,0);
/**/
/* Vy interpolate after interpolation weights */
for (m1=m1min;m1<=m1max;m1++)
for (m2=m2min;m2<=m2max;m2++)
	{
	/* Current node num, wt */
	m3=m1*ynumy+m2;
	ival=cn[m1-m1min][1]*cn[m2-m2min][0];
	eps[12]+=ival*vy[m3];
	}
/* End Vy interpolation ------------------------ */
/**/
}
/* Calculation of Vx,Vy, EPSxx*SIGxx,EPSyy*SIGyy,EPSxy*SIGxy by interpolation */
/**/
/**/
/* Calculation of T,T0 for current location by interpolation */
void allintert(double x, double y)
/* x,y - XY location of point for Vx,Vy calc */
{
/* Counters */
long int m1,m2,m3,m10,m20,m1min,m1max,m2min,m2max;
/* en - normalized distance */
double ival;
/**/
/**/
/* Check X,Y */
if(x<0) x=0; else if(x>xsize) x=xsize;
if(y<0) y=0; else if(y>ysize) y=ysize;
/**/
/**/
/**/
/* Up Left Node X,Y Num */
wn[0]=m10=m1serch(x);
wn[1]=m20=m2serch(y);
/**/
/**/
/**/
/* T interpolation ------------------------ */
/* Buffer clear */
eps[2]=eps[3]=0;
/* Horizontal, vertical limits for interpolation calc */
m1min=m10-intermod; if(m1min<0) m1min=0;
m1max=m10+1+intermod; if(m1max>xnumx-1) m1max=xnumx-1;
/**/
m2min=m20-intermod; if(m2min<0) m2min=0;
m2max=m20+1+intermod; if(m2max>ynumy-1) m2max=ynumy-1;
/**/
/* Interpolation weights calc after Fornberg (1996) */
nodewt(m1min,m1max,m2min,m2max,x,y,0,0);
/**/
/* T interpolate after interpolation weights */
for (m1=m1min;m1<=m1max;m1++)
for (m2=m2min;m2<=m2max;m2++)
	{
	/* Current node num, wt */
	m3=m1*ynumy+m2;
	ival=cn[m1-m1min][1]*cn[m2-m2min][0];
	eps[2]+=ival*tk[m3];
	eps[3]+=ival*tk2[m3];
	}
/**/
/* Wt for nodes save for T */
wn[2]=m1min; wn[3]=m1max;
for (m1=m1min;m1<=m1max;m1++)
	{
	cn[m1-m1min][5]=cn[m1-m1min][1];
	}
wn[4]=m2min; wn[5]=m2max;
for (m2=m2min;m2<=m2max;m2++)
	{
	cn[m2-m2min][4]=cn[m2-m2min][0];
	}
/* End T interpolation ------------------------ */
/**/
}
/* Calculation of T,T0 for current location by Interpolation */
/**/
/**/
/* Calculation of P on markers by interpolation from nodes */
void allinterp(double x, double y)
/* x,y - XY location of point for Vx,Vy calc */
{
/* Counters */
long int m1,m2,m3,m10,m20,m1min,m1max,m2min,m2max;
/* en - normalized distance */
double ival;
/**/
/**/
/* Check X,Y */
if(x<0) x=0; else if(x>xsize) x=xsize;
if(y<0) y=0; else if(y>ysize) y=ysize;
/**/
/**/
/**/
/* Up Left Node X,Y Num */
wn[0]=m10=m1serch(x);
wn[1]=m20=m2serch(y);
/**/
/**/
/**/
/* P interpolation ------------------------ */
/* Buffer clear */
eps[10]=0;
/* Horizontal, vertical limits for interpolation calc */
m1min=m10; if(x>(gx[m10]+gx[m10+1])/2.0) m1min+=1;
if(m1min<1) m1min=1; if(m1min>xnumx-2) m1min=xnumx-2;
m1max=m1min+1+intermod; if(m1max>xnumx-1) m1max=xnumx-1;
m1min=m1min-intermod; if(m1min<1) m1min=1;
/**/
m2min=m20; if(y>(gy[m20]+gy[m20+1])/2.0) m2min+=1;
if(m2min<1) m2min=1; if(m2min>ynumy-2) m2min=ynumy-2;
m2max=m2min+1+intermod; if(m2max>ynumy-1) m2max=ynumy-1;
m2min=m2min-intermod; if(m2min<1) m2min=1;
/**/
/* Interpolation weights calc after Fornberg (1996) */
nodewt(m1min,m1max,m2min,m2max,x,y,-1,-1);
/**/
/* EPSxx,EPSyy,P interpolate after interpolation weights */
for (m1=m1min;m1<=m1max;m1++)
for (m2=m2min;m2<=m2max;m2++)
	{
	/* Current node num, wt */
	m3=m1*ynumy+m2;
	ival=cn[m1-m1min][1]*cn[m2-m2min][0];
	eps[10]+=ival*pr[m3];
	}
/* End EPSxx,EPSyy,P interpolation ------------------------ */
/**/
}
/* Calculation of P by Interpolation */
/**/
/**/
/* Calculation of P from the depth below the surface */
void depthp(double x, double y)
/* x,y - XY location of point for Vx,Vy calc */
{
/* Counters */
long int m10;
/* en - normalized distance */
double e,n,ival;
/**/
/* Up Left Node X,Y Num */
m10=wn[0];
/**/
/* Check X,Y */
if(x<0) x=0; else if(x>xsize) x=xsize;
if(y<0) y=0; else if(y>ysize) y=ysize;
/**/
/**/
/**/
/* Pressure as function of depth from erosion level calc */
/* Relative normalized coord calc */
e=(x-gx[m10])/(gx[m10+1]-gx[m10]);
n=(e*ep[m10+1]+(1.0-e)*ep[m10]);
ival=3300.0*(y-n)*9.80665;
if (ival<1e+5) ival=1e+5;
eps[10]=ival;
/**/
/**/
}
/* Calculation of P from the depth below the surface */
/**/
/**/
/* Nu calc after rheological equation */
/* P-T-stress dependent rheology without/with brittle/ductile transition */
/* Rheological equations */
/* Stress>SScr */
/* Power law dislocation creep: SSii={NU0*EEii*exp[(E+PV)/RT]}^(1/n) */
/* Effective viscosity: NU=1/2*{NU0*exp[(E+PV)/RT]}^(1/n)*EEii^[(1-n)/n] */
/* Stress<SScr */
/* Newtonian diffusion creep: SSii=NU1*EEii*exp[(E+PV)/RT] */
/* Effective viscosity: NU=NU0/2*exp[(E+PV)/RT] */
/* NU1=NU0/SScr^(n-1) */
/* SScr - dislocation, diffusion transition stress */
/* SSii - second invariant of deviatoric stress tensor */
/* EEii - second invariant of strain rate tensor */
/* E - activation energy, J */
/* V - activation volume, J/bar */
/* R - gase constant 8.314 J/K */
/* Viscosity NU  calc after rheological equations */
/* NU=SSii/(2*EEii) */
/* Brittle - Ductile transition */
/* sbrit=MINV(0.85e+5*pb,60e+6+0.6e+5*pb)*lambda;  [Schott & Schmeling (1998)] */
/* sbrit=MINV(0.667e+5*pb,51.2e+6+0.512e+5*pb)*lambda; [Brace & Kohlstedt (1980)] */
double viscalc(double mtk, double mpb, double x, double y, long int mm1, int mm2)
/* mtk - T, K */
/* mpb - P, bar */
/* x,y - XY location of point for Vx,Vy calc */
/* mm1 - marker number */
/* mm2 - rock type */
{
/* Val buffer */
double e,n,rt=8.314*mtk,k1,e1,epsin,sduct,sbrit,nueff,smin,smax,nmin,nmax,strain,abrit,bbrit,nubrit,nunewt,nupowl;
/* Rheological Eq par */
double lamb,mnu2,A,sig0,p,q,Pepsin,nupeierls,siginnew,sigin;
/* Counters */
long int m1;
int n2;
/**/
/* Molten rocks */
if (mm2>20) return markn0[mm2];
/**/
/**/
/* Calc effective strain rate after second strain rate tensor invariant EEii=(1/2SUM(EPSik^2))^(1/2) */
epsin=markeii;
sigin=marksii;
/**/
/* Brittle "viscosity" calc: sbrit=A(eps)+LAM*B(eps)*P ------------------- */
/* Brittle/ductile transition stress calc */
lamb=markll[mm2];
abrit=marka0[mm2];
bbrit=markb0[mm2];
sbrit=abrit+bbrit*mpb*1e+5*lamb;
/* Constant limiting stress */
if(marks0[mm2]<=0 && sbrit>marks1[mm2]) sbrit=marks1[mm2];
/**/
/* Reduce plastic strength of impacted material for a certain time */
if(sbrit>0 && markhi[mm1]>0 && ((timesum-(markhi[mm1]*1.000e+6*365.250*24.000*3600.000))<(relax_time*1.000e+6*365.250*24.000*3600.000)))
	{
	sbrit=1.000e+6;     /* new, reduced strength of impacted material (pressure-independent) */
	}
/**/
/* Inverted value of Brittle "Mohr-Coulomb" viscosity calc */
nubrit=0;
if(epsin>0 && (bbrit>0 || abrit*lamb>0))
	{
	if(sbrit>0)
		{
		nubrit=1.0/(0.5*sbrit/epsin/markrii);
		}
	else
		{
		nubrit=1.0/markn0[mm2];
		}
	}
/* End brittle "viscosity" calc: sbrit=A(eps)+LAM*B(eps)*P ------------------- */
/**/
/**/
/**/
/* Peierls plasticity-creep mechanism */
nupeierls=0;
if(marks0[mm2]>0 && sigin>1e+7 && epsin>0)
	{
	/* A=pow(10.0,7.8)*1e-12=6.3e-5 Constant f(p,q), 1/s/MPa^2 */
	A=marks0[mm2];
	/* 9.1e+9  Dry Peierls stress at 0 K f(p,q), MPa [Evans & Goetze, JGR, 84, 5505-5524 (1979)] */
	/* 2.9e+9  Wet Peierls stress at 0 K f(p,q), MPa [Katayama & Karato, PEPI, 166, 57-66 (2008)] */
	sig0=marks1[mm2];
	/* Iterating Peierls stress */
	siginnew=sig0*0.95;
	for(n2=0;n2<5;n2++)
		{
		siginnew=pow(-rt/(markdh[mm2]+markdv[mm2]*mpb)*log(epsin/A/siginnew/siginnew),0.5);
		if(siginnew>1.0) siginnew=-siginnew;
		siginnew=(1.0-siginnew)*sig0;
		}
	if(siginnew>sig0) siginnew=sig0;
/**/
	nupeierls=1.0/(0.5*siginnew/epsin);
	}
/**/
/**/
/**/
/* Ductile viscosity calc -------------------------------------------*/
/* Inverted value of Newtonian NU set */
nunewt=0;
/**/
/* Inverted value of power-low NU set */
nupowl=0;
/**/
/* Check for the presence of ductile rheology */
if (marknu[mm2])
	{
	/* A)  Simple Newtonian rheology */
	/* Newtonian creep: SSii=NU0*2.0*EEii */
	/* Effective viscosity: NU=NU0 */
	/* Effective viscosity member in Stokes: NUs=NU */
	if(markdh[mm2]==0 && markdv[mm2]==0 && (markss[mm2]==0 || markmm[mm2]==1.0))
		{
		/* Inverted value of Newtonian NU calc */
		nunewt=1.0/marknu[mm2];
		}
	/**/
	/**/
	/**/
	/* B)  P-T dependent, stress independent Newtonian rheology */
	/* Newtonian diffusion creep: SSii=NU0*EEii*exp[(E+PV)/RT] */
	/* Effective viscosity: NU=NU0*exp[(E+PV)/RT] */
	if((markdh[mm2]!=0 || markdv[mm2]!=0) && (markss[mm2]==0 || markmm[mm2]==1.0))
		{
		/* Inverted value of Newtonian NU calc */
		e1=(markdh[mm2]+markdv[mm2]*mpb)/rt;
		if(e1>150.0) e1=150.0;
		/* Test creep Moresi & Solomatov (1995): SSii=NU0*exp[-a*(T-T0)] */
		if(markdh[mm2]<0 && markdv[mm2]>0)
			{
			e1=markdh[mm2]*(mtk-markdv[mm2]);
			if(e1<-150.0) e1=-150.0;
			}
		/* Test creep Turcotte & Schubert(1982): SSii=NU0*exp[E/RTo(1-(T-T0)/T0)] */
		if(markdh[mm2]<0 && markdv[mm2]<0)
			{
			e1=(-markdh[mm2])*(1.0-(mtk-(-markdv[mm2]))/(-markdv[mm2]))/8.314/(-markdv[mm2]);
			if(e1>150.0) e1=150.0;
			}
		nunewt=1.0/(marknu[mm2]*exp(e1));
		}
	/**/
	/**/
	/**/
	/* C)  P-T independent, stress dependent rheology without/with brittle/ductile transition */
	/* Stress>SScr */
	/* Power law creep: SSii={NU0*EEii}^(1/n) */
	/* Effective viscosity: NU=1/2*NU0^(1/n)*EEii^[(1-n)/n] */
	/* Effective viscosity member in Stoks: NUs=NU/n */
	/* Stress<SScr */
	/* Newtonian creep: SSii=NU1*EEii */
	/* Effective viscosity: NU=NU1/2 */
	/* Effective viscosity member in Stokes: NUs=NU */
	/* NU1=NU0/SScr^(n-1) */
	if((markdh[mm2]==0 && markdv[mm2]==0) && markss[mm2]!=0 && markmm[mm2]!=1.0)
		{
		/* Coef for stress-independent creep NU1 calc */
		k1=marknu[mm2]/pow(markss[mm2],markmm[mm2]-1.0);
		/**/
		/* Inverted value of Newtonian NU calc */
		nunewt=1.0/(0.5*k1);
		/**/
		/* Ductile power-law stress calc */
		sduct=pow(marknu[mm2]*epsin,1.0/markmm[mm2]);
		/**/
		/* Inverted value of power-law NU calc */
		if (epsin>0) nupowl=1.0/(0.5*sduct/epsin);
		}
	/**/
	/**/
	/**/
	/* D)  P-T-stress dependent rheology without/with brittle/ductile transition */
	/* Rheological equations */
	/* Stress>SScr */
	/* Power law dislocation creep: SSii={NU0*EEii*exp[(E+PV)/RT]}^(1/n) */
	/* Effective viscosity: NU=1/2*{NU0*exp[(E+PV)/RT]}^(1/n)*EEii^[(1-n)/n] */
	/* Effective viscosity member in Stoks: NUs=NU/n */
	/* Stress<SScr */
	/* Newtonian diffusion creep: SSii=NU1*EEii*exp[(E+PV)/RT] */
	/* Effective viscosity: NU=NU0/2*exp[(E+PV)/RT] */
	/* Effective viscosity member in Stokes: NUs=NU */
	/* NU1=NU0/SScr^(n-1) */
	if((markdh[mm2]!=0 || markdv[mm2]!=0) && markss[mm2]!=0 && markmm[mm2]!=1.0)
		{
		/* T-P exponent for effective NU calc */
		e1=(markdh[mm2]+markdv[mm2]*mpb)/rt;
		if(e1>150.0) e1=150.0;
		e1=exp(e1);
		/**/
		/* Coef for stress-independent creep NU1 calc */
		k1=marknu[mm2]/pow(markss[mm2],markmm[mm2]-1.0);
		/**/
		/* Inverted value of Newtonian NU calc */
		nunewt=1.0/(0.5*k1*e1);
		/**/
		/* Ductile power-law stress calc */
		sduct=pow(marknu[mm2]*e1*epsin,1.0/markmm[mm2]);
		/**/
		/* Inverted value of power law NU calc */
		if(epsin>0) nupowl=1.0/(0.5*sduct/epsin);
		}
	}
/* End ductile viscosity calc -------------------------------------------*/
/**/
/**/
/**/
/* Inverted value of effective viscosity calc, check */
nueff=nunewt+nupowl;
if(nubrit>nueff) nueff=nubrit;
if(nupeierls>nueff) nueff=nupeierls;
/**/
/* Inverted viscosity check */
if(nueff<=0) nueff=1.0/markn1[mm2];
/**/
/* Viscosity calc */
nueff=1.0/nueff;
/**/
/* Viscosity check */
if(nueff<markn0[mm2]) nueff=markn0[mm2];
if(nueff>markn1[mm2]) nueff=markn1[mm2];
if(nueff<nubeg) nueff=nubeg;
if(nueff>nuend) nueff=nuend;
/**/
/* Return calculated viscosity */
/**/
return nueff;
}
/* Nu calc after rheological equation */
/**/
/**/
/* Thermodynamic database use for ro, Cp */
void tdbasecalc(double mtk, double mpb, int mm2, long int mm1)
{
/* TD Database variables, dTK,dPB - TK, PB step for tabulation in TD database */
double H0,H1,H2,H3,R0,R1,R2,R3,G0,G1,G2,G3,W0,W1,W2,W3,n,e;
/* Val Buffers */
int n1,n2,mm3,ynpb;
double mhh0,mhh1,mdhh,maa,mwa,mnu,nueff,ival,dmwa,dmtk,dmpb,xmelt,wro,mro,mcp,mbb,mgg,mkt,mkt1,xold,kr01,kr1,kr10,xkr,krad;
long int m1=wn[0];
double sy1,e1;
/**/
/* Adiabate computing */
/**/
/**/
/* Reset TD variables */
eps[40]=eps[41]=eps[42]=eps[43]=eps[44]=eps[45]=0;
/**/
/* Thermal conductivity */
/* m895 Dry peridotite Fe=12 */
/* Olivine: Hoffmeister, 1999; Hoffmeister & Yuen, 2005 */
if(mpb<235000.0)
	{
	/* Lattice k */
	mkt1=(1.878+770.9/MINV(mtk,1200.0))*(1.0+4.26e-6*mpb);
	/* Radiative k 0.1 mm */
	kr01=pow(mtk/4000.0,3.0);
	/* Radiative k 1 mm */
	kr1=pow(mtk/1774.0,3.0);
	/* Radiative k 10 mm */
	xkr=pow(mtk/1636.0,10.0);
	xkr/=xkr+1.0; kr10=pow((mtk-1000.0*xkr)/1011.0,3.0)-0.7713*xkr;
	}
/* Perovskite: Hoffmeister, 1999; Hoffmeister & Yuen, 2005 */
else
	{
	/* Lattice k */
	mkt1=(1.291+1157.0/MINV(mtk,2100.0))*(1.0+2.50e-6*mpb);
	/* Radiative k 0.1 mm */
	kr01=pow(mtk/3591.0,3.0);
	/* Radiative k 1 mm */
	kr1=pow(mtk/2117.0,3.0);
	/* Radiative k 10 mm */
	xkr=pow(mtk/1500.0,4.0); xkr/=xkr+1.0;
	kr10=pow((mtk+4000.0*xkr)/5776.0,3.0)+2.822*xkr;
	}
krad=kr1;
/**/
/**/
/**/
/* Deep TD base type */
switch (mm2)
	{
	/* Mantle database */
	/* Silicate */
	case 5:
	case 6:
	case 25:
	case 26:
	mm3=0; break;
	/**/
	/* Crust database */
	case 11:
	case 12:
	case 31:
	case 32:
	mm3=1; break;
	/* Unknown type */
	default: {printf("Unknown rock type for TD database %d",mm2); exit(0);}
	}
/* ABCD-4Cell Number */
e=(mtk-tkmin)/tkstp;
if(e<0) e=0;
if(e>(double)(tknum-1)) e=(double)(tknum-1);
n=(mpb-pbmin)/pbstp;
if(n<0) n=0;
if(n>(double)(pbnum-1)) n=(double)(pbnum-1);
n1=(int)(e);
if(n1>tknum-2) n1=tknum-2;
n2=(int)(n);
if(n2>pbnum-2) n2=pbnum-2;
/* e,n Calc */
e=(e-(double)(n1));
n=(n-(double)(n2));
/* Ro H values */
/* 0 2 */
/* 1 3 */
R0=td[n1  ][n2  ][mm3][0]*1000.0;
R1=td[n1  ][n2+1][mm3][0]*1000.0;
R2=td[n1+1][n2  ][mm3][0]*1000.0;
R3=td[n1+1][n2+1][mm3][0]*1000.0;
H0=td[n1  ][n2  ][mm3][1]*1000.0*4.1837;
H1=td[n1  ][n2+1][mm3][1]*1000.0*4.1837;
H2=td[n1+1][n2  ][mm3][1]*1000.0*4.1837;
H3=td[n1+1][n2+1][mm3][1]*1000.0*4.1837;
/* Ro calc by interpolation */
mro=((R0*(1.0-n)+R1*n)*(1.0-e)+(R2*(1.0-n)+R3*n)*e);
/* Cp calc by interpolation */
mcp=((H2-H0)*(1.0-n)+(H3-H1)*n)/tkstp;
if(mcp<1e+2) mcp=1e+2; else if(mcp>5e+4) mcp=5e+4;
/* Effective adiabatic betta=1/V*dV/dT=ro/T*[-dH/dP+V] calc by interpolation */
mbb=(2.0/(R1+R0)-(H1-H0)/pbstp/1e+5)*(1.0-e)+(2.0/(R3+R2)-(H3-H2)/pbstp/1e+5)*e;
mbb*=mro/mtk;
if(mbb<-1e-2) mbb=-1e-2; else if(mbb>1e-2) mbb=1e-2;
/* Effective compressibility term alpha=1/ro*d(ro)/dP calc by interpolation */
maa=(2.0/(R1+R0)*(R1-R0)*(1.0-e)+2.0/(R3+R2)*(R3-R2)*e)/pbstp/1e+5;
if(maa<0) maa=0;
/* Thermal conductivity */
mkt=mkt1+krad;
/* Viscosity */
mnu=eps[48];
/**/
/* Partially molten rocks */
meltpart1(mtk,mpb,mm2);
if(eps[21]>0 && mm2<20) markt[mm1]+=20;
if(eps[21]<=0 && mm2>20) markt[mm1]-=20;
mm2=markt[mm1];
/* Check marker type */
if (mm2>20)
	{
	/* Calculate melt fraction */
	xmelt=eps[21];
	/**/
	/* Density */
	mro=mro*(1.0-xmelt)+xmelt*markro[mm2]*(1.0-markbb[mm2]*(mtk-298.15))*(1.0+markaa[mm2]*(mpb-1.0)*1e-3);
	/**/
	/* Viscosity */
	/* Effective NU calc check */
	if(xmelt>0.1)
		{
		/* Significant melt */
		/* Effective NU calc check [Pinkerton & Stevenson, Volc. Geotherm. Res., 53, 47-66 (1992) eq. (5)], bug fixed by Greg (20/12/2010) */
		nueff=marknu[mm2]*exp((2.500+pow(((1.000-xmelt)/xmelt),0.480))*(1.000-xmelt));
		/* Viscosity check */
		if(nueff<markn0[mm2]) nueff=markn0[mm2];
		if(nueff>markn1[mm2]) nueff=markn1[mm2];
		if(nueff<nubeg) nueff=nubeg;
		if(nueff>nuend) nueff=nuend;
		/**/
		mnu=nueff;
		}
	/**/
	/* Heat conductivity */
	mkt=((markkt[mm2]+markkf[mm2]/(mtk+77.0))*exp(markkp[mm2]*mpb))*xmelt+((markkt[mm2-20]+markkf[mm2-20]/(mtk+77.0))*exp(markkp[mm2-20]*mpb))*(1.0-xmelt);
	/**/
	/* Additional melting adiabatic term, heat capacity */
	if(xmelt>0 && xmelt<1.0)
		{
		/* Melting adiabatic term: alm=-ro*(dHlat/dP)/T */
		/* Numerical differentiation */
		dmpb=mpb*0.001;
		meltpart1(mtk,mpb-dmpb,mm2);
		ival=eps[22];
		meltpart1(mtk,mpb+dmpb,mm2);
		ival-=eps[22];
		ival*=mro/(mtk*2.0*dmpb*1e+5);  /* Bug fixed by Greg (01/06/2011) */
		mbb+=ival;
		/**/
		/* Melting heat capacity term: cpm=dHlat/dT */
		/* Numerical differentiation */
		dmtk=1.0;
		meltpart1(mtk+dmtk,mpb,mm2);
		ival=eps[22];
		meltpart1(mtk-dmtk,mpb,mm2);
		ival-=eps[22];
		ival/=2.0*dmtk;
		mcp+=ival;
		/**/
		}
	}
/* Save TD variables */
eps[41]=mro;
eps[43]=mcp;
eps[44]=mbb;
eps[45]=maa;
eps[47]=mkt;
eps[48]=mnu;
/**/
}
/* Thermodynamic database use for ro, Cp */


/* Melt fraction, latent heat calculation */
void meltpart1(double mtk, double mpb, int mm2)
/* mtk - T, K */
/* mpb - P, bar */
/* x,y - XY location of point for Vx,Vy calc */
/* mm1 - mark number */
/* mm2 - mark type */
/* yn  - type of calculation: 0 - Ro, 1 - Nu, 2 - Cp, 3 - kt */
{
/* Val buffer */
double xmelt=0,xmelt_fe=0,hlatent=0,hlatent_fe=0,ival;
long int m1;
double ykm=mpb*3e-3,ts=0,tl=0,tl_fe=0;
/**/
/**/
/* Calculate melt fraction using marker type */
switch(mm2)
	{
	/* Basalt, Gabbro: latent heat 380 kJ/kg [Bittner and Schmeling, Geophys. J. Int., 123, 59-70 (1995)] */
	/* Crust */
	case 11:
	case 12:
	case 31:
	case 32:
	/* Wet solidus [Schmidt and Poli (1998)]  */
	if (ykm<48.0)
		{
		ts=972.6-2111.0/(ykm+10.63)+70033.0/(ykm+10.63)/(ykm+10.63);
		}
	else
		{
		ts=935.4+0.1162*ykm+0.006937*ykm*ykm;
		}
	/* Dry Toleitic Basalt Liquidus [Hess, 1989] */
	tl=1423.15+3.5*ykm;
	hlatent=380000.0;
	break;
	/**/
	/* Dry Peridotite: latent heat 400 kJ/kg [Turcotte & Schubert, p. 171 (1982)] */
	case 5:
	case 6:
	case 25:
	case 26:
/**/
	if(si_melting == 0)
		{
		/* Dry peridotite solidus [Hirschmann, Geochem. Geophys. Geosyst., 1, 2000GC000070 (2000)] */
		if(mpb<10e+4)        /* P < (10e4 bar = 1e5 bar = 1e10 Pa =) 10 GPa */
			{
			ts=273.15+1120.661+132.899e-4*mpb-5.104e-8*mpb*mpb;
			}
		else
		/* Linear extrapolation to higher pressures */
			{
			ts=273.15+1939.251+30.819e-4*(mpb-10e+4);
			}
		}
/**/
	else if(si_melting == 1)
		{
		/* Dry peridotite solidus [Herzberg et al., Geochem. Geophys. Geosyst., 1, 2000GC000089 (2000)] */
		if(mpb<21.50e+4)        /* P < (21.5e4 bar = 2.15e5 bar = 2.15e10 Pa =) 21.5 GPa */
			{
			ts=273.15+1143.04342+58.2946423e-4*mpb+52.3439318e-8*mpb*mpb-16.3201032e-12*mpb*mpb*mpb+2.29886314e-16*pow(mpb,4.000)-0.180865486e-20*pow(mpb,5.000)+0.00815679773e-24*pow(mpb,6.000)-0.000197104325e-28*pow(mpb,7.000)+1.97908526e-38*pow(mpb,8.000);
			}
		else
		/* Linear extrapolation to higher pressures */
			{
			ts=273.15+2157.500+11.7297e-4*(mpb-21.50e+4);
			}
		}
/**/
	/* Martian mantle peridotite liquidus parametrized after [Wade and Wood, EPSL, 236, 78-95 (2005)] */
	tl=1973.000+28.570e-4*mpb;
	hlatent=400000.0;  /* in [J/kg] */
/**/
	break;
/**/
/* Iron melting by Gregor */
	case 7:
	case 8:
	case 9:
	case 10:
	case 17:
	case 18:
	case 19:
/**/
	if(fe_melting == 0)
		{
		/* Pure iron melting curve parametrized after [Boehler, Nature, 363, 534-536 (1993)] */
		/**/
		if(mpb>=0.000e6 && mpb<0.100e6)            /*  0 GPa <= P < 10 GPa */
			{
			tl_fe=1761.000+3.100e-3*mpb;
			}
		if(mpb>=0.100e6 && mpb<0.200e6)            /* 10 GPa <= P < 20 GPa */
			{
			tl_fe=1863.000+2.080e-3*mpb;
			}
		if(mpb>=0.200e6 && mpb<0.600e6)            /* 20 GPa <= P < 60 GPa */
			{
			tl_fe=2071.800+1.035e-3*mpb;
			}
		if(mpb>=0.600e6 && mpb<1.000e6)            /* 60 GPa <= P < 100 GPa */
			{
			tl_fe=2382.800+5.170e-4*mpb;
			}
		if(mpb>=1.000e6)                           /* P >= 100 GPa */
			{
			tl_fe=2900.000+(1.000/1.000e3)*mpb;
			}
		}
/**/
	else if(fe_melting == 1)
		{
		/* Eutectic Fe-FeS melting curve using data compilation by [Chudinovskikh & Boehler, EPSL, 257, 97-103 (2007)] */
		/**/
		if(mpb<=0.4180e6)        /* P < (0.418e6 bar = 4.18e5 bar = 4.18e10 Pa =) 41.8 GPa */
			{
			tl_fe=1260.1000+3.3171e-4*mpb-3.395e-9*mpb*mpb+2.660e-14*mpb*mpb*mpb-3.7688e-20*pow(mpb,4.000);
			}
		else
		/* Linear extrapolation to higher pressures */
			{
			tl_fe=1597.7000+4.2635e-4*(mpb-0.4180e6);
			}
		}
/**/
	hlatent_fe=240000.0;  /* Latent heat of pure iron in [J/kg] [Sramek et al., Geophys. J. Int., 181, 198-220 (2010)] */
	break;
	/**/
        /**/
	/* Other rocks - No melting */
	default:
	break;
	}
/**/
/* Silicate melt fraction, latent heat calculation */
eps[21]=eps[22]=0;
if(tl)
	{
	/* Silicate melt fraction calc, check */
	if(ts>=tl-100.0) ts=tl-100.0;
	xmelt=(mtk-ts)/(tl-ts);
	if(xmelt<0.0000) xmelt=0.0000;
	if(xmelt>1.0000) xmelt=1.0000;
/**/
/* Save silicate melt fraction and latent heat */
	eps[21]=xmelt;
	/* Latent heat calc for silicates */
	hlatent*=xmelt;
	eps[22]=hlatent;
	}
/**/
/* Iron melt fraction */
eps[50]=eps[51]=0;
if(tl_fe)
	{
	/* Iron melt fraction calc, check */
	if(mtk<tl_fe) xmelt_fe=0.0000;
	if(mtk>=tl_fe) xmelt_fe=1.0000; /* For eutectic iron solidus = liquidus [Andrault et al., High Pressure Res., 26, 267-276 (2006)] */
/**/
/* Save iron melt fraction */
	eps[50]=xmelt_fe;
	/* Latent heat calc for iron */
	hlatent_fe*=xmelt_fe;
	eps[51]=hlatent_fe;
	}
/**/
}
/* Melt fraction, latent heat calculation */


