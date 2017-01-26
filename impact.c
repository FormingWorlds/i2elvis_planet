/* New subroutine to introduce parametrized planetary accretion */
/* by Gregor J. Golabek (last modified: 19/01/2016) */
/**/
/**/
/* Global variables for impacts */
double hydr_frac,impact_time[1000],impact_angle[1000],impact_mass[1000],impact_vel[1000],iron_type[1000],spin_rate[1000],fe_frac,start_time;
int    no1;
long int impactnnn;
/**/
/**/
/* Read out all data needed about accretion */
int impactread()
{
long int m1;
/**/
	/* Read impact_history.t3c in every timestep, so code knows even after restart for which impact it has to look */
/**/
	fl1 = fopen("impact_no.t3c","rt");
        ffscanf1(); no1=atoi(sa);                          /* impactor no. [non-dim.] */
	fclose(fl1);
/**/
	/* Read data of all impact events */
	printf("Read data of all impactor bodies \n");
/**/
	fl1 = fopen("impact_history.t3c","rt");
	ffscanf1(); impactnnn  = atoi(sa);                 /* total no. of impactors [non-dim.] */
	ffscanf1(); start_time = atof(sa);                 /* start time of model after CAI formation [Ma] */
	ffscanf1(); fe_frac    = atof(sa);                 /* linear iron fraction of impactor bodies [non-dim.] */
/**/
	for(m1=1;m1<=impactnnn;m1++)
		{
		ffscanf1();impact_time[m1]   = atof(sa);   /* impact time [Ma] */
		ffscanf1();impact_angle[m1]  = atof(sa);   /* impact angle [degrees] */
		ffscanf1();impact_mass[m1]   = atof(sa);   /* mass of impactor [kg] */
		ffscanf1();impact_vel[m1]    = atof(sa);   /* impact velocity [m/s] */
		ffscanf1();spin_rate[m1]     = atof(sa);   /* spin rate of target body AFTER collision [1/a] */
		ffscanf1();iron_type[m1]     = atoi(sa);   /* iron marker number for each impact [non-dim.] */	
		}
	fclose(fl1);
/**/
/**/
	/* Read out how much mass was already accreted by the planetary body, what was it's initial mass */
	/* Useful that is is stored in external file as it makes sure that data are available after code restart */
/**/
	fl1 = fopen("impact_mass.t3c","rt");
	ffscanf1(); M_acc  = atof(sa);            /* atof reads out real value [kg] */
	fclose(fl1);
/**/
	/* Read out initial mass of the target body before starting */
/**/
	fl1 = fopen("start_mass.t3c","rt");
	ffscanf1(); M_init = atof(sa);            /* atof reads out real value [kg] */
	fclose(fl1);
/**/
	/* Read impactor core mass and the thickness of the resulting iron ponds */
	fl1 = fopen("core_stuff.t3c","rt");
	ffscanf1(); no_merger    = atoi(sa);        /* current number of completed iron core-iron diapir mergers [non-dim.] */
        ffscanf1(); core_radius  = atof(sa);        /* current radius of the preexistent iron core [m] */
        ffscanf1(); M_core       = atof(sa);        /* current mass of the preexistent iron core [kg] */
        ffscanf1(); M_diapir     = atof(sa);        /* mass of the iron diapirs from the latest impactor body [kg] */        
        ffscanf1(); M_diapir_ol  = atof(sa);        /* mass of iron diapirs from the previous impactor body [kg] */
        ffscanf1(); M_diapir_old = atof(sa);        /* mass of iron diapirs from second to last impactor body [kg] */
        ffscanf1(); d_iron       = atof(sa);        /* max. thickness of iron pond induced by latest impact [m] */
	ffscanf1(); d_iron_ol    = atof(sa);        /* max. thickness of iron pond induced by previous impact [m] */
        ffscanf1(); d_iron_old   = atof(sa);        /* max. thickness of iron pond induced by second to last impact [m] */
	fclose(fl1);
/**/	
return 0;
}
/* End read impact history */
/**/
/**/
/* Start save impact and planetary dynamo history */
int impactsave()
{
	/* Impact number will be set to current value in file "impact_no.t3c" */
	/* Makes sure that current number of impacts is not forgotten after restart */
	fl1 = fopen("impact_no.t3c","wt");
	fprintf(fl1,"%d \n",no1);                /* d stands for integers */
	fclose(fl1);
/**/
	/* Current mass of accreted material [kg] */
	fl1 = fopen("impact_mass.t3c","wt");
	fprintf(fl1,"% 25.0f",M_acc);            /* f stands for double precision */
	fclose(fl1);
/**/
	/* Save the current core radius of the target body [m] */
	fl1 = fopen("core_stuff.t3c","wt");
	fprintf(fl1,"%d %f %e %e %e %e %f %f %f \n",no_merger,core_radius,M_core,M_diapir,M_diapir_ol,M_diapir_old,d_iron,d_iron_ol,d_iron_old);
	fclose(fl1);
/**/
	/* Save new data regarding planetary dynamo history */
	/* Output data: */
        /* Current time [Ma] */ 
        /* Current core radius [m] */
        /* Mean core temperature [K] */ 
        /* Mean temperature of the mantle CMB [K] */ 
        /* Adiabatic core heat flux [W/m^2] */ 
        /* CMB heat flux at mantle side [W/m^2] */ 
        /* Magnetic Reynolds number [non-dim.] */ 
        /* Dynamo action parameter [non-dim.] */ 
        /* Magnetic dipole moment [A*m^2] */ 
        /* Local Rossby number [non-dim.] */
	fl1 = fopen("core_dynamo.t3c","a+");     /* a+ stands for adding new data at the end of preexistent file */
	fprintf(fl1,"%f %f %f %f %f %f %f %d %e %f \n",timesum/(3600.000*24.000*365.250*1.000e6),core_radius,t_mean,t_cmb_mean,q_core,q_mantle,reynolds,dynamo,dipl_mom,rossby);
	fclose(fl1);
/**/
	/* Save data regarding evolution of hydrous silicates */
	/* Output data: */
        /* Current time [Ma] */ 
        /* Current hydrous fraction [non-dim.] */
	fl1 = fopen("hydrous_silicates.t3c","a+");     /* a+ stands for adding new data at the end of preexistent file */
	fprintf(fl1,"%f %f \n",timesum/(3600.000*24.000*365.250*1.000e6),hydr_frac);
	fclose(fl1);
/**/
return 0;
}
/* End save impact and planetary dynamo history */
/**/
/**/
/*===========================================================================================================================*/
/**/
/* Now handle the impact event itself */
int impact()
{
/* Counters */
long int air_marker_no,count,count1,ii,ij,ij_max,m1,m2,m3,m4,m5,m6,m7;
long int mm1,mm2,mmm1,minmx=0,maxmx=0,minmy=0,maxmy=0,mm3,si_marker_no;
/* Nonstability for markers, m */
int nonstab;
double pival,wlen,wamp,wymin,wymid,wymax,xnonstab=0,ynonstab=0;
/* Distance val */
double x,y,x1,y1,x2,y2,x3,y3,x4,y4,dx,dy;
/* Initial Rock type for grid, rock body numbering, initial distribution */
double bx[4],by[4];
double cnst,koef,koef1,koef2,ival,ival1,akf,bkf;
long int m10,m11,m20,m21,nshift,nshift1,nshift2;
/**/
/* Define all additional parameters needed */
int    marker_counter,marker_counter_all,no_left,no_right,no_ic,no_it,no_ite,no_iter,no_iter_max;
int    stack_counter,status,status2,status3,status4,status_phi_left,status_phi_right;
int    status_layer_left,status_layer_right,status_pond_left,status_pond_right;
double angle,angle_far,angle1,angle4,angle5,angle6,angle11,angle22,corr,corr1,current_radius;
double c_p,crater_angle,crater_radius,crater_radius_corr,crit_density,crit_density_fe,crit_density_si;
double deltax,deltay,deltaxcirc,deltaycirc,defl_corr,deflection,d_layer,dist;
double e_al,e_fe,e_k,e_th,e_u,e_uu,e_diff,e_undiff,ejecta_left,ejecta_right,eros,eros1,eros2;
double escape_vel,f,fe_fract,form_time,g,G,height_left,height_right;
double impact1,impact2,impact3,impact4,impact5,impact6,impact_core_radius,impact_radius,impact_radius_corr;
double impact_far,impact_left,impact_right,impact_lleft,impact_rright,initial_radius,intervall;
double iso_core,isobaric_core,layer_min_left,layer_min_right,layer_max_left,layer_max_right,length_left,length_right;
double marker_temp_stack,marker_ratio,marker_mean_temp,mass_ejecta,mass_imp,mass_imp_core,mass_imp_mantle,mass_tot;
double offset,offsetcirc;
double phi_angle,phi_error,phi_angle_old,phi_max_left,phi_max_right,phi_min_left,phi_min_right;
double planet_angle,planet_angle1,planet_curr,planet_current[1000],planet_far;
double planet_leftmost,planet_rightmost,planet_lleftmost,planet_rrightmost;
double planet_low,planet_radius,planet_surf,planet_surface[1000],pond_angle_left,pond_angle_right;
double pond_min_left,pond_min_right,pond_max_left,pond_max_right;
double radius_new,radius_neww,radius_new_far,range,rho_blanket,rho_ejecta,rho_impact;
double rho_co,rho_ic,rho_icc,rho_ma,rho_planet,rho_Si,rho_Si_impact,rho_Si_melt,rho_Fe,rho_Fe_pond,rho_Fe_impact;
double search_radius,Si_melt,Si_melt_blanket,Si_melt_imp,stack,surface_left,surface_right;
double t0,t00,t_blanket,t_ejecta,t_ejectaa,t_form,t_form2,t_left,t_leftt,t_right,t_rightt,t_liq,t_sol,t_space,t_test;
double testradius,testradius_far,testradius1,testradius3,testradius4,testradius5,testradius11,testradius22,test_time;
double t_impactor_core,t_impactor_core_old,t_impactor_mantle,t_impactor_mantle_old,time_incr,v_ratio;
double xcoord_air_marker,ycoord_air_marker,xcoord_si_marker,ycoord_si_marker;
double xcoordcirc,ycoordcirc,xcoordiron,ycoordiron;
double xcoordsurf,ycoordsurf,xcoord_left,ycoord_left,xcoord_right,ycoord_right;
double xcoord_lleft,ycoord_lleft,xcoord_rright,ycoord_rright;
/**/
/**/
/* Set or read out important physical constants */
pival           = 3.141592654;   /* Set value of Pi */
G               = 6.67384e-11;   /* Gravitational constant [m^3/(kg*s^2)] */
rho_Fe          = markro[7];     /* Iron density from init.t3c file [kg/m^3] */
rho_Si          = markro[5];     /* Solid silicate density from init.t3c file [kg/m^3] */
rho_Si_melt     = markro[25];    /* Molten silicate density from init.t3c file [kg/m^3] */
c_p             = markcp[5];     /* Heat capacity for solid silicates [J/(kg*K)] (also valid for molten silicates and iron) */
t_space         = markk[1];      /* Temperature of space (sticky air markers), take for this marker no. 1, which is for sure in the sticky air zone [K] */
intervall       = 1.000e0;       /* Width of sector when setting ejecta blanket [degrees] */
range           = 0.100e0;       /* Allowed uncertainty for positioning of impact ejecta [degrees] */
crit_density    = 1500.000;      /* Critical density: Which grids are counted as part of the planetary body? */
crit_density_fe = 6000.000;      /* Critical density, grids with higher density are most likely in iron core of target body [kg/m^3] */
crit_density_si = 4500.000;      /* Critical density, grids with higher density are most likely iron enriched, not pure silicates [kg/m^3] */
/**/
/* Peridotite solidus temperature [K] assuming P = 0 GPa */
if(si_melting == 0)             /* [Hirschmann, Geochem. Geophys. Geosyst., 1, 2000GC000070 (2000)] */
	{
	t_sol = 1393.811;
	}
else if(si_melting == 1)        /* [Herzberg et al., Geochem. Geophys. Geosyst., 1, 2000GC000089 (2000)] */
	{
	t_sol = 1416.193;
	}
/* Peridotite liquidus temperature [K] parametrization for P = 0 GPa [Wade and Wood, EPSL, 236, 78-95 (2005)] */
t_liq = 1973.000;
/**/
/**/
/* Start the impact loop: */
/**/
	/* Include impact feature at the correct time into the model */
        if((timesum+timestep)>=(impact_time[no1]*1.000e+6*365.250*24.000*3600.000) && no1<=impactnnn)
	{
	printf("impact event! \n");
/*===========================================================================================================================*/
/**/
	/* Checking for the initial/current radii of the planetary body depending on the sector */
/**/
	/* Set initial values */
	ij             = 0;
	ij_max         = 0;
	current_radius = 0.000000;
	planet_radius  = 0.000000;
	planet_angle1  = 0.000000;
/**/
/**/
	/* Make sure not to count the zero position twice */
	for(planet_angle1=0.000000;planet_angle1<=359.000;planet_angle1=planet_angle1+intervall)
		{
		/* Start to count */
		ij = ij+1;
/**/
		/* Set initially the distance from the planetary center to the surface to zero */
       		planet_surface[ij] = 0.000000;   /* planet_surface is non-dim.! */
		testradius1        = 0.000000;
/**/
		/* Go through all markers */
		for(mm1=0;mm1<marknum;mm1++)
			{
			/* Check only solid and molten silicates and iron, stabilizing material and sticky air are excluded */
			/* To calculate initial radius of the planetary body */
        		if(markt[mm1]==5 || markt[mm1]==6 || markt[mm1]==7 || markt[mm1]==8 || markt[mm1]==9 || markt[mm1]==10 || markt[mm1]==17 || markt[mm1]==18 || markt[mm1]==19 || markt[mm1]==25 || markt[mm1]==26)
				{
				/* Impact angle 0 <= alpha < 90 degrees */
                       		if((planet_angle1>=0.000 && planet_angle1<90.000) && (markx[mm1]/xsize>=0.500) && (marky[mm1]/ysize<0.500))
					{
                               		angle1 = atan((markx[mm1]/xsize-0.500)/(0.500-marky[mm1]/ysize))*180.000/pival;
                               		if(ABSV(angle1-planet_angle1)<=range)
						{
						testradius1=pow(((markx[mm1]/xsize-0.500)*(markx[mm1]/xsize-0.500)+(marky[mm1]/ysize-0.500)*(marky[mm1]/ysize-0.500)),0.500);
                       				planet_surface[ij]=MAXV(testradius1,planet_surface[ij]);
						}
					}
				/* Impact angle 90 <= alpha < 180 degrees */
                       		else if((planet_angle1>=90.000 && planet_angle1<180.000) && (markx[mm1]/xsize>0.500) && (marky[mm1]/ysize>=0.500))
					{
                               		angle1 = 90.000+atan((marky[mm1]/ysize-0.500)/(markx[mm1]/xsize-0.500))*180.000/pival;
                               		if(ABSV(angle1-planet_angle1)<=range)
						{
						testradius1=pow(((markx[mm1]/xsize-0.500)*(markx[mm1]/xsize-0.500)+(marky[mm1]/ysize-0.500)*(marky[mm1]/ysize-0.500)),0.500);
                       				planet_surface[ij]=MAXV(testradius1,planet_surface[ij]);
						}
					}
/**/
				/* Impact angle 180 <= alpha < 270 degrees */
                       		else if((planet_angle1>=180.000 && planet_angle1<270.000) && (markx[mm1]/xsize<=0.500) && (marky[mm1]/ysize>0.500))
					{
                               		angle1 = 180.000+atan((0.500-markx[mm1]/xsize)/(marky[mm1]/ysize-0.500))*180.000/pival;
                               		if(ABSV(angle1-planet_angle1)<=range)
						{
						testradius1=pow(((markx[mm1]/xsize-0.500)*(markx[mm1]/xsize-0.500)+(marky[mm1]/ysize-0.500)*(marky[mm1]/ysize-0.500)),0.500);
                        			planet_surface[ij]=MAXV(testradius1,planet_surface[ij]);
						}
					}
				/* Impact angle 270 <= alpha < 360 degrees */
                        	else if((planet_angle1>=270.000 && planet_angle1<360.000) && (markx[mm1]/ysize<0.500) && (marky[mm1]/ysize<=0.500))
					{
                                	angle1 = 270.000+atan((0.500-marky[mm1]/ysize)/(0.500-markx[mm1]/xsize))*180.000/pival;
                                	if(ABSV(angle1-planet_angle1)<=range)
						{
						testradius1=pow(((markx[mm1]/xsize-0.500)*(markx[mm1]/xsize-0.500)+(marky[mm1]/ysize-0.500)*(marky[mm1]/ysize-0.500)),0.500);
                        			planet_surface[ij]=MAXV(testradius1,planet_surface[ij]);
						}
					}
				}
			}
		planet_radius+=planet_surface[ij];   /* Calculate sum of all radii */
		}
	ij_max = ij;
/**/
	if(no1==1)
		{
		initial_radius = planet_radius*xsize/(double)(ij_max);    /* Calculate mean (dim.!!!) radius of the planetary body */
        	current_radius = planet_radius/(double)(ij_max);          /* Calculate mean (non-dim.!!!) radius of the planetary body before first impact */
		}
/**/
	if(no1>=2)
		{
		current_radius = planet_radius/(double)(ij_max);          /* Calculate mean (non-dim.!!!) radius of the planetary body */
		}
/**/
/*===========================================================================================================================*/
/**/
	/* Find sticky air entrained into the target body BEFORE new impact takes place and convert it into silicates */    
	if(no1>=2)
		{
		for(mm1=0;mm1<marknum;mm1++)
			{
			/* Define search radius [m] */
			search_radius      = 25000.0000;
/**/
			/* Double the search radius, if result is unsatisfactory */
			change_radius: search_radius = 2.000*search_radius;
/**/
			/* Set up values */
			air_marker_no      = 0;            /* Set number of sticky air marker to zero */
			xcoord_air_marker  = 0.00000;      /* Set x coordinate of sticky air marker to zero */
			ycoord_air_marker  = 0.00000;      /* Set y coordinate of sticky air marker to zero */
			marker_counter     = 0;            /* Set marker counter for non-sticky air markers to zero */
			marker_counter_all = 0;            /* Set marker counter for all markers to zero */
			marker_temp_stack  = 0.00000;      /* Set temperature stack to zero */
			marker_mean_temp   = 0.00000;      /* Mean temperature of non sticky air markers within the search radius */
			marker_ratio       = 0.00000;      /* Ratio between silicate+iron and all markers in search area */
/**/
			if(markt[mm1]==0)  /* Sticky air marker is found */
				{
				/* Calculate distance of sticky air marker from center of planetary body */
				dist = pow(((markx[mm1]/xsize-0.500)*(markx[mm1]/xsize-0.500)+(marky[mm1]/ysize-0.500)*(marky[mm1]/ysize-0.500)),0.500);
/**/
				/* Take only markers into consideration which are BENEATH mean PRE-impact surface */
				if((dist-current_radius)<0.0000)
					{
					/* Probable entrained sticky air marker */
					xcoord_air_marker = (double)(markx[mm1])/xsize;  /* non-dim. */
					ycoord_air_marker = (double)(marky[mm1])/ysize;  /* non-dim. */
					air_marker_no     = (int)(mm1);
/**/
					/* Find neighbouring markers within the search radius */
					for(mmm1=0;mmm1<marknum;mmm1++)
						{
						x=(double)(markx[mmm1])/xsize;  /* non-dim. */
						y=(double)(marky[mmm1])/ysize;  /* non-dim. */
						x1=x-xcoord_air_marker;
						y1=(y-ycoord_air_marker)*ysize/xsize;
/**/
						/* Search region extends to a certain radius */
						if((x1*x1+y1*y1)<=((search_radius/xsize)*(search_radius/ysize)))
							{
							/* Lower boundary for the search is position of the sticky air marker */
							if((x1*x1+y1*y1)>=(0.0000*0.0000*ysize*ysize/xsize/xsize))
								{
								/* Count all markers in the surrounding */
								marker_counter_all+=1;
/**/
								/* Only for non-sticky air markers (solid and molten silicates or iron) */
								if(markt[mmm1]>=5)
									{
									marker_counter+=1;
									marker_temp_stack+=markk[mmm1];
									}
								}
							}
						}
/**/
					if(marker_counter_all<=1) goto change_radius;  /* If search radius is devoid of any other markers, increase search radius */
/**/
					if(marker_counter_all>1)
						{
						marker_ratio     = (double)(marker_counter)/(double)(marker_counter_all);
						}
/**/
					/* If the sticky air marker is dominantly surrounded by non-sticky air markers, convert it into silicates */
					if(marker_ratio>=0.8500)
						{
						/* Converting sticky air marker into silicate marker */
						markt[air_marker_no] = 6;
/**/
						/* Former sticky air marker is assigned the mean temperature of neighbour non-sticky air markers */
						marker_mean_temp = marker_temp_stack/(double)(marker_counter);
/**/
						markk[air_marker_no] = marker_mean_temp;
						}
					/* Close distance loop */
					}
				/* Close sticky air loop */
				}
			/* Close outer marker loop */
			}
		/* Close (no1>=2) loop */
		}        
/**/
/*===========================================================================================================================*/
/**/
	/* Checking for the current impact site */
/**/
	/* Set initially the distance from the planetary center to the surface to zero */
    	planet_surf = 0.000000;
	xcoordsurf  = 0.000000;
	ycoordsurf  = 0.000000;
	testradius  = 0.000000;
/**/
	/* Go through all markers */
	for(mm1=0;mm1<marknum;mm1++)
		{
/**/
		/* Check only solid and molten silicates and iron, stabilizing material and sticky air are excluded */
        	if(markt[mm1]==5 || markt[mm1]==6 || markt[mm1]==7 || markt[mm1]==8 || markt[mm1]==9 || markt[mm1]==10 || markt[mm1]==17 || markt[mm1]==18 || markt[mm1]==19 || markt[mm1]==25 || markt[mm1]==26)
			{
			/* Impact angle 0 <= alpha < 90 degrees */
                        if((impact_angle[no1]>=0.000 && impact_angle[no1]<90.000) && (markx[mm1]/xsize>=0.500) && (marky[mm1]/ysize<0.500))
				{
                                angle = atan((markx[mm1]/xsize-0.500)/(0.500-marky[mm1]/ysize))*180.000/pival;
                                if(ABSV(angle-impact_angle[no1])<=range)
					{
					testradius=pow(((markx[mm1]/xsize-0.500)*(markx[mm1]/xsize-0.500)+(marky[mm1]/ysize-0.500)*(marky[mm1]/ysize-0.500)),0.500);
                        		planet_surf=MAXV(testradius,planet_surf);
					if(testradius==planet_surf)
						{
						/* Global coordinates of impact site */
                                		xcoordsurf = markx[mm1]/xsize;
                                		ycoordsurf = marky[mm1]/ysize;
						}
					}
				}
			/* Impact angle 90 <= alpha < 180 degrees */
                        else if((impact_angle[no1]>=90.000 && impact_angle[no1]<180.000) && ((markx[mm1]/xsize)>0.500) && ((marky[mm1]/ysize)>=0.500))
				{
                                angle = 90.000+atan((marky[mm1]/ysize-0.500)/(markx[mm1]/xsize-0.500))*180.000/pival;
                                if(ABSV(angle-impact_angle[no1])<=range)
					{
					testradius=pow(((markx[mm1]/xsize-0.500)*(markx[mm1]/xsize-0.500)+(marky[mm1]/ysize-0.500)*(marky[mm1]/ysize-0.500)),0.500);
                        		planet_surf=MAXV(testradius,planet_surf);
					if(testradius==planet_surf)
						{
						/* Global coordinates of impact site */
                                		xcoordsurf = (double)(markx[mm1])/xsize;
                                		ycoordsurf = (double)(marky[mm1])/ysize;
						}
					}
				}
/**/
			/* Impact angle 180 <= alpha < 270 degrees */
                        else if((impact_angle[no1]>=180.000 && impact_angle[no1]<270.000) && (markx[mm1]/xsize<=0.500) && (marky[mm1]/ysize>0.500))
				{
                                angle = 180.000+atan((0.500-markx[mm1]/xsize)/(marky[mm1]/ysize-0.500))*180.000/pival;
                                if(ABSV(angle-impact_angle[no1])<=range)
					{
					testradius=pow(((markx[mm1]/xsize-0.500)*(markx[mm1]/xsize-0.500)+(marky[mm1]/ysize-0.500)*(marky[mm1]/ysize-0.500)),0.500);
                        		planet_surf=MAXV(testradius,planet_surf);
					if(testradius==planet_surf)
						{
						/* Global coordinates of impact site */
                                		xcoordsurf = (double)(markx[mm1])/xsize;
                                		ycoordsurf = (double)(marky[mm1])/ysize;
						}
					}
				}
			/* Impact angle 270 <= alpha < 360 degrees */
                        else if((impact_angle[no1]>=270.000 && impact_angle[no1]<360.000) && (markx[mm1]/xsize<0.500) && (marky[mm1]/ysize<=0.500))
				{
                                angle = 270.000+atan((0.500-marky[mm1]/ysize)/(0.500-markx[mm1]/xsize))*180.000/pival;
                                if(ABSV(angle-impact_angle[no1])<=range)
					{
					testradius=pow(((markx[mm1]/xsize-0.500)*(markx[mm1]/xsize-0.500)+(marky[mm1]/ysize-0.500)*(marky[mm1]/ysize-0.500)),0.500);
                        		planet_surf=MAXV(testradius,planet_surf);
				        if(testradius==planet_surf)
						{
						/* Global coordinates of impact site */
                                		xcoordsurf = (double)(markx[mm1])/xsize;
                                		ycoordsurf = (double)(marky[mm1])/ysize;
						}
					}
				}
			/* End of silicate/iron if loop */
			}
		/* End of marker loop */
		}
/**/
/*===========================================================================================================================*/
/**/
	/* Calculate density of whole target body */
/**/
/**/
	/* Calculate mean density of the planetary body BEFORE impact [kg/m^3] */
	stack         = 0.0000;
	stack_counter = 0;
/**/
	/* First the mean density of the iron core [kg/m^3] */
	for (m1=0;m1<xnumx;m1++)
	for (m2=0;m2<ynumy;m2++)
		{
		/* Node number */
		m3=m1*ynumy+m2;
		/**/
		if(ro[m3]>=crit_density_fe)
			{
			stack+=ro[m3];
			stack_counter+=1;
			}
		}
	rho_co = stack/(double)(stack_counter);
/**/
	/* Now the mean density of the silicate mantle [kg/m^3] */
	stack         = 0.0000;
	stack_counter = 0;
	/**/
	for(m1=0;m1<xnumx;m1++)
	for(m2=0;m2<ynumy;m2++)
		{
		/* Node number */
		m3=m1*ynumy+m2;
		/**/
		if(ro[m3]<crit_density_si && ro[m3]>crit_density)
			{
			stack+=ro[m3];
			stack_counter+=1;
			}
		}
	rho_ma = stack/(double)(stack_counter);
/**/
	/* Now compute the mean density of the planetary embryo [kg/m^3] */
	rho_planet = (rho_ma*(pow((current_radius*xsize),3.000)-pow(core_radius,3.000))+rho_co*pow(core_radius,3.000))/pow((current_radius*xsize),3.000);
/**/
	/* Calculate escape velocity of the target body (already squared!) [m^2/s^2] [Senshu et al., JGR, 107, 5118 (2002)] */
	escape_vel = 2.000*(4.000/3.000)*pival*rho_planet*G*current_radius*xsize*current_radius*xsize;
/**/
	/* Calculate ratio between impact velocity and escape velocity of target body [non-dim.] */
	v_ratio    = impact_vel[no1]/sqrt(escape_vel);
/**/
/**/
/*===========================================================================================================================*/
/**/
	/* Calculate first estimation of mean density and the core and mantle masses of impactor body */
	/**/
	/* Set parameters */
	rho_Fe_impact   = 0.0000;
	rho_Si_impact   = 0.0000;
	rho_impact      = 0.0000;
	mass_imp_core   = 0.0000;
	mass_imp_mantle = 0.0000;
	/**/
	/* Initial density of iron and silicates at temperature of space [kg/m^3] */
	rho_Fe_impact = rho_Fe*(1.000-markbb[7]*(t_space-298.000));
	rho_Si_impact = rho_Si*(1.000-markbb[5]*(t_space-298.000));
	/**/
	/* Mean density of agglomerated impactor body at temperature of space (neglecting pressure effects) [kg/m^3] */
	rho_impact      = rho_Fe_impact*pow(fe_frac,3.000)+rho_Si_impact*(1.000-pow(fe_frac,3.000));
	/**/
	/* Calculate initial radius of agglomerated impactor body [m] */
	impact_radius   = pow(((3.000*impact_mass[no1])/(4.000*pival*rho_impact)),(1.0000/3.0000));
	/**/
	/* Compute mass of the iron fraction [kg] */
	mass_imp_core   = (4.000/3.000)*pival*rho_Fe_impact*pow((fe_frac*impact_radius),3.000);
	/**/
	/* Compute mass of silicate fraction [kg] */
	mass_imp_mantle = impact_mass[no1]-mass_imp_core;
/**/
/*===========================================================================================================================*/
/**/
	/* Calculate density of future isobaric core region on target body */
/**/
	/* Use first mean density estimation for impactor body to estimate radius of isobaric core [Senshu et al., JGR, 107, 5118 (2002)] */
	iso_core = pow((3.000),(1.0000/3.0000))*(impact_radius/xsize);
/**/
	/* Calculate the mean density of material in future isobaric core region */
	no_ic    = 0;
	rho_ic   = 0.000;
	rho_icc  = 0.000;
/**/
	/* Computation on node level */
	for (m1=0;m1<xnumx;m1++)
	for (m2=0;m2<ynumy;m2++)
		{
		/* Node Num */
		m3=m1*ynumy+m2;
		/* Node X,Y General coordinates calc */
		x=gx[m1]/xsize;
		y=gy[m2]/ysize;
/**/
		/* Isobaric core centered on former surface */
		x1=x-xcoordsurf;
		y1=(y-ycoordsurf)*ysize/xsize;
/**/
		x4=pow(x1*x1+y1*y1,0.500);
		if(x4<=iso_core)
			{
			y4=pow(x1*x1+y1*y1,0.500);
			if(y4>=0.000*ysize/xsize)
				{
				/* Relative coordinates calc */
				x4+=y4;if(!x4) x4=1.000;
				x=y4/x4;
				x4=ABSV(360.000)-ABSV(0.000); if(!x4) x4=1.000;
				y=(ival-ABSV(0.000))/x4;
				/**/
				/* Stack densities in future isobaric core region */
				rho_icc+=ro[m3];
				no_ic+=1;
				}
			}
		}
	/**/
	/* Compute mean density of future isobaric core region [kg/m^3] */
	rho_ic = rho_icc/(double)(no_ic);
	/**/
/*===========================================================================================================================*/
	/**/
	/* Calculate upper limit for pre-impact thermal history of impactor body [K] */
	/**/
	/* Define parameters */
	e_undiff     = 0.0000;
	form_time    = 0.0000;
	test_time    = 0.0000;
	t_form       = 0.0000;
	t_form2      = 0.0000;
	/**/
	/* Time increment for test time [Ma] */
	time_incr = 0.00001;
	/**/
	/* Potential energy release due to accretion of undifferentiated impactor body [J] [Schubert et al., in: Satellites (1986)] */
	e_undiff  = (-3.000*G*impact_mass[no1]*impact_mass[no1])/(5.000*impact_radius);
	/**/
	for(no_it=1;no_it<=25000;no_it=no_it+1)
		{
		/* Erase old test values */
		test_time = 0.0000;
		e_al      = 0.0000;
		e_fe      = 0.0000;
		e_u       = 0.0000;
		e_uu      = 0.0000;
		e_k       = 0.0000;
		e_th      = 0.0000;
		t_test    = 0.0000;
		/**/
		/* Set new test time [Ma] */
		test_time = no_it*time_incr;
		/**/
		/* Radiogenic heating of impactor body until core separation happens */
		/* Decay constants in [1/s] */
		/**/
		/* Total energy in whole impactor until core separation provided by short-lived isotopes [J] */
		/* 26Al: t_1/2 = 0.716 Ma [Barr & Canup, Icarus, 198, 163-177 (2008)] */
        	e_al = (1.820e-7/3.0677e-14)*(exp(-3.0677e-14*start_time*1.000e+6*365.250*24.000*3600.000)-exp(-3.0677e-14*test_time*1.000e+6*365.250*24.000*3600.000))*impact_mass[no1];
		/**/
		/* New measurement: 60Fe: t_1/2 = 2.62 Ma [Rugel et al., Phys. Rev. Lett., 103, 072502 (2009)] */
		e_fe = (1.340e-9/8.3834e-15)*(exp(-8.3834e-15*start_time*1.000e+6*365.250*24.000*3600.000)-exp(-8.3834e-15*test_time*1.000e+6*365.250*24.000*3600.000))*impact_mass[no1];
		/**/
		/**/
		/* Total energy in whole impactor until core separation provided by long-lived isotopes [J] */
		/* 235U:  t_1/2 = 704 Ma [Barr & Canup, Icarus, 198, 163-177 (2008)] */
		e_u  = (2.990e-12/3.1200e-17)*(exp(-3.1200e-17*start_time*1.000e+6*365.250*24.000*3600.000)-exp(-3.1200e-17*test_time*1.000e+6*365.250*24.000*3600.000))*impact_mass[no1];
		/**/
		/* 238U:  t_1/2 = 4460 Ma [Barr & Canup, Icarus, 198, 163-177 (2008)] */
		e_uu = (1.600e-12/4.9248e-18)*(exp(-4.9248e-18*start_time*1.000e+6*365.250*24.000*3600.000)-exp(-4.9248e-18*test_time*1.000e+6*365.250*24.000*3600.000))*impact_mass[no1];
		/**/
		/*  40K:  t_1/2 = 1390 Ma [Barr & Canup, Icarus, 198, 163-177 (2008)] */
		/* Assume NO potassium depletion in impactor bodies */
		e_k  = (1.430e-11/1.5802e-17)*(exp(-1.5802e-17*start_time*1.000e+6*365.250*24.000*3600.000)-exp(-1.5802e-17*test_time*1.000e+6*365.250*24.000*3600.000))*impact_mass[no1];
		/**/
		/* 232Th: t_1/2 = 14000 Ma [Barr & Canup, Icarus, 198, 163-177 (2008)] */
		e_th = (1.000e-12/1.5689e-18)*(exp(-1.5689e-18*start_time*1.000e+6*365.250*24.000*3600.000)-exp(-1.5689e-18*test_time*1.000e+6*365.250*24.000*3600.000))*impact_mass[no1];
		/**/
		/**/
		/* Calculate expected temperature [K] */
		t_test = t_space+(((0.000-e_undiff)+e_al+e_fe+e_u+e_uu+e_k+e_th)/(c_p*impact_mass[no1]));
		/**/
		/* Iterate to consider effects of changing densities of iron and silicates and resulting change of radius of the body on expected temperature */
		for(no_ite=1;no_ite<=10;no_ite=no_ite+1)
			{
			/* Reset parameters */
			fe_fract           = 0.000;
			impact_core_radius = 0.000;
			impact_radius      = 0.000;
			rho_impact         = 0.000;
			rho_Fe_impact      = 0.000;
			rho_Si_impact      = 0.000;
			Si_melt_imp        = 0.000;
			/**/
			/* Compute new density of iron and silicates for new temperature [kg/m^3] */
			if(t_test<t_sol)                              /* Silicate mantle of impactor has sub-solidus temperature */
				{
				rho_Si_impact = rho_Si*(1.000-markbb[5]*(t_test-298.000));
				}
			else if((t_test>=t_sol) && (t_test<t_liq))    /* Silicate mantle of impactor is partially molten */
				{
				Si_melt_imp   = (t_test-t_sol)/(t_liq-t_sol);
				rho_Si_impact = rho_Si*(1.000-markbb[5]*(t_test-298.000))-Si_melt_imp*(rho_Si*(1.000-markbb[5]*(t_test-298.000))-rho_Si_melt*(1.000-markbb[25]*(t_test-298.000)));
				}
			else if(t_test>=t_liq)                        /* Silicate mantle of impactor is completely molten */
				{
				Si_melt_imp   = 1.000;
				rho_Si_impact = rho_Si_melt*(1.000-markbb[25]*(t_test-298.000));
				}
			/**/
			/* No density difference between solid and molten iron is assumed */
			rho_Fe_impact      = rho_Fe*(1.000-markbb[7]*(t_test-298.000));
			/**/
			/* Calculate new core radius (neglecting pressure effects) */
			impact_core_radius = pow((3.000*mass_imp_core/(4.000*pival*rho_Fe_impact)),(1.000/3.000));
			/**/
			/* Recalculate radius of differentiated impactor body [m] */
			impact_radius      = pow((3.000*mass_imp_mantle/(4.000*pival*rho_Si_impact)+pow(impact_core_radius,3.000)),(1.000/3.000));
			/**/
			/* Compute new linear fraction of iron in impactor body [non-dim.] */
			fe_fract           = impact_core_radius/impact_radius;
			/**/
			/* Compute new mean density of impactor body [kg/m^3] */
			rho_impact         = rho_Fe_impact*pow(fe_fract,3.000)+rho_Si_impact*(1.000-pow(fe_fract,3.000));
			/**/
			/* Reset parameters */
			e_undiff  = 0.0000;
			t_test    = 0.0000;
			/**/
			/* Recalculate potential energy release due to accretion and shrinking of undifferentiated impactor body [J] [Schubert et al., in: Satellites (1986)] */
			e_undiff = (-3.000*G*impact_mass[no1]*impact_mass[no1])/(5.000*impact_radius);
			/**/
			/* Recalculate expected temperature [K] */
			t_test = t_space+(((0.000-e_undiff)+e_al+e_fe+e_u+e_uu+e_k+e_th)/(c_p*impact_mass[no1]));
			/**/
			/* End of no_ite iteration loop */
			}
		/**/
		/* Save time as possible core formation time [Ma] as long as temperature is beneath or equal to the melting temperature of peridotites */
		if(t_test<=t_sol)
			{
			/* Possible core formation time [Ma] */
			form_time    = test_time;
			/**/
			/* Iron-silicate separation temperature [K] */
			t_form       = t_test;
			/**/
			/* Only the excess temperature [K] */
			t_form2      = t_test-t_space;
			}
		/**/
		if(t_test>t_sol)
			{
			/* Leave the iteration loop to save computational time */
			goto leave_loop;
			}
	/* End of no_it iteration loop */	
	}
	/**/
	/* Continuing here after obtaining minimum core formation time in impactor body */
	leave_loop: printf("Minimum core formation time in impactor body computed! \n");
	/**/
	/**/
	/* Now consider differentiation of impactor body when silicate melting temperature is reached */
	/**/
	/* Iterate to consider effects of changing densities of iron and silicates and resulting change of radius of the body on expected temperature */
	for(no_ite=1;no_ite<=10;no_ite=no_ite+1)
		{
		/* Reset parameters */
		e_diff                = 0.000;
		e_undiff              = 0.000;
		t_impactor_core_old   = 0.000;
		t_impactor_mantle_old = 0.000;
		/**/
		/* Calculate potential energy release due to accretion and change of radius of impactor body [J] [Schubert et al., in: Satellites (1986)] */
		e_undiff = (-3.000*G*impact_mass[no1]*impact_mass[no1])/(5.000*impact_radius);
		/**/
		/* Calculate potential energy release due to iron-silicate differentiation of the impactor body [J] [Schubert et al., in: Satellites (1986)] using the latest iron and silicate densities */
		e_diff=((-16.000/15.000)*pival*pival*G*pow(impact_radius,5.000))*(rho_Si_impact*rho_Si_impact+(5.000/2.000)*rho_Si_impact*(rho_impact-rho_Si_impact)+((3.000/2.000)*rho_Si_impact-rho_Fe_impact)*(rho_Si_impact-rho_Fe_impact)*pow(((rho_impact-rho_Si_impact)/(rho_Fe_impact-rho_Si_impact)),(5.000/3.000)));
		/**/
		/* Now compute maximum temperature increase of impactor core and mantle [K] after core formation assuming heat partioning into iron core */
		t_impactor_core_old   = t_form2+((e_undiff-e_diff)/(c_p*mass_imp_core));
		t_impactor_mantle_old = t_form2;
		/**/
		/* Reset parameters */
		fe_fract           = 0.000;
		impact_core_radius = 0.000;
		impact_radius      = 0.000;
		rho_impact         = 0.000;
		rho_Si_impact      = 0.000;
		rho_Fe_impact      = 0.000;
		Si_melt_imp        = 0.000;
		/**/
		/* Use first temperature estimations to compute corrected silicate and iron densities in impactor body [kg/m^3] (ignoring pressure-dependency of density) */
		if((t_impactor_mantle_old+t_space)<t_sol)                                         /* Silicate mantle of impactor has sub-solidus temperature */
			{
			rho_Si_impact = rho_Si*(1.000-markbb[5]*(t_impactor_mantle_old+t_space-298.000));
			}
		else if(((t_impactor_mantle_old+t_space)>=t_sol) && ((t_impactor_mantle_old+t_space)<t_liq))    /* Silicate mantle of impactor is partially molten */
			{
			Si_melt_imp   = (t_impactor_mantle_old+t_space-t_sol)/(t_liq-t_sol);
			rho_Si_impact = rho_Si*(1.000-markbb[5]*(t_impactor_mantle_old+t_space-298.000))-Si_melt_imp*(rho_Si*(1.000-markbb[5]*(t_impactor_mantle_old+t_space-298.000))-rho_Si_melt*(1.000-markbb[25]*(t_impactor_mantle_old+t_space-298.000)));
			}
		else if((t_impactor_mantle_old+t_space)>=t_liq)                                   /* Silicate mantle of impactor is completely molten */
			{
			Si_melt_imp   = 1.000;
			rho_Si_impact = rho_Si_melt*(1.000-markbb[25]*(t_impactor_mantle_old+t_space-298.000));
			}
		/**/
		/* No density difference between solid and molten iron is assumed */
		rho_Fe_impact         = rho_Fe*(1.000-markbb[7]*(t_impactor_core_old+t_space-298.000));
		/**/
		/* Calculate new core radius (neglecting pressure effects) */
		impact_core_radius    = pow((3.000*mass_imp_core/(4.000*pival*rho_Fe_impact)),(1.000/3.000));
		/**/
		/* Recalculate radius of differentiated impactor body [m] */
		impact_radius         = pow((3.000*mass_imp_mantle/(4.000*pival*rho_Si_impact)+pow(impact_core_radius,3.000)),(1.000/3.000));
		/**/
		/* Compute new linear fraction of iron in impactor body [non-dim.] */
		fe_fract              = impact_core_radius/impact_radius;
		/**/
		/* Compute new mean density of impactor body [kg/m^3] */
		rho_impact            = rho_Fe_impact*pow(fe_fract,3.000)+rho_Si_impact*(1.000-pow(fe_fract,3.000));
		/**/
		/* End of no_ite loop */
		}
	/**/
	printf("Post-differentiation temperatures in impactor body computed! \n");
	/**/
	/**/
	/* Now consider thermal evolution of differentiated impactor body until impact occurs */
	/* Reset parameters */
	e_al   = 0.000;
	e_fe   = 0.000;
	e_u    = 0.000;
	e_uu   = 0.000;
	e_k    = 0.000;
	e_th   = 0.000;
	t0     = 0.000;
	/**/
	/* Radiogenic heating of impactor body from core separation until impact happens */
	/* Decay constants in [1/s] */
	/**/
	/* Total energy in whole impactor until impact event provided by short-lived isotopes [J] */
	/* 26Al: t_1/2 = 0.716 Ma [Barr & Canup, Icarus, 198, 163-177 (2008)] */
        e_al = (1.820e-7/3.0677e-14)*(exp(-3.0677e-14*form_time*1.000e+6*365.250*24.000*3600.000)-exp(-3.0677e-14*impact_time[no1]*1.000e+6*365.250*24.000*3600.000))*impact_mass[no1];
	/**/
	/* New measurement: 60Fe: t_1/2 = 2.62 Ma [Rugel et al., Phys. Rev. Lett., 103, 072502 (2009)] */
	e_fe = (1.340e-9/8.3834e-15)*(exp(-8.3834e-15*form_time*1.000e+6*365.250*24.000*3600.000)-exp(-8.3834e-15*impact_time[no1]*1.000e+6*365.250*24.000*3600.000))*impact_mass[no1];
	/**/
	/**/
	/* Total energy in whole impactor until impact event provided by long-lived isotopes [J] */
	/* 235U:  t_1/2 = 704 Ma [Barr & Canup, Icarus, 198, 163-177 (2008)] */
	e_u  = (2.990e-12/3.1200e-17)*(exp(-3.1200e-17*form_time*1.000e+6*365.250*24.000*3600.000)-exp(-3.1200e-17*impact_time[no1]*1.000e+6*365.250*24.000*3600.000))*impact_mass[no1];
	/**/
	/* 238U:  t_1/2 = 4460 Ma [Barr & Canup, Icarus, 198, 163-177 (2008)] */
	e_uu = (1.600e-12/4.9248e-18)*(exp(-4.9248e-18*form_time*1.000e+6*365.250*24.000*3600.000)-exp(-4.9248e-18*impact_time[no1]*1.000e+6*365.250*24.000*3600.000))*impact_mass[no1];
	/**/
	/*  40K:  t_1/2 = 1390 Ma [Barr & Canup, Icarus, 198, 163-177 (2008)] */
	/* Assume NO potassium depletion in impactor bodies */
	e_k  = (1.430e-11/1.5802e-17)*(exp(-1.5802e-17*form_time*1.000e+6*365.250*24.000*3600.000)-exp(-1.5802e-17*impact_time[no1]*1.000e+6*365.250*24.000*3600.000))*impact_mass[no1];
	/**/
	/* 232Th: t_1/2 = 14000 Ma [Barr & Canup, Icarus, 198, 163-177 (2008)] */
	e_th = (1.000e-12/1.5689e-18)*(exp(-1.5689e-18*form_time*1.000e+6*365.250*24.000*3600.000)-exp(-1.5689e-18*impact_time[no1]*1.000e+6*365.250*24.000*3600.000))*impact_mass[no1];
	/**/
	/**/
	/* Calculate thermal anomaly due to the impact [K] */
	/* Takes into account that mean densities of impactor, target body and isobaric core can be different, impact velocity is not necessary equal to escape velocity */
	/* Hemispherical melt model is considered as in later implementation of heating [Tonks & Melosh, JGR, 98, 5319-5333 (1993)] */
	/* Model slightly overestimates mass of isobaric core region due to assumption of planar geometry */
	t0   = 8.000*pival*gamma_eff*pow(v_ratio,2.000)*rho_impact*rho_planet*G*pow((current_radius*xsize),2.000)/(9.000*2.700*rho_ic*c_p);
	/*t0   = (4.000*pival/9.000)*(gamma_eff/2.700)*((G*rho_impact*(current_radius)*(current_radius)*xsize*xsize)/c_p); */
	/**/
	/* Reset parameters */
	t_impactor_core   = 0.000;
	t_impactor_mantle = 0.000;
	/**/
	/* Now compute maximum temperature increase of impactor core and mantle [K] after collisional heating taking heat partitioning into iron phase into account */
	t_impactor_core   = t_impactor_core_old+t0+((e_fe+e_k)/(c_p*mass_imp_core));
	t_impactor_mantle = t_impactor_mantle_old+t0+((e_al+e_u+e_uu+e_th)/(c_p*mass_imp_mantle));
	/**/
	/* Iterate to consider effects of changing density on release of potential energy and temperature */
	for(no_ite=1;no_ite<=10;no_ite=no_ite+1)
		{
		/* Reset parameters */
		fe_fract           = 0.000;
		impact_core_radius = 0.000;
		impact_radius      = 0.000;
		rho_impact         = 0.000;
		rho_Si_impact      = 0.000;
		rho_Fe_impact      = 0.000;
		Si_melt_imp        = 0.000;
		/**/
		/* Use first temperature estimations to compute corrected silicate and iron densities in impactor body [kg/m^3] (ignoring pressure-dependency of density) */
		if((t_impactor_mantle+t_space)<t_sol)                                         /* Silicate mantle of impactor has sub-solidus temperature */
			{
			rho_Si_impact = rho_Si*(1.000-markbb[5]*(t_impactor_mantle+t_space-298.000));
			}
		else if(((t_impactor_mantle+t_space)>=t_sol) && ((t_impactor_mantle+t_space)<t_liq))    /* Silicate mantle of impactor is partially molten */
			{
			Si_melt_imp   = (t_impactor_mantle+t_space-t_sol)/(t_liq-t_sol);
			rho_Si_impact = rho_Si*(1.000-markbb[5]*(t_impactor_mantle+t_space-298.000))-Si_melt_imp*(rho_Si*(1.000-markbb[5]*(t_impactor_mantle+t_space-298.000))-rho_Si_melt*(1.000-markbb[25]*(t_impactor_mantle+t_space-298.000)));
			}
		else if((t_impactor_mantle+t_space)>=t_liq)                                   /* Silicate mantle of impactor is completely molten */
			{
			Si_melt_imp   = 1.000;
			rho_Si_impact = rho_Si_melt*(1.000-markbb[25]*(t_impactor_mantle+t_space-298.000));
			}
		/**/
		/* No density difference between solid and molten iron is assumed */
		rho_Fe_impact      = rho_Fe*(1.000-markbb[7]*(t_impactor_core+t_space-298.000));
		/**/
		/* Calculate new core radius (neglecting pressure effects) */
		impact_core_radius = pow((3.000*mass_imp_core/(4.000*pival*rho_Fe_impact)),(1.000/3.000));
		/**/
		/* Recalculate radius of differentiated impactor body [m] */
		impact_radius      = pow((3.000*mass_imp_mantle/(4.000*pival*rho_Si_impact)+pow(impact_core_radius,3.000)),(1.000/3.000));
		/**/
		/* Compute new linear iron fraction in impactor body [non-dim.] */
		fe_fract           = impact_core_radius/impact_radius;
		/**/
		/* Compute mean density of impactor body [kg/m^3] */
		rho_impact         = rho_Fe_impact*pow(fe_fract,3.000)+rho_Si_impact*(1.000-pow(fe_fract,3.000));
		/**/
		/* End of no_ite loop */
		}
	/**/
	printf("Temperatures in impactor body at collision time computed! \n");
/**/
/*===========================================================================================================================*/
/**/
	/* Calculate the accreted mass AFTER impact [kg] */
	M_acc         = M_acc+impact_mass[no1];
/**/
        /* Now compute the latest mass of the iron diapir material contributed by the impactor body [kg] and shift positions of old values */
        M_diapir_old = M_diapir_ol;
        M_diapir_ol  = M_diapir;
/**/
	/* Mass of newly added iron material added to the impactor mantle [kg] */
        M_diapir     = mass_imp_core;
/**/
/*===========================================================================================================================*/
/**/
	/* Non-dim. radius of isobaric core */
	isobaric_core = pow((3.000),(1.0000/3.0000))*(impact_radius/xsize);
/**/
	/* Calculate radius of transient crater [m] [Senshu et al., JGR, 107, 5118 (2002)] */
	crater_radius=1.160*pow((rho_impact/rho_planet),(1.0000/3.0000))*pow((impact_vel[no1]*impact_vel[no1]/escape_vel),0.220)*pow(current_radius*xsize,0.220)*pow(impact_radius,0.780);   /* dim. !!! [m] */
/**/
	crater_angle=(180.000*crater_radius)/(pival*current_radius*xsize);
/**/
/*===========================================================================================================================*/
/**/
	/* Iron ponds in the crater: Volume of crater remains unchanged */
/**/
	if(znumz==1)    /* Model is 2D cylindrical */
		{
		corr  = (1.000/2.000)*(crater_radius/xsize);    /* non-dim. */
		corr1 = corr;
		}
/**/
	if(znumz>=2)    /* Model is 3D spherical */
		{
		corr  = (1.000/2.000)*(crater_radius/xsize);   /* non-dim. */
		corr1 = corr;
		}
/**/
/**/
	/* The center of the smaller ejecta zone will be shifted by factor deflection */
	deflection = crater_radius/(2.000*xsize)-corr;   /* non-dim. */
/**/
	/* Correction factor in degrees */
	defl_corr  = (180.000*deflection)/(pival*current_radius);
/**/
/**/
	/* 2D/3D correction factor only when we use 2D setup */
	if((corr2d3d==0) && (znumz==1))
		{
		eros = 1.0000;       /* When problem is treated as 2D, no correction is applied */
		}
	else if((corr2d3d==1) && (znumz==1))
		{
		/* Calculate correction factor, so the radius of the planet in 2D cylindrical geometry fits the expected one in 3D spherical geometry */
                eros1 = pow((current_radius*current_radius*current_radius*xsize*xsize*xsize+impact_radius*impact_radius*impact_radius),(2.0000/3.0000)); /* Square of expected new radius AFTER impact in 3D spherical geometry [m^2] */
		eros2 = current_radius*current_radius*xsize*xsize;          /* Square of measured mean target radius BEFORE impact [m^2] */
		eros  = (eros1-eros2)/(impact_radius*impact_radius);        /* Expected correction factor e [non-dim.] */
		}
/**/
	/* Make sure no negative value is later used for Newton-Raphson iteration */
	/* Tests show that this correction gives exactly the same result except sign */
	if(eros<0.000)
		{
		eros = ABSV(eros);
		}
/**/
/**/
/*===========================================================================================================================*/
/**/
	/* Look for the coordinates of the centers of the ejecta circles */
/**/
	/* Predefine values for ejecta boundaries */
	xcoord_left       = 0.000000;
	ycoord_left       = 0.000000;
	xcoord_right      = 0.000000;
	ycoord_right      = 0.000000;
	planet_leftmost   = 0.000000;
	planet_rightmost  = 0.000000;
	testradius3       = 0.000000;
	testradius4       = 0.000000;
/**/
	/* Predefine values for crater boundaries */
	xcoord_lleft      = 0.000000;
	ycoord_lleft      = 0.000000;
	xcoord_rright     = 0.000000;
	ycoord_rright     = 0.000000;
	planet_lleftmost  = 0.000000;
	planet_rrightmost = 0.000000;
	testradius11      = 0.000000;
	testradius22      = 0.000000;
/**/
	/* Boundaries of the crater */
	impact1 = impact_angle[no1]+crater_angle-2.000*defl_corr;  /* Right corner of right ejecta zone. takes, if needed, into account that iron will be later put into the mantle */
	impact2 = impact_angle[no1]-crater_angle+2.000*defl_corr;  /* Left corner of left ejecta zone. takes, if needed, into account that iron will be later put into the mantle */
	status  = 0;
/**/
	/* Ejecta zones are located around the impact site */
	impact3 = impact_angle[no1]+(crater_angle/2.000)-defl_corr;  /* For right ejecta circle center */
	impact4 = impact_angle[no1]-(crater_angle/2.000)+defl_corr;  /* For left ejecta circle center */
/**/
	/* Make sure impact angle is always in the range 0 to 360 degrees */
	if((impact_angle[no1]+crater_angle-2.000*defl_corr)>=360.000)
		{
		impact1 = (impact_angle[no1]+crater_angle-2.000*defl_corr)-360.000;
		status  = 1;
		}
	if((impact_angle[no1]-crater_angle+2.000*defl_corr)<0.000)
		{
		impact2 = 360.000+(impact_angle[no1]-crater_angle+2.000*defl_corr);  /* Subtract negative angle from 360 degrees */
		status  = 1;
		}
/**/
	/* Make sure impact angle is always in the range 0 to 360 degrees */
	if((impact_angle[no1]+(crater_angle/2.000)-defl_corr)>=360.000)
		{
		impact3 = (impact_angle[no1]+(crater_angle/2.000)-defl_corr)-360.000;
		}
	if((impact_angle[no1]-(crater_angle/2.000)+defl_corr)<0.000)
		{
		impact4 = 360.000+(impact_angle[no1]-(crater_angle/2.000)+defl_corr);  /* Subtract negative angle from 360 degrees */
		}
/**/
	/* Go through all markers to find the centres of the ejecta zones */
/**/
	for (mm1=0;mm1<marknum;mm1++)
		{
		/* Find center of left ejecta zone */
/**/
		/* Check only solid molten/solid silicates and iron, stabilizing material and sticky air is excluded */
        	if(markt[mm1]==5 || markt[mm1]==6 || markt[mm1]==7 || markt[mm1]==8 || markt[mm1]==9 || markt[mm1]==10 || markt[mm1]==17 || markt[mm1]==18 || markt[mm1]==19 || markt[mm1]==25 || markt[mm1]==26)
			{
/**/
			/* Left corner at 0 < alpha < 90 degrees */
			if((impact4>=0.000 && impact4<90.000) && (markx[mm1]/xsize>=0.500) && (marky[mm1]/ysize<0.500))
				{
				angle4 = atan((markx[mm1]/xsize-0.500)/(0.500-marky[mm1]/ysize))*180.000/pival;
				if(ABSV(angle4-impact4)<=range)
					{
					testradius3=pow(((markx[mm1]/xsize-0.500)*(markx[mm1]/xsize-0.500)+(marky[mm1]/ysize-0.500)*(marky[mm1]/ysize-0.500)),0.500);
					planet_leftmost=MAXV(testradius3,planet_leftmost);
					if(testradius3==planet_leftmost)
						{
						/* Global coordinates of left boundary of impact crater */
						xcoord_left = (double)(markx[mm1])/xsize;
						ycoord_left = (double)(marky[mm1])/ysize;
						}
					}
				}
			/* Left corner at 90 < alpha < 180 degrees */
			else if((impact4>=90.000 && impact4<180.000) && (markx[mm1]/xsize>0.500) && (marky[mm1]/ysize>=0.500))
				{
				angle4 = 90.000+atan((marky[mm1]/ysize-0.500)/(markx[mm1]/xsize-0.500))*180.000/pival;
				if (ABSV(angle4-impact4)<=range)
					{
					testradius3=pow(((markx[mm1]/xsize-0.500)*(markx[mm1]/xsize-0.500)+(marky[mm1]/ysize-0.500)*(marky[mm1]/ysize-0.500)),0.500);
	        			planet_leftmost=MAXV(testradius3,planet_leftmost);
					if(testradius3==planet_leftmost)
						{
						/* Global coordinates of left boundary of impact crater */
                           			xcoord_left = (double)(markx[mm1])/xsize;
                           			ycoord_left = (double)(marky[mm1])/ysize;
						}
					}
				}
			/* Left corner at 180 < alpha < 270 degrees */
                	else if((impact4>=180.000 && impact4<270.000) && (markx[mm1]/xsize<=0.500) && (marky[mm1]/ysize>0.500))
				{
                        	angle4 = 180.000+atan((0.500-markx[mm1]/xsize)/(marky[mm1]/ysize-0.500))*180.000/pival;
                        	if (ABSV(angle4-impact4)<=range)
					{
					testradius3=pow(((markx[mm1]/xsize-0.500)*(markx[mm1]/xsize-0.500)+(marky[mm1]/ysize-0.500)*(marky[mm1]/ysize-0.500)),0.500);
                       			planet_leftmost=MAXV(testradius3,planet_leftmost);
					if(testradius3==planet_leftmost)
						{
						/* Global coordinates of left boundary of impact crater */
                           			xcoord_left = (double)(markx[mm1])/xsize;
                           			ycoord_left = (double)(marky[mm1])/ysize;
						}
					}
				}
			/* Left corner at 270 < alpha < 360 degrees */
                	else if((impact4>=270.000 && impact4<360.000) && (markx[mm1]/xsize<0.500) && (marky[mm1]/ysize<=0.500))
				{
                        	angle4 = 270.000+atan((0.500-marky[mm1]/ysize)/(0.500-markx[mm1]/xsize))*180.000/pival;
                        	if (ABSV(angle4-impact4)<=range)
					{
					testradius3=pow(((markx[mm1]/xsize-0.500)*(markx[mm1]/xsize-0.500)+(marky[mm1]/ysize-0.500)*(marky[mm1]/ysize-0.500)),0.500);
                       			planet_leftmost=MAXV(testradius3,planet_leftmost);
					if(testradius3==planet_leftmost)
						{
						/* Global coordinates of left boundary of impact crater */
                           			xcoord_left = (double)(markx[mm1])/xsize;
                           			ycoord_left = (double)(marky[mm1])/ysize;
						}
					}
				}
/**/
/*---------------------------------------------------------------------------------------------------------------------------*/
/**/
			/* Find center of right ejecta zone */
/**/
			/* Right corner at 0 < alpha < 90 degrees */
                	if((impact3>=0.000 && impact3<90.000) && (markx[mm1]/xsize>=0.500) && (marky[mm1]/ysize<0.500))
				{
                        	angle5 = atan((markx[mm1]/xsize-0.500)/(0.500-marky[mm1]/ysize))*180.000/pival;
                        	if(ABSV(angle5-impact3)<=range)
					{
					testradius4=pow(((markx[mm1]/xsize-0.500)*(markx[mm1]/xsize-0.500)+(marky[mm1]/ysize-0.500)*(marky[mm1]/ysize-0.500)),0.500);
                       			planet_rightmost=MAXV(testradius4,planet_rightmost);
					if(testradius4==planet_rightmost)
						{
						/* Global coordinates of right boundary of impact crater */
                               			xcoord_right = (double)(markx[mm1])/xsize;
                              			ycoord_right = (double)(marky[mm1])/ysize;
						}
					}
				}
			/* Right corner at 90 < alpha < 180 degrees */
                	else if((impact3>=90.000 && impact3<180.000) && (markx[mm1]/xsize>0.500) && (marky[mm1]/ysize>=0.500))
				{
                        	angle5 = 90.000+atan((marky[mm1]/ysize-0.500)/(markx[mm1]/xsize-0.500))*180.000/pival;
                        	if(ABSV(angle5-impact3)<=range)
					{
					testradius4=pow(((markx[mm1]/xsize-0.500)*(markx[mm1]/xsize-0.500)+(marky[mm1]/ysize-0.500)*(marky[mm1]/ysize-0.500)),0.500);
                       			planet_rightmost=MAXV(testradius4,planet_rightmost);
					if(testradius4==planet_rightmost)
						{
						/* Global coordinates of right boundary of impact crater */
                               			xcoord_right = (double)(markx[mm1])/xsize;
                               			ycoord_right = (double)(marky[mm1])/ysize;
						}
					}
				}
			/* Right corner at 180 < alpha < 270 degrees */
	                else if((impact3>=180.000 && impact3<270.000) && (markx[mm1]/xsize<=0.500) && (marky[mm1]/ysize>0.500))
				{
                	       	angle5 = 180.000+atan((0.500-markx[mm1]/xsize)/(marky[mm1]/ysize-0.500))*180.000/pival;
                        	if(ABSV(angle5-impact3)<=range)
					{
					testradius4=pow(((markx[mm1]/xsize-0.500)*(markx[mm1]/xsize-0.500)+(marky[mm1]/ysize-0.500)*(marky[mm1]/ysize-0.500)),0.500);
                       			planet_rightmost=MAXV(testradius4,planet_rightmost);
					if(testradius4==planet_rightmost)
						{
						/* Global coordinates of right boundary of impact crater */
                               			xcoord_right = (double)(markx[mm1])/xsize;
                               			ycoord_right = (double)(marky[mm1])/ysize;
						}
					}
				}
			/* Right corner at 270 < alpha < 360 degrees */
                	else if((impact3>=270.000 && impact3<360.000) && (markx[mm1]/xsize<0.500) && (marky[mm1]/ysize<=0.500))
				{
                        	angle5 = 270.000+atan((0.500-marky[mm1]/ysize)/(0.500-markx[mm1]/xsize))*180.000/pival;
                        	if(ABSV(angle5-impact3)<=range)
					{
					testradius4=pow(((markx[mm1]/xsize-0.500)*(markx[mm1]/xsize-0.500)+(marky[mm1]/ysize-0.500)*(marky[mm1]/ysize-0.500)),0.500);
                       			planet_rightmost=MAXV(testradius4,planet_rightmost);
					if(testradius4==planet_rightmost)
						{
						/* Global coordinates of right boundary of impact crater */
                               			xcoord_right = (double)(markx[mm1])/xsize;
                               			ycoord_right = (double)(marky[mm1])/ysize;
						}
					}
				}
/**/
			/* End of material loop */
			}
		/* End of marker loop */
		}
/**/
/**/
/*===========================================================================================================================*/
/**/
	/* Calculate the mean excess temperature in both future ejecta zones before the impact takes place */
	no_left = 0;
	t_left  = 0.000;
	t_leftt = 0.000;
/**/
	/* Done on the marker level! */
	for(mm1=0;mm1<marknum;mm1++)
		{
		x=(double)(markx[mm1])/xsize;  /* non-dim. */
		y=(double)(marky[mm1])/ysize;  /* non-dim. */
		x1=x-xcoord_left;
		y1=(y-ycoord_left)*ysize/xsize;
/**/
		/* Upper boundary for the crater is the bottom of the impact crater */
		if((x1*x1+y1*y1)<=(crater_radius/(2.000*xsize)*crater_radius/(2.000*xsize)))
			{
			/* Lower boundary for the crater is the protoplanet's surface */
			if((x1*x1+y1*y1)>=(0.0000*0.0000*ysize*ysize/xsize/xsize))
				{
				/* Calculate angles involved */
				ival=pow(x1*x1+y1*y1,0.500); if(!ival) ival=1.000;
				ival=x1/ival; ival=asin(ABSV(ival))/pival*180.000;
				if(x1>=0 && y1>0) ival=180.000-ival;
				if(x1<0 && y1>0) ival+=180.000;
				if(x1<0 && y1<=0) ival=360.000-ival;
				/**/
				/* Works out fine, when impact crater is not around 360 degrees */
				if(ival>=ABSV(0.0000) && ival<=ABSV(360.000))
					{
					/**/
					/* Only silicate solids and silicate melts are taken into account, iron and sticky air are not taken into account */
					if(markt[mm1]==5 || markt[mm1]==6 || markt[mm1]==25 || markt[mm1]==26)
						{
						/* Takes into account that not all of the future transient crater has to be within the isobaric core */
						if((x1*x1+y1*y1)<=(isobaric_core*isobaric_core*ysize*ysize/xsize/xsize))
							{
							t00 = t0;
							}
						else if((x1*x1+y1*y1)>(isobaric_core*isobaric_core*ysize*ysize/xsize/xsize))
							{
							t00 = t0*pow((isobaric_core/pow((x1*x1+y1*y1),0.500)),4.40000);
							}
						/**/
						/* Count and measure temperature of silicate markers in future transient impact crater  */
						t_left+=(markk[mm1]+t00-t_space); /* Excess temperature of affected material after impact heating */
						no_left=no_left+1;                /* Count particles involved for mena calculation */
						/**/
						t_leftt+=markk[mm1];              /* Count total temperature for PRE-impact density calculation */
						}
					}
				}
			}
		}
	t_left   = t_left/(double)(no_left);    /* Calculate mean excess temperature in future left ejecta zone */
	t_leftt  = t_leftt/(double)(no_left);   /* Calculate mean total temperature in future left ejecta zone */
/**/
/**/
	no_right = 0;
	t_right  = 0.000;
	t_rightt = 0.000;
/**/
	for(mm1=0;mm1<marknum;mm1++)
		{
		x=(double)(markx[mm1])/xsize;  /* non-dim. */
		y=(double)(marky[mm1])/ysize;  /* non-dim. */
		x1=x-xcoord_right;
		y1=(y-ycoord_right)*ysize/xsize;
/**/
		/* Upper boundary for the crater is the bottom of the impact crater */
		if((x1*x1+y1*y1)<=(crater_radius/(2.000*xsize)*crater_radius/(2.000*xsize)))
			{
			/* Lower boundary for the crater is the protoplanet's surface */
			if((x1*x1+y1*y1)>=(0.0000*0.0000*ysize*ysize/xsize/xsize))
				{
				/* Calculate angles involved */
				ival=pow(x1*x1+y1*y1,0.500); if(!ival) ival=1.000;
				ival=x1/ival; ival=asin(ABSV(ival))/pival*180.000;
				if(x1>=0 && y1>0) ival=180.000-ival;
				if(x1<0 && y1>0) ival+=180.000;
				if(x1<0 && y1<=0) ival=360.000-ival;
				/**/
				/* Works out fine, when impact crater is not around 360 degrees */
				if(ival>=ABSV(0.0000) && ival<=ABSV(360.000))
					{
					/**/
					/* Only silicate solids and silicate melts are taken into account, iron and sticky air are not taken into account */
					if(markt[mm1]==5 || markt[mm1]==6 || markt[mm1]==25 || markt[mm1]==26)
						{
						/* Takes into account that not all of the future transient crater has to be within the isobaric core */
						if((x1*x1+y1*y1)<=(isobaric_core*isobaric_core*ysize*ysize/xsize/xsize))
							{
							t00 = t0;
							}
						else if((x1*x1+y1*y1)>(isobaric_core*isobaric_core*ysize*ysize/xsize/xsize))
							{
							t00 = t0*pow((isobaric_core/pow((x1*x1+y1*y1),0.500)),4.40000);
							}
						/**/
						/* Count and measure temperature of silicate markers in future transient impact crater */
						t_right+=(markk[mm1]+t00-t_space); /* Excess temperature of affected material after impact heating */
						no_right=no_right+1;               /* Count particles involved for mean calculation */
						/**/
						t_rightt+=markk[mm1];               /* Count total temperature for PRE-impact density calculation */
						}
					}
				}
			}
		}
	t_right    = t_right/(double)(no_right);    /* Calculate mean excess temperature of future right ejecta zone */
        t_rightt   = t_rightt/(double)(no_right);   /* Calculate mean total temperature of future right ejecta zone */
/**/
	/* Calculate (arithmetic) mean excess temperature of the whole future ejecta zones */
	t_ejecta   = (t_left+t_right)/2.000;
/**/
	/* Calculate (arithmetic) mean total temperature of the whole future ejecta zones */
	t_ejectaa  = (t_leftt+t_rightt)/2.000;
/**/
	/* Calculate expected density of ejecta material BEFORE impact assuming pure silicate composition and ignoring pressure effects [kg/m^3] */
	/* Take peridotite solidus and liquidus T [Herzberg et al., Geochem. Geophys. Geosyst., 1 (2000); Wade and Wood, EPSL, 236, 78-95 (2005)] assuming P = 0 GPa */
/**/
	Si_melt    = 0.000;     /* Define silicate melt fraction in ejecta zone */
/**/
	if(t_ejectaa<t_sol)                                                /* PRE-impact ejecta zone material has sub-solidus temperature */
		{
		rho_ejecta = rho_Si*(1.000-markbb[5]*(t_ejectaa-298.000));
		}
	else if((t_ejectaa>=t_sol) && (t_ejectaa<t_liq))                   /* PRE-impact ejecta zone material is partially molten */
		{
		Si_melt    = (t_ejectaa-t_sol)/(t_liq-t_sol);
		rho_ejecta = rho_Si*(1.000-markbb[5]*(t_ejectaa-298.000))-Si_melt*(rho_Si*(1.000-markbb[5]*(t_ejectaa-298.000))-rho_Si_melt*(1.000-markbb[25]*(t_ejectaa-298.000)));
		}
	else if(t_ejectaa>=t_liq)                                          /* PRE-impact ejecta zone material is completely molten */
		{
		rho_ejecta = rho_Si_melt*(1.000-markbb[25]*(t_ejectaa-298.000));
		}
/**/
	/* Calculate mass distribution between ejecta zone and silicate mantle of impactor */
	if(znumz>=1)  /* Model is 3D spherical or at least simulates it, calculate masses */
		{
		/* Mass of impactor mantle [kg] */
		mass_imp     = (4.000/3.000)*pival*rho_Si_impact*pow(impact_radius,3.000)*(1.000-pow(fe_fract,3.000));
/**/
		/* Mass of ejecta from target body, approximated as half torus [kg] [Bronstein, Taschenbuch der Mathematik (2001)] using PRE-impact density of affected material */
		mass_ejecta  = (1.000/4.000)*pival*pival*rho_ejecta*pow(crater_radius,3.000);
		}
/**/
        mass_tot    = mass_imp+mass_ejecta;    /* Total mass of ejecta material [kg] */
/**/
	/* Calculate weighted mean excess temperature (compared to space) of ejecta mixture [K] */
	t_blanket   = (mass_ejecta/mass_tot)*t_ejecta+(mass_imp/mass_tot)*t_impactor_mantle;
/**/
	/* Calculate expected density of silicate blanket AFTER thermal equilibration [kg/m^3] */
	/* Take peridotite solidus and liquidus T [Herzberg et al., Geochem. Geophys. Geosyst., 1 (2000); Wade and Wood, EPSL, 236, 78-95 (2005)] assuming P = 0 GPa */
/**/
	Si_melt_blanket = 0.000;     /* Define silicate melt fraction in blanket */
/**/
	if((memory_si*t_blanket+t_space)<t_sol)                                                     /* Blanket material has sub-solidus temperature */
		{
		rho_blanket = rho_Si*(1.000-markbb[5]*((memory_si*t_blanket+t_space)-298.000));
		}
	else if(((memory_si*t_blanket+t_space)>=t_sol) && ((memory_si*t_blanket+t_space)<t_liq)) /* Blanket material is partially molten */
		{
		Si_melt_blanket = ((memory_si*t_blanket+t_space)-t_sol)/(t_liq-t_sol);
		rho_blanket = rho_Si*(1.000-markbb[5]*((memory_si*t_blanket+t_space)-298.000))-Si_melt_blanket*(rho_Si*(1.000-markbb[5]*((memory_si*t_blanket+t_space)-298.000))-rho_Si_melt*(1.000-markbb[25]*((memory_si*t_blanket+t_space)-298.000)));
		}
	else if((memory_si*t_blanket+t_space)>=t_liq)                                               /* Blanket material is completely molten */
		{
		rho_blanket = rho_Si_melt*(1.000-markbb[25]*((memory_si*t_blanket+t_space)-298.000));
		}
/**/
	/* Calculate new radii AFTER the impact by adding the thermally equilibrated and (virtually) squeezed blanket, so it has the impactor's mean density [m] */
	impact_radius_corr = pow((rho_blanket/rho_impact),(1.0000/3.0000))*pow((rho_Si_impact/rho_blanket),(1.0000/3.0000))*impact_radius;
	crater_radius_corr = pow((rho_blanket/rho_impact),(1.0000/3.0000))*pow((rho_ejecta/rho_blanket),(1.0000/3.0000))*crater_radius;
/**/
/**/
/*===========================================================================================================================*/
	/* Add thermal anomaly of impact on node level */
/**/
	/* Impact structure using parametrization [Monteux et al., GRL, 34, L24201 (2007)] */
	for (m1=0;m1<xnumx;m1++)
	for (m2=0;m2<ynumy;m2++)
		{
		/* Node Num */
		m3=m1*ynumy+m2;
		/* Node X,Y General coordinates calc */
		x=gx[m1]/xsize;
		y=gy[m2]/ysize;
/**/
		/* Isobaric core centered on former surface */
		x1=x-xcoordsurf;
		y1=(y-ycoordsurf)*ysize/xsize;
/**/
		x4=pow(x1*x1+y1*y1,0.500);
		if(x4<=isobaric_core)
			{
			y4=pow(x1*x1+y1*y1,0.500);
			if(y4>=0.000*ysize/xsize)
				{
				/* Relative coordinates calc */
				x4+=y4;if(!x4) x4=1.000;
				x=y4/x4;
				x4=ABSV(360.000)-ABSV(0.000); if(!x4) x4=1.000;
				y=(ival-ABSV(0.000))/x4;
				/* Thermal anomaly */
				tk[m3]+=t0;
				}
			}
		}
/**/
	/* Decaying thermal anomaly around impact structure */
	for (m1=0;m1<xnumx;m1++)
	for (m2=0;m2<ynumy;m2++)
		{
		/* Node Num */
		m3=m1*ynumy+m2;
		/* Node X,Y General coordinates calc */
		x=gx[m1]/xsize;
		y=gy[m2]/ysize;
/**/
		/* Decaying thermal anomaly centered on former surface */
		x1=x-xcoordsurf;
		y1=(y-ycoordsurf)*ysize/xsize;
/**/
		x4=pow((x1*x1+y1*y1),0.500);
		if(x4>=isobaric_core)
			{
			/* Calculate and check distance from second center */
			y4=pow((x1*x1+y1*y1),0.500);
			if(y4>=isobaric_core)
				{
				/* Relative coordinates calc */
				x4+=y4;if(!x4) x4=1.000;
				x=y4/x4;
				x4=ABSV(360.000)-ABSV(0.000); if(!x4) x4=1.000;
				y=(ival-ABSV(0.000))/x4;
				/* Thermal anomaly */
				tk[m3]+=t0*pow((isobaric_core/y4),4.40000);
				}
			}
		}
/**/
/*===========================================================================================================================*/
	/* Calculate PRE-impact radius of target body at location at maximum distance from previous impact sites */
	/* These region should be largely unaffected by relaxation processes */
/**/
	/* Set initially the distance from the planetary center to the surface to zero */
	planet_far     = 0.000000;     /* planet_far is non-dim.! */
	testradius_far = 0.000000;
	impact_far     = 0.000000;
/**/
	/* Location most distant from previous impacts [degrees] */
	/**/
	/* For very first impact not relevant */
	if(no1==1)
		{
		impact_far = impact_angle[no1];
		}
	/**/
	/* For second impact at maximum distance from first impact site */
	if(no1==2)
		{
		impact_far = impact_angle[no1-1]-180.000;
		}
	/**/
	/* Third or later impact most distant from last two previous impacts [degrees] for later impacts */
	if(no1>=3)
		{
		impact_far = (ABSV(impact_angle[no1-1]+impact_angle[no1-2])/2.000)+180.000;
		}
	/**/
	/* Make sure value is smaller than 360 degrees */
	if(impact_far>=360.000)
		{
		impact_far = impact_far-360.000;
		}
/**/
	/* Go through all markers */
	for(mm1=0;mm1<marknum;mm1++)
		{
		/* Check only solid silicates and iron, stabilizing material and molten silicates, only sticky air is excluded */
        	if(markt[mm1]==5 || markt[mm1]==6 || markt[mm1]==7 || markt[mm1]==8 || markt[mm1]==9 || markt[mm1]==10 || markt[mm1]==17 || markt[mm1]==18 || markt[mm1]==19 || markt[mm1]==25 || markt[mm1]==26)
			{
			/* Far angle 0 <= alpha < 90 degrees */
				if((impact_far>=0.000 && impact_far<90.000) && (markx[mm1]/xsize>=0.500) && (marky[mm1]/ysize<0.500))
				{
				angle_far = atan((markx[mm1]/xsize-0.500)/(0.500-marky[mm1]/ysize))*180.000/pival;
                		if(ABSV(angle_far-impact_far)<=range)
					{
					testradius_far=pow(((markx[mm1]/xsize-0.500)*(markx[mm1]/xsize-0.500)+(marky[mm1]/ysize-0.500)*(marky[mm1]/ysize-0.500)),0.500);
               				planet_far=MAXV(testradius_far,planet_far);
					}
				}
/**/
			/* Far angle 90 <= alpha < 180 degrees */
               		else if((impact_far>=90.000 && impact_far<180.000) && ((markx[mm1]/xsize)>0.500) && ((marky[mm1]/ysize)>=0.500))
				{
                        	angle_far = 90.000+atan((marky[mm1]/ysize-0.500)/(markx[mm1]/xsize-0.500))*180.000/pival;
                        	if(ABSV(angle_far-impact_far)<=range)
					{
					testradius_far=pow(((markx[mm1]/xsize-0.500)*(markx[mm1]/xsize-0.500)+(marky[mm1]/ysize-0.500)*(marky[mm1]/ysize-0.500)),0.500);
               				planet_far=MAXV(testradius_far,planet_far);
					}
				}
/**/
			/* Far angle 180 <= alpha < 270 degrees */
                	else if((impact_far>=180.000 && impact_far<270.000) && (markx[mm1]/xsize<=0.500) && (marky[mm1]/ysize>0.500))
				{
                        	angle_far = 180.000+atan((0.500-markx[mm1]/xsize)/(marky[mm1]/ysize-0.500))*180.000/pival;
                        	if(ABSV(angle_far-impact_far)<=range)
					{
					testradius_far=pow(((markx[mm1]/xsize-0.500)*(markx[mm1]/xsize-0.500)+(marky[mm1]/ysize-0.500)*(marky[mm1]/ysize-0.500)),0.500);
                       			planet_far=MAXV(testradius_far,planet_far);
					}
				}
/**/
			/* Far angle 270 <= alpha < 360 degrees */
                       	else if((impact_far>=270.000 && impact_far<360.000) && (markx[mm1]/xsize<0.500) && (marky[mm1]/ysize<=0.500))
				{
                        	angle_far = 270.000+atan((0.500-marky[mm1]/ysize)/(0.500-markx[mm1]/xsize))*180.000/pival;
                        	if(ABSV(angle_far-impact_far)<=range)
					{
					testradius_far=pow(((markx[mm1]/xsize-0.500)*(markx[mm1]/xsize-0.500)+(marky[mm1]/ysize-0.500)*(marky[mm1]/ysize-0.500)),0.500);
                        		planet_far=MAXV(testradius_far,planet_far);
					}
				}
			/* End of silicate/iron if loop */
			}
			/* End of marker loop */
		}
/**/
/**/
/*===========================================================================================================================*/
	/* Checking for the current radii of the planetary body in different sectors */
/**/
/**/
	/* Set initial value for the angle */
	planet_angle   = 0.000000;
/**/
/**/
	/* Set counter */
	ii = 0;
/**/
	for(planet_angle=0.000000;planet_angle<=360.000;planet_angle=planet_angle+intervall)
		{
		/* Start to count */
		ii = ii+1;
/**/
		/* Set initially the distance from the planetary center to the surface to zero */
       		planet_current[ii] = 0.000000;   /* planet_current is non-dim.! */
		radius_neww        = 0.000000;   /* uncorrected POST-impact radius */
		radius_new         = 0.000000;   /* corrected POST-impact radius */
		testradius5        = 0.000000;
/**/
		/* Go through all markers */
		for(mm1=0;mm1<marknum;mm1++)
			{
/**/
			/* Check only solid silicates and iron, stabilizing material and molten silicates, only sticky air is excluded */
        		if(markt[mm1]==5 || markt[mm1]==6 || markt[mm1]==7 || markt[mm1]==8 || markt[mm1]==9 || markt[mm1]==10 || markt[mm1]==17 || markt[mm1]==18 || markt[mm1]==19 || markt[mm1]==25 || markt[mm1]==26)
				{
				/* Impact angle 0 <= alpha < 90 degrees */
                       		if((planet_angle>=0.000 && planet_angle<90.000) && (markx[mm1]/xsize>=0.500) && (marky[mm1]/ysize<0.500))
					{
                               		angle6 = atan((markx[mm1]/xsize-0.500)/(0.500-marky[mm1]/ysize))*180.000/pival;
                               		if(ABSV(angle6-planet_angle)<=range)
						{
						testradius5=pow(((markx[mm1]/xsize-0.500)*(markx[mm1]/xsize-0.500)+(marky[mm1]/ysize-0.500)*(marky[mm1]/ysize-0.500)),0.500);
                       				planet_current[ii]=MAXV(testradius5,planet_current[ii]);
						}
					}
				/* Impact angle 90 <= alpha < 180 degrees */
                       		else if((planet_angle>=90.000 && planet_angle<180.000) && ((markx[mm1]/xsize)>0.500) && ((marky[mm1]/ysize)>=0.500))
					{
                               		angle6 = 90.000+atan((marky[mm1]/ysize-0.500)/(markx[mm1]/xsize-0.500))*180.000/pival;
                               		if(ABSV(angle6-planet_angle)<=range)
						{
						testradius5=pow(((markx[mm1]/xsize-0.500)*(markx[mm1]/xsize-0.500)+(marky[mm1]/ysize-0.500)*(marky[mm1]/ysize-0.500)),0.500);
                       				planet_current[ii]=MAXV(testradius5,planet_current[ii]);
						}
					}
/**/
				/* Impact angle 180 <= alpha < 270 degrees */
                       		else if((planet_angle>=180.000 && planet_angle<270.000) && (markx[mm1]/xsize<=0.500) && (marky[mm1]/ysize>0.500))
					{
                               		angle6 = 180.000+atan((0.500-markx[mm1]/xsize)/(marky[mm1]/ysize-0.500))*180.000/pival;
                               		if(ABSV(angle6-planet_angle)<=range)
						{
						testradius5=pow(((markx[mm1]/xsize-0.500)*(markx[mm1]/xsize-0.500)+(marky[mm1]/ysize-0.500)*(marky[mm1]/ysize-0.500)),0.500);
                        			planet_current[ii]=MAXV(testradius5,planet_current[ii]);
						}
					}
				/* Impact angle 270 <= alpha < 360 degrees */
                        	else if((planet_angle>=270.000 && planet_angle<360.000) && (markx[mm1]/ysize<0.500) && (marky[mm1]/ysize<=0.500))
					{
                                	angle6 = 270.000+atan((0.500-marky[mm1]/ysize)/(0.500-markx[mm1]/xsize))*180.000/pival;
                                	if(ABSV(angle6-planet_angle)<=range)
						{
						testradius5=pow(((markx[mm1]/xsize-0.500)*(markx[mm1]/xsize-0.500)+(marky[mm1]/ysize-0.500)*(marky[mm1]/ysize-0.500)),0.500);
                        			planet_current[ii]=MAXV(testradius5,planet_current[ii]);
						}
					}
				/* End of silicate/iron if loop */
				}
			/* End of marker loop */
			}
/**/
/**/
/**/
		/* As first add ejecta ring all over the planetary body */
/**/
		/* First go into this loop when two values for planetary radius are available */
		if(ii>=2)
			{
			if(ii<361)
				{
				/* Take maximum radius in sector for radius calculation */
				planet_curr = MAXV(planet_current[ii],planet_current[ii-1]);
/**/
				/* Take minimum radius in sector for lower boundary */
				planet_low  = MINV(planet_current[ii],planet_current[ii-1]);
				}

			/* Special case: 360 degrees (0 degrees => i=1) */
			else if(ii==361)
				{
				/* Take maximum radius in sector for radius calculation */
				planet_curr = MAXV(planet_current[ii-360],planet_current[ii-1]);
/**/
				/* Take minimum radius in sector for lower boundary */
				planet_low  = MINV(planet_current[ii-360],planet_current[ii-1]);
				}
/**/
			/* Calculate external boundary of the ejecta layer in the current sector */
/**/
			if(znumz==1) /* Model is 2D cylindrical */
				{
				radius_neww=pow((eros*pow((impact_radius_corr/xsize),2.000)*(1.000-pow(fe_fract,2.000))+(1.000/4.000)*pow((crater_radius_corr/xsize),2.000)+pow(planet_curr,2.000)),0.500);
				}
/**/
			else if(znumz>=2) /* Model is 3D spherical */
				{
				radius_neww=pow(((pow((impact_radius/xsize),3.000)*(1.000-pow(fe_fract,3.000))+(3.000/16.000)*pival*pow(corr,3.000))+pow(planet_curr,3.000)),(1.0000/3.0000));
				}
/**/
			/* Take into account that blanket material has in reality NOT the mean density of the impactor */
            		if(ii<361)
                		{
                		radius_new = pow(((rho_impact/rho_blanket)*(pow(radius_neww,3.000)-pow(planet_current[ii],3.000))+pow(planet_current[ii],3.000)),(1.0000/3.0000));
                		}
            		else if(ii==361)
                		{
                		radius_new = pow(((rho_impact/rho_blanket)*(pow(radius_neww,3.000)-pow(planet_current[ii-360],3.000))+pow(planet_current[ii-360],3.000)),(1.0000/3.0000));                
                		}
/**/
/**/
/* //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// */
/**/
			/* Calculate coordinates for rounding circles on the left and right rim of the crater (needed due to numerical reasons) */
/**/
			/* Calculate ejecta thickness at left crater rim */
			if((ABSV(ii-impact2)<=range) && ((ii-impact2)<=0.000))
				{
				/* Ejecta thickness (non-dim.) */
				ejecta_left = ABSV(radius_new-planet_low);
				/* Angle of left centre point */
				impact_left = impact2-(180.000*ejecta_left/(pival*planet_low));
/**/
				/* Check angles before proceeding with circle calculation */
/**/
				/* Define status */
				status3 = 0;
				/* Check left angles */
				/* Make sure angle is always in the range 0 to 360 degrees */
				if(impact_left>=360.000)
					{
					impact_left = impact_left-360.000;
					status3 = 0;
					}
				if(impact_left<0.000)
					{
					impact_left = 360.000+impact_left;              /* subtract negative angle from 360 degrees */
					status3 = 1;
					}
				if((impact_left+90.000)>=360.000)
					{
					impact_lleft = (impact_left+90.000)-360.000;
					status3 = 1;
					}
/**/
				/* Find coordinates of left rounding center */
				for (mm1=0;mm1<marknum;mm1++)
					{
/**/
					/* Check only solid molten/solid silicates and iron, stabilizing material and sticky air is excluded */
        				if(markt[mm1]==5 || markt[mm1]==6 || markt[mm1]==7 || markt[mm1]==8 || markt[mm1]==9 || markt[mm1]==10 || markt[mm1]==17 || markt[mm1]==18 || markt[mm1]==19 || markt[mm1]==25 || markt[mm1]==26)
						{
						/* Find left center for rounding */
/**/
						/* Left corner at 0 < alpha < 90 degrees */
						if((impact_left>=0.000 && impact_left<90.000) && (markx[mm1]/xsize>=0.500) && (marky[mm1]/ysize<0.500))
							{
							angle22 = atan((markx[mm1]/xsize-0.500)/(0.500-marky[mm1]/ysize))*180.000/pival;
							if((ABSV(angle22-impact_left)<=range) && ((angle22-impact_left)<=0.000))
								{
								testradius22=pow(((markx[mm1]/xsize-0.500)*(markx[mm1]/xsize-0.500)+(marky[mm1]/ysize-0.500)*(marky[mm1]/ysize-0.500)),0.500);
								planet_lleftmost=MAXV(testradius22,planet_lleftmost);
								if(testradius22==planet_lleftmost)
									{
									/* Global coordinates of left rounding center */
									xcoord_lleft = (double)(markx[mm1])/xsize;
									ycoord_lleft = (double)(marky[mm1])/ysize;
									}
								}
							}
						/* Left corner at 90 < alpha < 180 degrees */
						else if((impact_left>=90.000 && impact_left<180.000) && (markx[mm1]/xsize>0.500) && (marky[mm1]/ysize>=0.500))
							{
							angle22 = 90.000+atan((marky[mm1]/ysize-0.500)/(markx[mm1]/xsize-0.500))*180.000/pival;
							if ((ABSV(angle22-impact_left)<=range) && ((angle22-impact_left)<=0.000))
								{
								testradius22=pow(((markx[mm1]/xsize-0.500)*(markx[mm1]/xsize-0.500)+(marky[mm1]/ysize-0.500)*(marky[mm1]/ysize-0.500)),0.500);
	        						planet_lleftmost=MAXV(testradius22,planet_lleftmost);
								if(testradius22==planet_lleftmost)
									{
									/* Global coordinates of left rounding center */
                           						xcoord_lleft = (double)(markx[mm1])/xsize;
                           						ycoord_lleft = (double)(marky[mm1])/ysize;
									}
								}
							}
						/* Left corner at 180 < alpha < 270 degrees */
                				else if((impact_left>=180.000 && impact_left<270.000) && (markx[mm1]/xsize<=0.500) && (marky[mm1]/ysize>0.500))
							{
                        				angle22 = 180.000+atan((0.500-markx[mm1]/xsize)/(marky[mm1]/ysize-0.500))*180.000/pival;
                        				if ((ABSV(angle22-impact_left)<=range) && ((angle22-impact_left)<=0.000))
								{
								testradius22=pow(((markx[mm1]/xsize-0.500)*(markx[mm1]/xsize-0.500)+(marky[mm1]/ysize-0.500)*(marky[mm1]/ysize-0.500)),0.500);
                       						planet_lleftmost=MAXV(testradius22,planet_lleftmost);
								if(testradius22==planet_lleftmost)
									{
									/* Global coordinates of left rounding center */
                           						xcoord_lleft = (double)(markx[mm1])/xsize;
                           						ycoord_lleft = (double)(marky[mm1])/ysize;
									}
								}
							}
						/* Left corner at 270 < alpha < 360 degrees */
                				else if((impact_left>=270.000 && impact_left<360.000) && (markx[mm1]/xsize<0.500) && (marky[mm1]/ysize<=0.500))
							{
                        				angle22 = 270.000+atan((0.500-marky[mm1]/ysize)/(0.500-markx[mm1]/xsize))*180.000/pival;
                        				if ((ABSV(angle22-impact_left)<=range) && ((angle22-impact_left)<=0.000))
								{
								testradius22=pow(((markx[mm1]/xsize-0.500)*(markx[mm1]/xsize-0.500)+(marky[mm1]/ysize-0.500)*(marky[mm1]/ysize-0.500)),0.500);
                       						planet_lleftmost=MAXV(testradius22,planet_lleftmost);
								if(testradius22==planet_lleftmost)
									{
									/* Global coordinates of left rounding center */
                           						xcoord_lleft = (double)(markx[mm1])/xsize;
                           						ycoord_lleft = (double)(marky[mm1])/ysize;
									}
								}
							}
						/* End of material loop */
						}
					/* End of marker loop */
					}
				/* End of left center loop */
				}
/**/
/*---------------------------------------------------------------------------------------------------------------------------*/
/**/
			/* Calculate ejecta thickness at right crater rim */
			if((ABSV(ii-impact1)<=range) && ((ii-impact1)>0.000))
				{
				/* Ejecta thickness (non-dim.) */
				ejecta_right = ABSV(radius_new-planet_low);
				/* Angle of right centre point */
				impact_right = impact1+(180.000*ejecta_right/(pival*planet_low));
/**/
				/* Check angles before proceeding with position calculation */
/**/
				/* Define status */
				status4 = 0;
				/* Check right angles */
				/* Make sure angle is always in the range 0 to 360 degrees */
				if(impact_right>=360.000)
					{
					impact_right = impact_right-360.000;
					status4 = 0;
					}
				if((impact_right<0.000) && (impact_right-90.000)<0.000)
					{
					impact_right = 360.000+impact_right;            /* subtract negative angle from 360 degrees */
					impact_rright = 360.000+(impact_right-90.000);  /* subtract negative angle from 360 degrees */
					status4 = 0;
					}
				if(impact_right>0.000 && (impact_right-90.000)<0.000)
					{
					impact_rright = 360.000+(impact_right-90.000);  /* subtract negative angle from 360 degrees */
					status4 = 1;
					}
/**/
				/* Find coordinates of left rounding center */
				for(mm1=0;mm1<marknum;mm1++)
					{
/**/
					/* Check only solid molten/solid silicates and iron, stabilizing material and sticky air is excluded */
        				if(markt[mm1]==5 || markt[mm1]==6 || markt[mm1]==7 || markt[mm1]==8 || markt[mm1]==9 || markt[mm1]==10 || markt[mm1]==17 || markt[mm1]==18 || markt[mm1]==19 || markt[mm1]==25 || markt[mm1]==26)
						{
						/* Find right center for rounding */
/**/
						/* Right corner at 0 < alpha < 90 degrees */
                				if((impact_right>=0.000 && impact_right<90.000) && (markx[mm1]/xsize>=0.500) && (marky[mm1]/ysize<0.500))
						{
                        			angle11 = atan((markx[mm1]/xsize-0.500)/(0.500-marky[mm1]/ysize))*180.000/pival;
                        			if((ABSV(angle11-impact_right)<=range) && ((ii-impact1)>0.000))
							{
							testradius11=pow(((markx[mm1]/xsize-0.500)*(markx[mm1]/xsize-0.500)+(marky[mm1]/ysize-0.500)*(marky[mm1]/ysize-0.500)),0.500);
                       					planet_rrightmost=MAXV(testradius11,planet_rrightmost);
							if(testradius11==planet_rrightmost)
								{
								/* Global coordinates of right rounding center */
                               					xcoord_rright = (double)(markx[mm1])/xsize;
                              					ycoord_rright = (double)(marky[mm1])/ysize;
								}
							}
						}
						/* Right corner at 90 < alpha < 180 degrees */
                				else if((impact_right>=90.000 && impact_right<180.000) && (markx[mm1]/xsize>0.500) && (marky[mm1]/ysize>=0.500))
							{
                        				angle11 = 90.000+atan((marky[mm1]/ysize-0.500)/(markx[mm1]/xsize-0.500))*180.000/pival;
                        				if((ABSV(angle11-impact_right)<=range) && ((ii-impact1)>0.000))
								{
								testradius11=pow(((markx[mm1]/xsize-0.500)*(markx[mm1]/xsize-0.500)+(marky[mm1]/ysize-0.500)*(marky[mm1]/ysize-0.500)),0.500);
                       						planet_rrightmost=MAXV(testradius11,planet_rrightmost);
								if(testradius11==planet_rrightmost)
									{
									/* Global coordinates of right rounding center */
                               						xcoord_rright = (double)(markx[mm1])/xsize;
                               						ycoord_rright = (double)(marky[mm1])/ysize;
									}
								}
							}
						/* Right corner at 180 < alpha < 270 degrees */
	                			else if((impact_right>=180.000 && impact_right<270.000) && (markx[mm1]/xsize<=0.500) && (marky[mm1]/ysize>0.500))
							{
                	       				angle11 = 180.000+atan((0.500-markx[mm1]/xsize)/(marky[mm1]/ysize-0.500))*180.000/pival;
                        				if((ABSV(angle11-impact_right)<=range) && ((ii-impact1)>0.000))
								{
								testradius11=pow(((markx[mm1]/xsize-0.500)*(markx[mm1]/xsize-0.500)+(marky[mm1]/ysize-0.500)*(marky[mm1]/ysize-0.500)),0.500);
								planet_rrightmost=MAXV(testradius11,planet_rrightmost);
								if(testradius11==planet_rrightmost)
									{
									/* Global coordinates of right rounding center */
                               						xcoord_rright = (double)(markx[mm1])/xsize;
                               						ycoord_rright = (double)(marky[mm1])/ysize;
									}
								}
							}
						/* Right corner at 270 < alpha < 360 degrees */
                				else if((impact_right>=270.000 && impact_right<360.000) && (markx[mm1]/xsize<0.500) && (marky[mm1]/ysize<=0.500))
							{
                        				angle11 = 270.000+atan((0.500-marky[mm1]/ysize)/(0.500-markx[mm1]/xsize))*180.000/pival;
                        				if((ABSV(angle11-impact_right)<=range) && ((ii-impact1)>0.000))
								{
								testradius11=pow(((markx[mm1]/xsize-0.500)*(markx[mm1]/xsize-0.500)+(marky[mm1]/ysize-0.500)*(marky[mm1]/ysize-0.500)),0.500);
                       						planet_rrightmost=MAXV(testradius11,planet_rrightmost);
								if(testradius11==planet_rrightmost)
									{
									/* Global coordinates of right rounding center */
                               						xcoord_rright = (double)(markx[mm1])/xsize;
                               						ycoord_rright = (double)(marky[mm1])/ysize;
									}
								}
							}
						/* End of material loop */
						}
					/* End of marker loop */
					}
				/* End of right center loop */
				}
/**/
/* ///////////////////////////////////////////////////////////////////////////////////////////////////////////////// */
/**/
			/* Now go through all markers to create ejecta blanket */
/**/
			for(mm1=0;mm1<marknum;mm1++)
				{
				x=(double)(markx[mm1])/xsize;
				y=(double)(marky[mm1])/ysize;
				x1=x-0.50000;
				y1=(y-0.50000)*ysize/xsize;
/**/
				/* Upper boundary for the the global ejecta blanket */
				if((x1*x1+y1*y1)<=((radius_new)*(radius_new)))
					{
					/* Lower boundary for the crater is the planetary surface minus an additional term to make sure that no sticky air gets entrained into the planetary mantle */
					if((x1*x1+y1*y1)>=((planet_low)*(planet_low)*ysize*ysize/xsize/xsize))
						{
						/* Calculate angles involved */
						ival=pow(x1*x1+y1*y1,0.500); if(!ival) ival=1.000;
						ival=x1/ival; ival=asin(ABSV(ival))/pival*180.000;
						if(x1>=0 && y1>0) ival=180.000-ival;
						if(x1<0 && y1>0) ival+=180.000;
						if(x1<0 && y1<=0) ival=360.000-ival;
						/**/
						/* Put the ring in small steps all over the planetary body */
						if(ival>=ABSV((ii-2)*intervall) && ival<=ABSV((ii-1)*intervall))
							{
							/* Sticky air markers become solid silicates */
							markt[mm1]=5;
							/**/
							/* Note when the material is weakened to make it again stiff after relax_time passes */
							markhi[mm1]=impact_time[no1];  /* in [Ma] */
							}
						}
					}
				}
/**/
			/* Temperature of silicate ejecta blanket (retains certain part of thermal history of impactor body and preimpact target mantle) */
			for (m1=0;m1<xnumx;m1++)
			for (m2=0;m2<ynumy;m2++)
				{
				/* Node Num */
				m3=m1*ynumy+m2;
				/* Node X,Y General coordinates calc */
				x=gx[m1]/xsize;
				y=gy[m2]/ysize;
/**/
				x1=x-0.50000;
				y1=(y-0.50000)*ysize/xsize;
/**/
				x4=pow((x1*x1+y1*y1),0.500);
				if(x4<=radius_new)
					{
					/* Calculate and check distance from second center */
					y4=pow((x1*x1+y1*y1),0.500);
					if(y4>=planet_low*ysize/xsize)
						{
						/* Calculate and check angle relatively first or second center */
						x3=x1;y3=y1; if((ii-1)*intervall<0) {x3=x2; y3=y2;} 
						ival=pow(x3*x3+y3*y3,0.500); if(!ival) ival=1.000;
						ival=x3/ival; ival=asin(ABSV(ival))/pival*180.000;
						if(x3>=0 && y3>0) ival=180.000-ival;
						if(x3<0 && y3>0) ival+=180.000;
						if(x3<0 && y3<=0) ival=360.000-ival;
						/* Put the thermal ring in small steps all over the planetary body */
						if(ival>=ABSV((ii-2)*intervall) && ival<=ABSV((ii-1)*intervall))
							{
							/* Relative coordinates calc */
							x4+=y4;if(!x4) x4=1.000;
							x=y4/x4;
							x4=ABSV((ii-1)*intervall)-ABSV((ii-2)*intervall); if(!x4) x4=1.000;
							/* Thermal anomaly */
							y=(ival-ABSV((ii-1)*intervall))/x4;
                                                        /* First reset temperature on specified node to zero */
                                                        tk[m3]=0.0000;
                                                        /* Now set temperature of new ejecta layer on top of target body to its value [K] */
							tk[m3]=(memory_si*t_blanket)+t_space;
							}
						}
					}
				}
			/* End of angle if loop */
			}
		/* End of angle for loop */
		}
/**/
	/* Go through all markers to remove the the ejecta blanket just above the crater again */
/**/
	for(mm1=0;mm1<marknum;mm1++)
		{
		x=(double)(markx[mm1])/xsize;
		y=(double)(marky[mm1])/ysize;
		x1=x-0.50000;
		y1=(y-0.50000)*ysize/xsize;
/**/
		/* Upper boundary for the crater is the protoplanet's surface plus a small additional factor to make sure that all silicates/iron previosly entrained in the sticky air or added before as ejecta material on top of the crater are removed! */
		if((x1*x1+y1*y1)<=((1.30*radius_new)*(1.30*radius_new)))
			{
			/* Lower boundary for the crater is the original surface before the impact */
			if((x1*x1+y1*y1)>=(planet_surf*planet_surf*ysize*ysize/xsize/xsize))
				{
				/* Calculate angles involved */
				ival=pow(x1*x1+y1*y1,0.500); if(!ival) ival=1.000;
				ival=x1/ival; ival=asin(ABSV(ival))/pival*180.000;
				if(x1>=0 && y1>0) ival=180.000-ival;
				if(x1<0 && y1>0) ival+=180.000;
				if(x1<0 && y1<=0) ival=360.000-ival;
				/**/
				/* Works out fine, when impact crater is not around 360 degrees */
				if(status == 0 && (ival>=ABSV(impact2) && ival<=ABSV(impact1)))
					{
					/* Markers in impact crater become sticky air */
					markt[mm1]=0;
					}
				/* Special case when impact crater is around 360 degrees */
				/* First fill the part 3xx to 360 degrees with markers */
				if(status == 1 && (ival>=ABSV(impact2) && ival<=360.000))
					{
					/* Markers in impact crater become sticky air */
					markt[mm1]=0;
					}
				/* Now fill the part 0 to xx degrees with markers */
				if(status == 1 && (ival>=0.000 && ival<=ABSV(impact1)))
					{
					/* Markers in impact crater become sticky air */
					markt[mm1]=0;
					}
				}
			}
		}
/**/
/* /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// */
/**/
/* radius_new for distant measurement region */
	radius_new_far = pow(((rho_impact/rho_blanket)*(pow(radius_neww,3.000)-pow(planet_far,3.000))+pow(planet_far,3.000)),(1.0000/3.0000));
/**/
/* Thickness of ejecta layer on top of excavation zone (PRE-impact radius planet_far measured at maximum distance from previous impacts, which should be largely unaffected by relaxation processes, used) */
	d_layer = ABSV(radius_new_far-planet_far);
/**/
/*===========================================================================================================================*/
/**/
	/* Add impact induced ejecta zones and fill them with sticky air material */
/**/
/**/
	/* Go through all markers to create the two ejecta zones */
/**/
	/* Left ejecta zone */
	for(mm1=0;mm1<marknum;mm1++)
		{
		x=(double)(markx[mm1])/xsize;  /* non-dim. */
		y=(double)(marky[mm1])/ysize;  /* non-dim. */
		x1=x-xcoord_left;
		y1=(y-ycoord_left)*ysize/xsize;
/**/
		/* Upper boundary for the crater is the protoplanet's surface */
		if((x1*x1+y1*y1)<=(corr*corr))
			{
			/* Lower boundary for the crater is the bottom of the impact crater */
			if((x1*x1+y1*y1)>=(0.0000*0.0000*ysize*ysize/xsize/xsize))
				{
				/* Calculate angles involved */
				ival=pow(x1*x1+y1*y1,0.500); if(!ival) ival=1.000;
				ival=x1/ival; ival=asin(ABSV(ival))/pival*180.000;
				if(x1>=0 && y1>0) ival=180.000-ival;
				if(x1<0 && y1>0) ival+=180.000;
				if(x1<0 && y1<=0) ival=360.000-ival;
				/**/
				if(ival>=ABSV(0.0000) && ival<=ABSV(360.000))
					{
					/* Markers in impact crater become sticky air */
					markt[mm1]=0;
					}
				}
			}
		}
/**/
	/* Right ejecta zone */
	for(mm1=0;mm1<marknum;mm1++)
		{
		x=(double)(markx[mm1])/xsize;  /* non-dim. */
		y=(double)(marky[mm1])/ysize;  /* non-dim. */
		x1=x-xcoord_right;
		y1=(y-ycoord_right)*ysize/xsize;
/**/
		/* Upper boundary for the crater is the protoplanet's surface */
		if((x1*x1+y1*y1)<=(corr*corr))
			{
			/* Lower boundary for the crater is the bottom of the impact crater */
			if((x1*x1+y1*y1)>=(0.0000*0.0000*ysize*ysize/xsize/xsize))
				{
				/* Calculate angles involved */
				ival=pow(x1*x1+y1*y1,0.500); if(!ival) ival=1.000;
				ival=x1/ival; ival=asin(ABSV(ival))/pival*180.000;
				if(x1>=0 && y1>0) ival=180.000-ival;
				if(x1<0 && y1>0) ival+=180.000;
				if(x1<0 && y1<=0) ival=360.000-ival;
				/**/
				if(ival>=ABSV(0.0000) && ival<=ABSV(360.000))
					{
					/* Markers in impact crater become sticky air */
					markt[mm1]=0;
					}
				}
			}
		}
/**/
	/* Round the central peak of the excavation zone (due to numerical reasons) */
/**/
	offsetcirc = crater_radius/(3.000*xsize);   /* Set the center of the new circle in the middle of the central "hill" */
/**/
	/* Predefine values */
	xcoordcirc = 0.000000;
	ycoordcirc = 0.000000;
	deltaxcirc = 0.000000;
	deltaycirc = 0.000000;
/**/
	/* Calculate coordinates of the circle center to round central peak */
/**/
	if(impact_angle[no1]>=0.000 && impact_angle[no1]<=90.000)
		{
		/* Coordinate correction */
		deltaxcirc = ABSV(sin(impact_angle[no1]/180.000*pival)*offsetcirc);
		deltaycirc = ABSV(cos(impact_angle[no1]/180.000*pival)*offsetcirc);
		/* Rounding circle coordinates */
		xcoordcirc = xcoordsurf-deltaxcirc;
		ycoordcirc = ycoordsurf+deltaycirc;
		}
/**/
	else if(impact_angle[no1]>90.000 && impact_angle[no1]<=180.000)
		{
		/* Coordinate correction */
		deltaxcirc = ABSV(cos((impact_angle[no1]-90.000)/180.000*pival)*offsetcirc);
		deltaycirc = ABSV(sin((impact_angle[no1]-90.000)/180.000*pival)*offsetcirc);
		/* Rounding circle coordinates */
		xcoordcirc = xcoordsurf-deltaxcirc;
		ycoordcirc = ycoordsurf-deltaycirc;
		}
/**/
	else if(impact_angle[no1]>180.000 && impact_angle[no1]<=270.000)
		{
		/* Coordinate correction */
		deltaxcirc = ABSV(sin((impact_angle[no1]-180.000)/180.000*pival)*offsetcirc);
		deltaycirc = ABSV(cos((impact_angle[no1]-180.000)/180.000*pival)*offsetcirc);
		/* Rounding circle coordinates */
		xcoordcirc = xcoordsurf+deltaxcirc;
		ycoordcirc = ycoordsurf-deltaycirc;
		}
/**/
	else if(impact_angle[no1]>270.000 && impact_angle[no1]<360.000)
		{
		/* Coordinate correction */
		deltaxcirc = ABSV(cos((impact_angle[no1]-270.000)/180.000*pival)*offsetcirc);
		deltaycirc = ABSV(sin((impact_angle[no1]-270.000)/180.000*pival)*offsetcirc);
		/* Rounding circle coordinates */
		xcoordcirc = xcoordsurf+deltaxcirc;
		ycoordcirc = ycoordsurf+deltaycirc;
		}
/**/
	/* Check angles involved */
	impact5 = impact_angle[no1]-45.000;
        impact6 = impact_angle[no1]+45.000;
	status2 = 0;
/**/
        if((impact_angle[no1]+45.000)>=360.000)
		{
		impact6 = (impact_angle[no1]+45.000)-360.000;  /* substract 360 degrees from the too large angle */
		status2 = 1;
		}
        if((impact_angle[no1]-45.000)<0.000)
		{
		impact5 = 360.000+(impact_angle[no1]-45.000);  /* subtract negative angle from 360 degrees */
		status2 = 1;
		}
/**/
	/* Go through all markers to round the central peak */
	for(mm1=0;mm1<marknum;mm1++)
		{
		/* Calculate and check distance from first center */
		x=(double)(markx[mm1])/xsize;
		y=(double)(marky[mm1])/ysize;
/**/
		x1=x-xcoordcirc;
		y1=(y-ycoordcirc)*ysize/xsize;
/**/
		/* Upper boundary for the rounding, use a little bit more than 1/3 to make sure that all silicate/iron markers are removed! */
		if((x1*x1+y1*y1)<=((crater_radius/xsize*0.360)*(crater_radius/xsize*0.360)))
			{
/**/
			/* Lower boundary for the rounding */
			if((x1*x1+y1*y1)>=((crater_radius/xsize*0.100*ysize/xsize)*(crater_radius/xsize*0.100*ysize/xsize)))
				{
				/* Calculate angles involved */
				ival=pow(x1*x1+y1*y1,0.500); if(!ival) ival=1.000;
				ival=x1/ival; ival=asin(ABSV(ival))/pival*180.000;
				if(x1>=0 && y1>0) ival=180.000-ival;
				if(x1<0 && y1>0) ival+=180.000;
				if(x1<0 && y1<=0) ival=360.000-ival;
				/**/
				/* Works out fine, when impact crater is not around 360 degrees */
				if(status2 == 0 && (ival>=ABSV(impact5) && ival<=ABSV(impact6)))
					{
					/* Markers in rounding area become sticky air */
					markt[mm1]=0;
					}
				/* Special case when impact crater is around 360 degrees */
				/* First fill the part 3xx to 360 degrees with markers */
				if(status2 == 1 && (ival>=ABSV(impact5) && ival<=360.000))
					{
					/* Markers in rounding area become sticky air */
					markt[mm1]=0;
					}
				/* Now fill the part 0 to xx degrees with markers */
				if(status2 == 1 && (ival>=0.000 && ival<=ABSV(impact6)))
					{
					/* Markers in rounding area become sticky air */
					markt[mm1]=0;
					}
				}
			}
		}
/**/
/*===========================================================================================================================*/
/**/
	/* Proceed with corner rounding */
/**/
	/* Go through all markers to add circles to round up the left and right corner of the ejecta zones after ejecta is already added on top */
/**/
	/* Left rounding circle */
	for(mm1=0;mm1<marknum;mm1++)
		{
		/* Calculate and check distance from first center */
		x=(double)(markx[mm1])/xsize;
		y=(double)(marky[mm1])/ysize;
/**/
		x1=x-xcoord_lleft;
		y1=(y-ycoord_lleft)*ysize/xsize;
/**/
		/* Upper boundary for the rounding, use a little bit more to make sure that all silicate/iron markers are removed! */
		if((x1*x1+y1*y1)<=((2.000*ejecta_left)*(2.000*ejecta_left)))
			{
			/* Lower boundary for the rounding */
			if((x1*x1+y1*y1)>=((1.000*ejecta_left*ysize/xsize)*(1.000*ejecta_left*ysize/xsize)))
				{
				/* Calculate angles involved */
				ival=pow(x1*x1+y1*y1,0.500); if(!ival) ival=1.000;
				ival=x1/ival; ival=asin(ABSV(ival))/pival*180.000;
				if(x1>=0 && y1>0) ival=180.000-ival;
				if(x1<0 && y1>0) ival+=180.000;
				if(x1<0 && y1<=0) ival=360.000-ival;
				/**/
				/* Works out fine, when impact crater is not around 360 degrees */
				if(status3 == 0 && (ival>=ABSV(impact_left) && ival<=ABSV(impact_left+90.000)))
					{
					/* Markers in rounding area become sticky air */
					markt[mm1]=0;
					}
				/* Special case when rounding area is around 360 degrees */
				/* First fill the part 3xx to 360 degrees with markers */
				if(status3 == 1 && (ival>=ABSV(impact_left) && ival<=360.000))
					{
					/* Markers in rounding area become sticky air */
					markt[mm1]=0;
					}
				/* Now fill the part 0 to xx degrees with markers */
				if(status3 == 1 && (ival>=0.000 && ival<=ABSV(impact_lleft)))
					{
					/* Markers in rounding area become sticky air */
					markt[mm1]=0;
					}
				}
			}
		}
/**/
	/* Right rounding circle */
	for(mm1=0;mm1<marknum;mm1++)
		{
		/* Calculate and check distance from first center */
		x=(double)(markx[mm1])/xsize;
		y=(double)(marky[mm1])/ysize;
/**/
		x1=x-xcoord_rright;
		y1=(y-ycoord_rright)*ysize/xsize;
/**/
/* Upper boundary for the rounding, use a little bit more to make sure that all silicate/iron markers are removed! */
		if((x1*x1+y1*y1)<=((2.000*ejecta_right)*(2.000*ejecta_right)))
			{
			/* Lower boundary for the rounding */
			if((x1*x1+y1*y1)>=((1.000*ejecta_right*ysize/xsize)*(1.000*ejecta_right*ysize/xsize)))
				{
				/* Calculate angles involved */
				ival=pow(x1*x1+y1*y1,0.500); if(!ival) ival=1.000;
				ival=x1/ival; ival=asin(ABSV(ival))/pival*180.000;
				if(x1>=0 && y1>0) ival=180.000-ival;
				if(x1<0 && y1>0) ival+=180.000;
				if(x1<0 && y1<=0) ival=360.000-ival;
				/**/
				/* Works out fine, when impact crater is not around 360 degrees */
				if(status4 == 0 && (ival>=ABSV(impact_right-90.000) && ival<=ABSV(impact_right)))
					{
					/* Markers in rounding area become sticky air */
					markt[mm1]=0;
					}
				/* Special case when rounding is around 360 degrees */
				/* First fill the part 3xx to 360 degrees with markers */
				if(status4 == 1 && (ival>=ABSV(impact_rright) && ival<=360.000))
					{
					/* Markers in rounding area become sticky air */
					markt[mm1]=0;
					}
				/* Now fill the part 0 to xx degrees with markers */
				if(status4 == 1 && (ival>=0.000 && ival<=ABSV(impact_right)))
					{
					/* Markers in rounding area become sticky air */
					markt[mm1]=0;
					}
				}
			}
		}
/**/
/*===========================================================================================================================*/
/**/
	/* Put the iron of the impactor core as ponds into the crater */
/**/
/**/
	/* Recalculate iron density after impact as temperature changes */
	rho_Fe_pond = rho_Fe*(1.000-markbb[7]*((memory_fe*t_impactor_core+t_space)-298.000));
/**/
	/* Recalculate size of iron core as material changes temperature when emplaced after impact */
	fe_fract    = fe_fract*pow((rho_Fe_impact/rho_Fe_pond),(1.0000/3.0000));
/**/
	if(fe_fract<=0.000)
		{
		fe_fract = 0.000;
		}
	if(fe_fract>=1.000)
		{
		fe_fract = 1.000;
		}
/**/
	/* Define minimum value for phi_angle */
	phi_angle = 0.000;
/**/
	/* Compute phi_angle only if iron was added to impactor mantle */
	if(fe_fract>0.000)
		{
		/* Calculate angle of iron pond by using Newton-Raphson method and Taylor expansion of the sine function to degree 6 */
		/**/
        	/* Initial values for angle phi and number of iterations */
		phi_angle   = 100.000;
		no_iter_max = 10;
		/**/
		/* Adapt parameters, if necessary */
		change_angle: phi_angle = phi_angle/2.000;   /* Halve the initial guess for the angle phi */
		reiterate: no_iter_max  = 2*no_iter_max;     /* Double the amount of iterations, if result is unsatisfactory */
/**/
		for(no_iter=1;no_iter<=no_iter_max;no_iter=no_iter+1)
			{
			f = 1.000/6.000*pow((pival*phi_angle/180.000),3.000)-1.000/120.000*pow((pival*phi_angle/180.000),5.000)+1.000/5040.000*pow((pival*phi_angle/180.000),7.000)-1.000/362880.000*pow((pival*phi_angle/180.000),9.000)+1.000/39916800.000*pow((pival*phi_angle/180.000),11.000)-((4.000*pival*eros*fe_fract*fe_fract*impact_radius*impact_radius/xsize/xsize)/(2.000*2.000*corr*corr));
/**/
			g = 1.000/2.000*pow((pival/180.000),3.000)*pow(phi_angle,2.000)-1.000/24.000*pow((pival/180.000),5.000)*pow(phi_angle,4.000)+1.000/720.000*pow((pival/180.000),7.000)*pow(phi_angle,6.000)-1.000/40320.000*pow((pival/180.000),9.000)*pow(phi_angle,8.000)+1.000/3628800.000*pow((pival/180.000),11.000)*pow(phi_angle,10.000);
/**/
			/* Save the old phi angle value */
			phi_angle_old = phi_angle;
/**/
			/* Calculate new phi value */
			phi_angle = phi_angle_old-(f/g);
			}
/**/
		/* Is the angle phi varying significantly? */
		phi_error = phi_angle-phi_angle_old;
/**/
		if(phi_error>=5.000e-2) goto reiterate;                          /* If result is not stable yet, increase number of iterations */
		if((phi_angle<0.000) || (phi_angle>360.000)) goto change_angle;  /* If result is unphysical, reduce the initial guess for the angle phi */
		}
/**/
/**/
	/* Calculate minimum and maximum angle for including left iron pond */
	phi_min_left    = impact4+180.000-phi_angle/2.000;
	phi_max_left    = impact4+180.000+phi_angle/2.000;
	status_phi_left = 0;
/**/
	/* Make sure angles are between 0 and 360 degrees */
	if((phi_min_left>=360.000) && (phi_max_left>=360.000))
		{
		phi_min_left    = phi_min_left-360.000;
		phi_max_left    = phi_max_left-360.000;
		status_phi_left = 0;
		}
	if((phi_max_left>=360.000) && (phi_min_left<360.000))
		{
		phi_max_left    = phi_max_left-360.000;
		status_phi_left = 1;
		}
/**/
	/* Add left iron pond */
	for(mm1=0;mm1<marknum;mm1++)
		{
		x=(double)(markx[mm1])/xsize;  /* non-dim. */
		y=(double)(marky[mm1])/ysize;  /* non-dim. */
		x1=x-xcoord_left;
		y1=(y-ycoord_left)*ysize/xsize;
/**/
		/* Upper boundary for the crater is the bottom of the ejecta zone */
		if((x1*x1+y1*y1)<=(corr*corr))
			{
			/* Lower boundary for the crater is the former surface */
			if((x1*x1+y1*y1)>=(0.0000*0.0000*ysize*ysize/xsize/xsize))
				{
				/* Calculate angles involved */
				ival=pow(x1*x1+y1*y1,0.500); if(!ival) ival=1.000;
				ival=x1/ival; ival=asin(ABSV(ival))/pival*180.000;
				if(x1>=0 && y1>0) ival=180.000-ival;
				if(x1<0 && y1>0) ival+=180.000;
				if(x1<0 && y1<=0) ival=360.000-ival;
				/**/
				/* Works out fine, when impact crater is not around 360 degrees */
				if((status_phi_left==0) && (ival>=ABSV(phi_min_left) && ival<=ABSV(phi_max_left)))
					{
					/* Markers in iron pond */
					markt[mm1]=(int)(iron_type[no1]);
					}
				/* Special case when iron pond is around 360 degrees */
				/* First fill the part 3xx to 360 degrees with markers */
				if((status_phi_left==1) && (ival>=ABSV(phi_min_left) && ival<=360.000))
					{
					/* Markers in iron pond */
					markt[mm1]=(int)(iron_type[no1]);
					}
				/* Now fill the part 0 to xx degrees with markers */
				if((status_phi_left==1) && (ival>=0.000 && ival<=ABSV(phi_max_left)))
					{
					/* Markers in iron pond */
					markt[mm1]=(int)(iron_type[no1]);
					}
				}
			}
		}
/**/
	/* Calculate minimum and maximum angle for including right iron pond */
	phi_min_right    = impact3+180.000-phi_angle/2.000;
	phi_max_right    = impact3+180.000+phi_angle/2.000;
	status_phi_right = 0;
/**/
	/* Make sure angles are between 0 and 360 degrees */
	if((phi_min_right>=360.000) && (phi_max_right>=360.000))
		{
		phi_min_right    = phi_min_right-360.000;
		phi_max_right    = phi_max_right-360.000;
		status_phi_right = 0;
		}
	if((phi_max_right>=360.000) && (phi_min_right<360.000))
		{
		phi_max_right    = phi_max_right-360.000;
		status_phi_right = 1;
		}
/**/
	/* Add right iron pond */
	for(mm1=0;mm1<marknum;mm1++)
		{
		x=(double)(markx[mm1])/xsize;  /* non-dim. */
		y=(double)(marky[mm1])/ysize;  /* non-dim. */
		x1=x-xcoord_right;
		y1=(y-ycoord_right)*ysize/xsize;
/**/
		/* Upper boundary for the crater is the bottom of the ejecta zone */
		if((x1*x1+y1*y1)<=(corr*corr))
			{
			/* Lower boundary for the crater is the former surface */
			if((x1*x1+y1*y1)>=(0.0000*0.0000*ysize*ysize/xsize/xsize))
				{
				/* Calculate angles involved */
				ival=pow(x1*x1+y1*y1,0.500); if(!ival) ival=1.000;
				ival=x1/ival; ival=asin(ABSV(ival))/pival*180.000;
				if(x1>=0 && y1>0) ival=180.000-ival;
				if(x1<0 && y1>0) ival+=180.000;
				if(x1<0 && y1<=0) ival=360.000-ival;
				/**/
				/* Works out fine, when impact crater is not around 360 degrees */
				if((status_phi_right==0) && (ival>=ABSV(phi_min_right) && ival<=ABSV(phi_max_right)))
					{
					/* Markers in iron pond */
					markt[mm1]=(int)(iron_type[no1]);
					}
				/* Special case when iron pond is around 360 degrees */
				/* First fill the part 3xx to 360 degrees with markers */
				if((status_phi_right==1) && (ival>=ABSV(phi_min_right) && ival<=360.000))
					{
					/* Markers in iron pond */
					markt[mm1]=(int)(iron_type[no1]);
					}
				/* Now fill the part 0 to xx degrees with markers */
				if((status_phi_right==1) && (ival>=0.000 && ival<=ABSV(phi_max_right)))
					{
					/* Markers in iron pond */
					markt[mm1]=(int)(iron_type[no1]);
					}
				}
			}
		}
/**/
	/* Temperature of LEFT iron pond (retains certain part of thermal history of impactor body) */
	for (m1=0;m1<xnumx;m1++)
	for (m2=0;m2<ynumy;m2++)
		{
		/* Node Num */
		m3=m1*ynumy+m2;
		/* Node X,Y General coordinates calc */
		x=gx[m1]/xsize;
		y=gy[m2]/ysize;
/**/
		x1=x-xcoord_left;
		y1=(y-ycoord_left)*ysize/xsize;
/**/
		x4=pow((x1*x1+y1*y1),0.500);
		if(x4<=corr)
			{
			/* Calculate and check distance from second center */
			y4=pow((x1*x1+y1*y1),0.500);
			if(y4>=0.000*ysize/xsize)
				{
				/* Relative coordinates calc */
				x4+=y4;if(!x4) x4=1.000;
				x=y4/x4;
				x4=ABSV(360.000)-ABSV(0.000); if(!x4) x4=1.000;
				y=(ival-ABSV(0.000))/x4;
                                /* First erase the temperature on that specific node */
                                tk[m3]=0.0000;
				/* Now put the temperature of the newly introduced LEFT iron pond on the specific node */
				tk[m3]=(memory_fe*t_impactor_core)+t_space;
				}
			}
		}
/**/
	/* Temperature of RIGHT iron pond (retains certain part of thermal history of impactor body) */
	for (m1=0;m1<xnumx;m1++)
	for (m2=0;m2<ynumy;m2++)
		{
		/* Node Num */
		m3=m1*ynumy+m2;
		/* Node X,Y General coordinates calc */
		x=gx[m1]/xsize;
		y=gy[m2]/ysize;
/**/
		x1=x-xcoord_right;
		y1=(y-ycoord_right)*ysize/xsize;
/**/
		x4=pow((x1*x1+y1*y1),0.500);
		if(x4<=corr)
			{
			/* Calculate and check distance from second center */
			y4=pow((x1*x1+y1*y1),0.500);
			if(y4>=0.000*ysize/xsize)
				{
				/* Relative coordinates calc */
				x4+=y4;if(!x4) x4=1.000;
				x=y4/x4;
				x4=ABSV(360.000)-ABSV(0.000); if(!x4) x4=1.000;
				y=(ival-ABSV(0.000))/x4;
                                /* First erase the temperature on that specific node */
                                tk[m3]=0.0000;
				/* Now put the temperature of the newly introduced RIGHT iron pond on the specific node */
				tk[m3]=(memory_fe*t_impactor_core)+t_space;
				}
			}
		}
/**/
/**/
	/* Calculate distance between center of planetary body and surface of magma pond (for left ejecta zone) */
	height_left  = pow(((xcoord_left-0.500)*(xcoord_left-0.500)+(ycoord_left-0.500)*(ycoord_left-0.500)),0.500)-ABSV(corr*cos((phi_angle/2.000)*(pival/180.000)));
/**/
	/* Distance pre-impact surface - planetary center */
	surface_left = pow(((xcoord_left-0.500)*(xcoord_left-0.500)+(ycoord_left-0.500)*(ycoord_left-0.500)),0.500);
/**/
	/* Half length of the surface of the iron layer */
	length_left      = ABSV(corr*sin(phi_angle/(2.000*180.000)*pival));
/**/
	/* Half angle of pond as seen from planetary center */
	pond_angle_left  = ABSV(atan(length_left/height_left)*180.000/pival);
/**/
	pond_min_left    = impact4-pond_angle_left;
	pond_max_left    = impact4+pond_angle_left;
	status_pond_left = 0;
/**/
	if(pond_min_left<=0.000)
		{
		pond_min_left = pond_min_left+360.000;
		status_pond_left = 1;
		}
/**/
	if(pond_max_left>=360.000)
		{
		pond_max_left = pond_max_left-360.000;
		status_pond_left = 1;
		}
/**/
	/* Remove the unnecessary part of the iron plus additional factor to make sure that all the iron in this zone is gone */
	for(mm1=0;mm1<marknum;mm1++)
		{
		x=(double)(markx[mm1])/xsize;  /* non-dim. */
		y=(double)(marky[mm1])/ysize;  /* non-dim. */
		x1=x-0.50000;
		y1=(y-0.50000)*ysize/xsize;
/**/
		/* Upper boundary for the erasure of unneeded material is above the protoplanet's surface */
		if((x1*x1+y1*y1)<=(1.50*surface_left*1.50*surface_left))
			{
			/* Lower boundary for the crater is the surface of the iron pond */
			if((x1*x1+y1*y1)>=(height_left*height_left*ysize*ysize/xsize/xsize))
				{
				/* Calculate angles involved */
				ival=pow(x1*x1+y1*y1,0.500); if(!ival) ival=1.000;
				ival=x1/ival; ival=asin(ABSV(ival))/pival*180.000;
				if(x1>=0 && y1>0) ival=180.000-ival;
				if(x1<0 && y1>0) ival+=180.000;
				if(x1<0 && y1<=0) ival=360.000-ival;
				/**/
				/* Works out fine, when iron pond is not around 360 degrees */
				if(status_pond_left == 0 && (ival>=ABSV(pond_min_left) && ival<=ABSV(pond_max_left)))
					{
					/* Unneeded markers in iron pond are transformed into sticky air */
					markt[mm1]=0;
					}
				/* Special case when iron pond is around 360 degrees */
				/* First fill the part 3xx to 360 degrees with markers */
				if(status_pond_left == 1 && (ival>=ABSV(pond_min_left) && ival<=360.000))
					{
					/* Unneeded markers in iron pond are transformed into sticky air */
					markt[mm1]=0;
					}
				/* Now fill the part 0 to xx degrees with markers */
				if(status_pond_left == 1 && (ival>=0.000 && ival<=ABSV(pond_max_left)))
					{
					/* Unneeded markers in iron pond are transformed into sticky air */
					markt[mm1]=0;
					}
				}
			}
		}
/**/
/* Calculate angles to put silicate ejecta layer on top of the iron pond */
	layer_min_left    = impact4-(crater_angle/2.000);
	layer_max_left    = impact4+(crater_angle/2.000);
	status_layer_left = 0;
/**/
	if(layer_min_left<=0.000)
		{
		layer_min_left = layer_min_left+360.000;
		status_layer_left = 1;
		}
/**/
	if(layer_max_left>=360.000)
		{
		layer_max_left = layer_max_left-360.000;
		status_layer_left = 1;
		}
/**/
	/* Add silicate ejecta layers only if iron was added to the impactor mantle */
	if(fe_fract>0.000)
		{
		/* Add silicate ejecta material on top of the LEFT iron pond */
		for(mm1=0;mm1<marknum;mm1++)
			{
			x=(double)(markx[mm1])/xsize;  /* non-dim. */
			y=(double)(marky[mm1])/ysize;  /* non-dim. */
			x1=x-0.50000;
			y1=(y-0.50000)*ysize/xsize;
/**/
			/* Upper boundary for the crater is the new protoplanet surface */
                	if((x1*x1+y1*y1)<=((height_left+d_layer)*(height_left+d_layer)))
				{
				/* Lower boundary for the crater is the surface of the iron pond */
				if((x1*x1+y1*y1)>=(height_left*height_left*ysize*ysize/xsize/xsize))
					{
					/* Calculate angles involved */
					ival=pow(x1*x1+y1*y1,0.500); if(!ival) ival=1.000;
					ival=x1/ival; ival=asin(ABSV(ival))/pival*180.000;
					if(x1>=0 && y1>0) ival=180.000-ival;
					if(x1<0 && y1>0) ival+=180.000;
					if(x1<0 && y1<=0) ival=360.000-ival;
					/**/
					/* Works out fine, when iron pond is not around 360 degrees */
					if(status_layer_left == 0 && (ival>=ABSV(layer_min_left) && ival<=ABSV(layer_max_left)))
						{
						/* Silicate ejecta material on top of iron pond */
						markt[mm1]=5;
						}
					/* Special case when iron pond is around 360 degrees */
					/* First fill the part 3xx to 360 degrees with markers */
					if(status_layer_left == 1 && (ival>=ABSV(layer_min_left) && ival<=360.000))
						{
						/* Silicate ejecta material on top of iron pond */
						markt[mm1]=5;
						}
					/* Now fill the part 0 to xx degrees with markers */
					if(status_layer_left == 1 && (ival>=0.000 && ival<=ABSV(layer_max_left)))
						{
						/* Silicate ejecta material on top of iron pond */
						markt[mm1]=5;
						}
					}
				}
			}
/**/
		/* Temperature of LEFT ejecta layer (retains certain part of thermal history of impactor body) */
		for (m1=0;m1<xnumx;m1++)
		for (m2=0;m2<ynumy;m2++)
			{
			/* Node Num */
			m3=m1*ynumy+m2;
			/* Node X,Y General coordinates calc */
			x=gx[m1]/xsize;
			y=gy[m2]/ysize;
/**/
			x1=x-0.50000;
			y1=(y-0.50000)*ysize/xsize;
/**/
			x4=pow((x1*x1+y1*y1),0.500);
			if(x4<=(height_left+d_layer))
				{
				/* Calculate and check distance from second center */
				y4=pow((x1*x1+y1*y1),0.500);
				if(y4>=height_left*ysize/xsize)
					{
					/* Calculate and check angle relatively first or second center */
					x3=x1;y3=y1; if(layer_max_left<0) {x3=x2; y3=y2;}
					ival=pow(x3*x3+y3*y3,0.500); if(!ival) ival=1.000;
					ival=x3/ival; ival=asin(ABSV(ival))/pival*180.000;
					if(x3>=0 && y3>0) ival=180.000-ival;
					if(x3<0 && y3>0) ival+=180.000;
					if(x3<0 && y3<=0) ival=360.000-ival;
                                	/**/
					/* Put the silicate ejecta layer over the LEFT iron pond */
					if(status_layer_left == 0 && (ival>=ABSV(layer_min_left) && ival<=ABSV(layer_max_left)))
						{
						/* Relative coordinates calc */
						x4+=y4;if(!x4) x4=1.000;
						x=y4/x4;
						x4=ABSV(layer_max_left)-ABSV(layer_min_left); if(!x4) x4=1.000;
						y=(ival-ABSV(layer_max_left))/x4;
                                        	/* First erase the temperature on that specific node */
                                        	tk[m3]=0.0000;
						/* Now put the temperature of the newly introduced LEFT silicate ejecta layer on the specific node */
						tk[m3]=(memory_si*t_blanket)+t_space;
						}
					/* Special case when ejecta layer is around 360 degrees */
					/* First change the temperature in section 3xx to 360 degrees */
					if(status_layer_left == 1 && (ival>=ABSV(layer_min_left) && ival<=360.000))
						{
						/* Relative coordinates calc */
						x4+=y4;if(!x4) x4=1.000;
						x=y4/x4;
						x4=ABSV(layer_max_left)-ABSV(layer_min_left); if(!x4) x4=1.000;
						y=(ival-ABSV(layer_max_left))/x4;
                                        	/* First erase the temperature on that specific node */
                                        	tk[m3]=0.0000;
						/* Now put the temperature of the newly introduced LEFT silicate ejecta layer on the specific node */
						tk[m3]=(memory_si*t_blanket)+t_space;
						}
					/* Now change temperature in section 0 to xx degrees */
					if(status_layer_left == 1 && (ival>=0.000 && ival<=ABSV(layer_max_left)))
						{
 						/* Relative coordinates calc */
						x4+=y4;if(!x4) x4=1.000;
						x=y4/x4;
						x4=ABSV(layer_max_left)-ABSV(layer_min_left); if(!x4) x4=1.000;
						y=(ival-ABSV(layer_max_left))/x4;
                                        	/* First erase the temperature on that specific node */
                                        	tk[m3]=0.0000;
						/* Now put the temperature of the newly introduced LEFT silicate ejecta layer on the specific node */
						tk[m3]=(memory_si*t_blanket)+t_space;
						}
					}
				}	
			}
		/* End of (fe_fract>0) if loop */
		}
/**/
/**/
	/* Calculate distance between center of planetary body and surface of magma pond (for right ejecta zone) */
	height_right  = pow(((xcoord_right-0.500)*(xcoord_right-0.500)+(ycoord_right-0.500)*(ycoord_right-0.500)),0.500)-ABSV(corr*cos((phi_angle/2.000)*(pival/180.000)));
/**/
	/* Distance pre-impact surface - planetary center */
	surface_right = pow(((xcoord_right-0.500)*(xcoord_right-0.500)+(ycoord_right-0.500)*(ycoord_right-0.500)),0.500);
/**/
	/* Half length of the surface of the iron layer */
	length_right      = ABSV(corr*sin(phi_angle/(2.000*180.000)*pival));
/**/
	/* Half angle of pond as seen from planetary center */
	pond_angle_right  = ABSV(atan(length_right/height_right)*180.000/pival);
/**/
	pond_min_right    = impact3-pond_angle_right;
	pond_max_right    = impact3+pond_angle_right;
	status_pond_right = 0;
/**/
	if(pond_min_right<=0.000)
		{
		pond_min_right = pond_min_right+360.000;
		status_pond_right = 1;
		}
/**/
	if(pond_max_right>=360.000)
		{
		pond_max_right = pond_max_right-360.000;
		status_pond_right = 1;
		}
/**/
	/* Remove the unnecessary part of the iron plus additional factor to make sure that all the iron in this zone is gone */
	for(mm1=0;mm1<marknum;mm1++)
		{
		x=(double)(markx[mm1])/xsize;  /* non-dim. */
		y=(double)(marky[mm1])/ysize;  /* non-dim. */
		x1=x-0.50000;
		y1=(y-0.50000)*ysize/xsize;
/**/
		/* Upper boundary for the erasure of unneeded material is above the protoplanet's surface */
		if((x1*x1+y1*y1)<=(1.50*surface_right*1.50*surface_right))
			{
			/* Lower boundary for the crater is the bottom of the impact crater */
			if((x1*x1+y1*y1)>=(height_right*height_right*ysize*ysize/xsize/xsize))
				{
				/* Calculate angles involved */
				ival=pow(x1*x1+y1*y1,0.500); if(!ival) ival=1.000;
				ival=x1/ival; ival=asin(ABSV(ival))/pival*180.000;
				if(x1>=0 && y1>0) ival=180.000-ival;
				if(x1<0 && y1>0) ival+=180.000;
				if(x1<0 && y1<=0) ival=360.000-ival;
				/**/
				/* Works out fine, when iron pond is not around 360 degrees */
				if(status_pond_right == 0 && (ival>=ABSV(pond_min_right) && ival<=ABSV(pond_max_right)))
					{
					/* Unneeded markers in iron pond are transformed into sticky air */
					markt[mm1]=0;
					}
				/* Special case when iron pond is around 360 degrees */
				/* First fill the part 3xx to 360 degrees with markers */
				if(status_pond_right == 1 && (ival>=ABSV(pond_min_right) && ival<=360.000))
					{
					/* Unneeded markers in iron pond are transformed into sticky air */
					markt[mm1]=0;
					}
				/* Now fill the part 0 to xx degrees with markers */
				if(status_pond_right == 1 && (ival>=0.000 && ival<=ABSV(pond_max_right)))
					{
					/* Unneeded markers in iron pond are transformed into sticky air */
					markt[mm1]=0;
					}
				}
			}
		}
/**/
/**/
/* Calculate angles to put silicate ejecta layer on top of the iron pond */
	layer_min_right    = impact3-(crater_angle/2.000);
	layer_max_right    = impact3+(crater_angle/2.000);
	status_layer_right = 0;
/**/
	if(layer_min_right<=0.000)
		{
		layer_min_right = layer_min_right+360.000;
		status_layer_right = 1;
		}
/**/
	if(layer_max_right>=360.000)
		{
		layer_max_right = layer_max_right-360.000;
		status_layer_right = 1;
		}
/**/
	/* Only add the ejecta layer when iron was contributed to the impactor mantle */
	if(fe_fract>0.000)
		{
		/* Add silicate ejecta material on top of the RIGHT iron pond */
		for(mm1=0;mm1<marknum;mm1++)
			{
			x=(double)(markx[mm1])/xsize;  /* non-dim. */
			y=(double)(marky[mm1])/ysize;  /* non-dim. */
			x1=x-0.50000;
			y1=(y-0.50000)*ysize/xsize;
/**/
			/* Upper boundary for the crater is the new protoplanet surface */
			if((x1*x1+y1*y1)<=((height_right+d_layer)*(height_right+d_layer)))
				{
				/* Lower boundary for the crater is the surface of the iron pond */
				if((x1*x1+y1*y1)>=(height_right*height_right*ysize*ysize/xsize/xsize))
					{
					/* Calculate angles involved */
					ival=pow(x1*x1+y1*y1,0.500); if(!ival) ival=1.000;
					ival=x1/ival; ival=asin(ABSV(ival))/pival*180.000;
					if(x1>=0 && y1>0) ival=180.000-ival;
					if(x1<0 && y1>0) ival+=180.000;
					if(x1<0 && y1<=0) ival=360.000-ival;
					/**/
					/* Works out fine, when iron pond is not around 360 degrees */
					if(status_layer_right == 0 && (ival>=ABSV(layer_min_right) && ival<=ABSV(layer_max_right)))
						{
						/* Silicate ejecta material on top of iron pond */
						markt[mm1]=5;
						}
					/* Special case when iron pond is around 360 degrees */
					/* First fill the part 3xx to 360 degrees with markers */
					if(status_layer_right == 1 && (ival>=ABSV(layer_min_right) && ival<=360.000))
						{
						/* Silicate ejecta material on top of iron pond */
						markt[mm1]=5;
						}
					/* Now fill the part 0 to xx degrees with markers */
					if(status_layer_right == 1 && (ival>=0.000 && ival<=ABSV(layer_max_right)))
						{
						/* Silicate ejecta material on top of iron pond */
						markt[mm1]=5;
						}
					}
				}
			}
/**/
		/* Temperature of RIGHT ejecta layer (retains certain part of thermal history of impactor body) */
		for (m1=0;m1<xnumx;m1++)
		for (m2=0;m2<ynumy;m2++)
			{
			/* Node Num */
			m3=m1*ynumy+m2;
			/* Node X,Y General coordinates calc */
			x=gx[m1]/xsize;
			y=gy[m2]/ysize;
/**/
			x1=x-0.50000;
			y1=(y-0.50000)*ysize/xsize;
/**/
			x4=pow((x1*x1+y1*y1),0.500);
			if(x4<=(height_right+d_layer))
				{
				/* Calculate and check distance from second center */
				y4=pow((x1*x1+y1*y1),0.500);
				if(y4>=height_right*ysize/xsize)
					{
					/* Calculate and check angle relatively first or second center */
					x3=x1;y3=y1; if(layer_max_right<0) {x3=x2; y3=y2;}
					ival=pow(x3*x3+y3*y3,0.500); if(!ival) ival=1.000;
					ival=x3/ival; ival=asin(ABSV(ival))/pival*180.000;
					if(x3>=0 && y3>0) ival=180.000-ival;
					if(x3<0 && y3>0) ival+=180.000;
					if(x3<0 && y3<=0) ival=360.000-ival;
                                	/**/
					/* Put the silicate ejecta layer over the RIGHT iron pond */
					if(status_layer_right == 0 && (ival>=ABSV(layer_min_right) && ival<=ABSV(layer_max_right)))
						{
						/* Relative coordinates calc */
						x4+=y4;if(!x4) x4=1.000;
						x=y4/x4;
						x4=ABSV(layer_max_right)-ABSV(layer_min_right); if(!x4) x4=1.000;
						y=(ival-ABSV(layer_max_right))/x4;
                                        	/* First erase the temperature on that specific node */
                                        	tk[m3]=0.0000;
						/* Now put the temperature of the newly introduced RIGHT silicate ejecta layer on the specific node */
						tk[m3]=(memory_si*t_blanket)+t_space;
						}
					/* Special case when ejecta layer is around 360 degrees */
					/* First change the temperature in section 3xx to 360 degrees */
					if(status_layer_right == 1 && (ival>=ABSV(layer_min_right) && ival<=360.000))
						{
						/* Relative coordinates calc */
						x4+=y4;if(!x4) x4=1.000;
						x=y4/x4;
						x4=ABSV(layer_max_right)-ABSV(layer_min_right); if(!x4) x4=1.000;
						y=(ival-ABSV(layer_max_right))/x4;
                                        	/* First erase the temperature on that specific node */
                                        	tk[m3]=0.0000;
						/* Now put the temperature of the newly introduced RIGHT silicate ejecta layer on the specific node */
						tk[m3]=(memory_si*t_blanket)+t_space;
						}
					/* Now change temperature in section 0 to xx degrees */
					if(status_layer_right == 1 && (ival>=0.000 && ival<=ABSV(layer_max_right)))
						{
 						/* Relative coordinates calc */
						x4+=y4;if(!x4) x4=1.000;
						x=y4/x4;
						x4=ABSV(layer_max_right)-ABSV(layer_min_right); if(!x4) x4=1.000;
						y=(ival-ABSV(layer_max_right))/x4;
                                        	/* First erase the temperature on that specific node */
                                        	tk[m3]=0.0000;
						/* Now put the temperature of the newly introduced RIGHT silicate ejecta layer on the specific node */
						tk[m3]=(memory_si*t_blanket)+t_space;
						}
					}
				}
			}
		/* End of (fe_fract>0) loop */
		}
/**/
/*===========================================================================================================================*/
/**/
	/* Material in impacted region is weakened for a certain time */
	for(mm1=0;mm1<marknum;mm1++)
		{
		x=(double)(markx[mm1])/xsize;  /* non-dim. */
		y=(double)(marky[mm1])/ysize;  /* non-dim. */
		x1=x-xcoordsurf;
		y1=(y-ycoordsurf)*ysize/xsize;
/**/
		/* Weakened region encompasses both ejecta zones and extends 1.5 times the crater radius (assuming pond mode) */
		if((x1*x1+y1*y1)<=((3.000*corr1)*(3.000*corr1)))
			{
			/* Lower boundary for the crater is the bottom of the impact crater */
			if((x1*x1+y1*y1)>=(0.0000*0.0000*ysize*ysize/xsize/xsize))
				{
				/* Calculate angles involved */
				ival=pow(x1*x1+y1*y1,0.500); if(!ival) ival=1.000;
				ival=x1/ival; ival=asin(ABSV(ival))/pival*180.000;
				if(x1>=0 && y1>0) ival=180.000-ival;
				if(x1<0 && y1>0) ival+=180.000;
				if(x1<0 && y1<=0) ival=360.000-ival;
				/**/
				/* Only silicates are weakened */
        			if(markt[mm1]==5 || markt[mm1]==6)
					{
					/* Note when the material was weakened to make it again stiff after relax_time passes */
					markhi[mm1]=impact_time[no1];  /* weakening time in [Ma] */
					}
				}
			}
		}
/**/
/*===========================================================================================================================*/
/**/
        /* Calculate new maximum thickness of iron pond [m] to be used in merger calculations in core.c and move positions of old values */
        d_iron_old = d_iron_ol;
        d_iron_ol  = d_iron;
/**/
	/* Give d_iron only a new positive value when new iron was added to the impactor mantle */	 
	d_iron     = 0.000;
/**/
	if(fe_fract>0.000)
		{
		d_iron = (height_left-pow(((xcoord_left-0.5000)*(xcoord_left-0.5000)+(ycoord_left-0.5000)*(ycoord_left-0.5000)),0.5000)+corr)*xsize;
		}
/**/
/*===========================================================================================================================*/
/**/
	/* Make sure code knows timing of next impact event! */
	/* Set impact counter to the next value */
        no1=no1+1;
/**/
/**/
	/* Recalculate temperature for markers */
	for (mm1=0;mm1<marknum;mm1++)
		/* Check markers out of grid */
		if ((double)(markx[mm1])>0 && (double)(marky[mm1])>0 && (double)(markx[mm1])<xsize && (double)(marky[mm1])<ysize)
		{
		allintert((double)(markx[mm1]),(double)(marky[mm1]));
		markk[mm1]=(float)(eps[2]);
		}
/**/
	/* Recalculate densities and viscosities */
	ronurecalc();
/**/
/**/
	/* Now close the big accretion loop */
	}
/**/
/**/
/* Hydration/dehydration routine */
/**/
/* Define counters and current hydrated fraction */
count     = 0;
count1    = 0;
hydr_frac = 10.000;  /* In case of malfunction or deactivation easy to spot */
/**/
/* Make sure that hydration/dehydration is only active when no olivine grain growth is computed at the same time and no initial porosity is considered */
if(growth_model==0 && por_init<=0.001)
        {
        /* Check for each timestep whether silicate markers would present rock+ice, hydrated silicates or dehydrated silicates */
        for(mm1=0;mm1<marknum;mm1++)
                {
                /* If temperature is sufficienty high silicate material will hydrate */
                if(markt[mm1]==5 && markk[mm1]>273.15)
                        {
                        markt[mm1]=6;
                        }
                /* If temperature is sufficiently high silicate material will dehydrate */
                if(markt[mm1]==6 && markk[mm1]>273.15+950.00)
                        {
                        markt[mm1]=5;
                        }
                /**/
                /* Count how many silicate markers are currently hydrous */
                if(markt[mm1]==6)
                        {
                        count+=1;
                        }
                /**/
                /* Count how many markers are currently chondritic or silicatic */
                if(markt[mm1]==5 || markt[mm1]==6 || markt[mm1]==25 || markt[mm1]==26)
                        {
                        count1+=1;
                        }
                }
        /**/
        /* Compute current fraction of chondritic markers */
        if(count1>0)
                {
                hydr_frac = (double)(count)/(double)(count1);
                }
        }
/**/
return 0;
/**/
/* End subroutine impact */
}
