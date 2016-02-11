/* Treat iron core evolution and planetary dynamo during planetary accretion */
/* by Greg (last updated: 08/05/2012) */
/**/
int core()
{
	/* Define parameters */
        int merger;
	long int ic,ic_max,m1,m2,m3,mm1,ml[1000],mll[1000],no_core,no_diapir,no_iron,stack_count;
	double alpha_core,alpha_mantle,b_dip,b_rms,core_angle,cP_core,cP_mantle,crit_dens1,crit_dens2,crit_dens3,diff_fe,distance[1000],ek,eta_fe;
	double g_core,G,ival,iron_layer,k_core,k_mantle[1000],l_nor,l_nor1,M_iron;
	double magn_diff,magn_perm,max_radius_core[1000],max_radius_diapir,min_radius_diapir,max_radius_planet,mean_radius_core,min_radius_mantle_local[1000];
	double pi,ppu,prandtl,pr_center,prm,q_a,q_adv,q_c,q_cond,q_adv_stack,q_cond_stack,ra_q,rho_core,rho_mantle,stack,step,t_above[1000],t_ad_core,t_ad_diapir;
        double t_core,t_cmb[1000],t_cmb_stack,t_diapir,t_mean_core,t_mean_diapir,t_melting,test_angle,test_center,test_center2,test_distance,uncert;
	double variation,v_marker,v_rad,vx_marker,vy_marker,vx_rad,vy_rad,x_dist,x_nor,x_nor1[1000],y_dist,y_nor,y_nor1[1000],x,x1,y,y1,year_to_sec,zeta;
	/**/
	/* Set constants and define parameters */
	pi          = 3.141592654;              /* Define Pi */
	G           = 6.67384e-11;              /* Gravitational constant [m^3/(kg*s^2)] */
	variation   = 0.500e0;                  /* Allowed deviation between core_angle and test_angle [degrees] */
        crit_dens1  = 6000.000;                 /* Critical density, grids with higher density are most likely in iron core of target body [kg/m^3] */
        crit_dens2  = 4500.000;                 /* Critical density, grids with higher density are most likely iron enriched, not pure silicates [kg/m^3] */
        crit_dens3  = 1500.000;                 /* Critical density, grids with lower density are most likely sticky air material [kg/m^3] */
	year_to_sec = 365.250*24.000*3600.000;  /* Conversion factor year to seconds */
	/**/
        /* Distance (uncert*mean_radius_core) is used to check whether merger happened and to modify marker type when merger happens */
	/* First collisions may move planetary body in numerical box, thus the criterion is less strict in the beginning */
        if(no1<=2)
		{
        	uncert = 1.150;
		}
        else if(no1>2)
		{
		uncert = 1.100;
		}
	/**/
	/* Obtain various physical parameters */
	alpha_core = markbb[8];    /* Thermal expansivity of iron material [1/K] */
	cP_mantle  = markcp[5];    /* Heat capacity of solid/molten silicate material [J/(K*kg)] */  
	cP_core    = markcp[8];    /* Heat capacity of iron material [J/(K*kg)] */
        k_core     = 46.000;       /* Thermal conductivity of iron material [W/(m*K)] (C.A. Jones, Treatise on Geophysics, Vol. 8, 131-185, 2007 => T.M. Rogers et al., Physical Review E 67, 026315, 2003) */
	/**/
	step       = 1.000e0;      /* Search sector width [degrees] */
	/**/
	/* Predefine heat flux parameters [W/m^2] */
	q_mantle   = 0.00000;
	q_core     = 0.00000;
	/**/
	/* Assume as standard that the planetary dynamo is inactive */
	dynamo     = 0;
	/**/
	/* Additional iron core properties taken from (C.A. Jones, Treatise on Geophysics, Vol. 8, 131-185, 2007) */
	/* Magnetic diffusivity of molten iron [m^2/s] */
	magn_diff  = 2.000;
	/**/
 	/* Magnetic permeability of molten iron [H/m = N/A^2] (core is hot, above the Curie point => the magnetic permeability is essentially that of free space) */
	magn_perm  = 4.000*pi*1.000e-7;
	/**/
	/* Diffusivity in molten iron [m^2/s] */
	diff_fe    = 5.000e-6;
	/**/
	/* Viscosity of molten iron [Pa s] (D.C. Rubie et al., Treatise on Geophysics Vol. 9, 51-90, 2007) */
	eta_fe     = 0.010;
	/**/
	/* Assume as standard a magnetic dipole moment of zero [A*m^2] */
	dipl_mom   = 0.00000;
	/**/
	/* Assume as standard a magnetic Reynolds number of zero [non-dim.] */
	reynolds   = 0.00000;
	/**/
	/* Assume as standard a local Rossby number of zero [non-dim.] */
	rossby     = 0.00000;
	/**/
/*==============================================================================================================================================================================================*/
	/**/
	/* Check geodynamo activity only for time>0.0000 */
	if(timesum>0.0000)
	{
	/**/
	/* Calculate MEAN density of the target body iron core [kg/m^3] (considering compressional effects) */
	stack       = 0.0000;
	stack_count = 0;
	/**/
	for(m1=0;m1<xnumx;m1++)
	for(m2=0;m2<ynumy;m2++)
		{
		/* Node number */
		m3=m1*ynumy+m2;
		/**/
		if(ro[m3]>crit_dens1)
			{
			stack+=ro[m3];
			stack_count=stack_count+1;
			}
		}
	rho_core = stack/(double)(stack_count);
	/**/
	/* Calculate MEAN density of the target body silicate mantle [kg/m^3] (considering compressional effects) */
	stack       = 0.0000;
	stack_count = 0;
	/**/
	for(m1=0;m1<xnumx;m1++)
	for(m2=0;m2<ynumy;m2++)
		{
		/* Node number */
		m3=m1*ynumy+m2;
		/**/
		if(ro[m3]<crit_dens2 && ro[m3]>crit_dens3)
			{
			stack+=ro[m3];
			stack_count=stack_count+1;
			}
		}
	rho_mantle = stack/(double)(stack_count);
	/**/
	/**/
	/* I.) Find current mean radius of the iron core [non-dim.] */
        /**/
	/* Set counter and other parameters initially to zero */
	ic       = 0;          /* Counter */
	g_core   = 0.000000;   /* Used to measure gravity acceleration of the upper bound of preexistent iron core [m/s^2] */
        stack    = 0.000000;   /* Used to stack maximum core radii from the different sectors [non-dim.] */
	/**/
	/* Make sure not to count the zero position twice */
	for(core_angle=0.000000;core_angle<=(360.000-step);core_angle=core_angle+step)
		{
		/* Reset counter */
		ic = ic+1;
		/**/
                /* Set maximum iron core radius for each sector initially to zero */
	        max_radius_core[ic] = 0.000000;
                /**/
                /* Reset distance test value for each sector */
                test_center = 0.000000;
		/**/
		/* Go through all markers */
		for(mm1=0;mm1<marknum;mm1++)
			{
                        /* Reset test_angle always to zero */
                        test_angle  = 0.000000;
                        /**/
			/* Check only preexistent iron core markers, all other materials are excluded */
			/* Estimate radius of preexistent iron core [non-dim.] */
        		if(markt[mm1]==10)
				{
				/* Angle 0 <= alpha < 90 degrees */
                       		if((core_angle>=0.000 && core_angle<90.000) && (markx[mm1]/xsize>=0.500) && (marky[mm1]/ysize<0.500))
					{
                               		test_angle = atan((markx[mm1]/xsize-0.500)/(0.500-marky[mm1]/ysize))*180.000/pi;
                               		if(ABSV(test_angle-core_angle)<=variation)
						{
						test_center=pow(((markx[mm1]/xsize-0.500)*(markx[mm1]/xsize-0.500)+(marky[mm1]/ysize-0.500)*(marky[mm1]/ysize-0.500)),0.500);
						/**/
                       				max_radius_core[ic]=MAXV(test_center,max_radius_core[ic]);
						}
					}
				/* Angle 90 <= alpha < 180 degrees */
                       		else if((core_angle>=90.000 && core_angle<180.000) && (markx[mm1]/xsize>0.500) && (marky[mm1]/ysize>=0.500))
					{
                               		test_angle = 90.000+atan((marky[mm1]/ysize-0.500)/(markx[mm1]/xsize-0.500))*180.000/pi;
                               		if(ABSV(test_angle-core_angle)<=variation)
						{
						test_center=pow(((markx[mm1]/xsize-0.500)*(markx[mm1]/xsize-0.500)+(marky[mm1]/ysize-0.500)*(marky[mm1]/ysize-0.500)),0.500);
						/**/
                       				max_radius_core[ic]=MAXV(test_center,max_radius_core[ic]);
						}
					}
				/* Angle 180 <= alpha < 270 degrees */
                       		else if((core_angle>=180.000 && core_angle<270.000) && (markx[mm1]/xsize<=0.500) && (marky[mm1]/ysize>0.500))
					{
                               		test_angle = 180.000+atan((0.500-markx[mm1]/xsize)/(marky[mm1]/ysize-0.500))*180.000/pi;
                               		if(ABSV(test_angle-core_angle)<=variation)
						{
						test_center=pow(((markx[mm1]/xsize-0.500)*(markx[mm1]/xsize-0.500)+(marky[mm1]/ysize-0.500)*(marky[mm1]/ysize-0.500)),0.500);
						/**/
                       				max_radius_core[ic]=MAXV(test_center,max_radius_core[ic]);
						}
					}
				/* Angle 270 <= alpha < 360 degrees */
                        	else if((core_angle>=270.000 && core_angle<360.000) && (markx[mm1]/xsize<0.500) && (marky[mm1]/ysize<=0.500))
					{
                                	test_angle = 270.000+atan((0.500-marky[mm1]/ysize)/(0.500-markx[mm1]/xsize))*180.000/pi;
                                	if(ABSV(test_angle-core_angle)<=variation)
						{
						test_center=pow(((markx[mm1]/xsize-0.500)*(markx[mm1]/xsize-0.500)+(marky[mm1]/ysize-0.500)*(marky[mm1]/ysize-0.500)),0.500);
						/**/
                       				max_radius_core[ic]=MAXV(test_center,max_radius_core[ic]);
						}
					}
				}
                        }
		/**/
		/* Maximum number of steps */
		ic_max = ic;
                /**/
                /* Stack maximum core radii from the different sectors [non-dim.] */
                stack+=max_radius_core[ic];
		}
        /**/
	/* Compute current mean iron core radius [non-dim.] */
        mean_radius_core = stack/(double)(ic_max);
        /**/
	/* Compute approximate gravity acceleration at CMB [m/s^2] */
	g_core = (4.000/3.000)*pi*rho_core*G*mean_radius_core*xsize;
        /**/
        /* Compute the initial mass of the preexistent iron core before the first impact happens [kg] */
        if(no1==1)
		{
		M_core = (4.000/3.000)*pi*rho_core*pow((mean_radius_core*xsize),3.000);
		}
	/**/
	/* II.) Find minimum and maximum iron diapir radius [non-dim.] */
        /**/
	/* Set counter and other parameters initially to zero */
	max_radius_diapir = 0.000000;
	min_radius_diapir = 1.000000;      /* Initially set equal to non-dim. box size */
	no_iron           = 0;             /* Counter for possible iron markers in the mantle */
	/**/
	/* Make sure not to count the zero position twice */
	for(core_angle=0.000000;core_angle<=(360.000-step);core_angle=core_angle+step)
		{
		/**/
                /* Reset distance test value for each sector */
                test_center  = 0.000000;
                test_center2 = 1.000000;
                /**/
		/* Go through all markers */
		for(mm1=0;mm1<marknum;mm1++)
			{
			/**/
                        /* Reset test_angle always to zero */
                        test_angle  = 0.000000;
			/**/
			/* Now calculate the distance of possible iron diapirs from the preexistent iron core [non-dim.] */
        		if(markt[mm1]==(int)(iron_type[no_merger+1]))
				{
				/* Reset counter */
				no_iron+=1;
				/**/
				/* Angle 0 <= alpha < 90 degrees */
                       		if((core_angle>=0.000 && core_angle<90.000) && (markx[mm1]/xsize>=0.500) && (marky[mm1]/ysize<0.500))
					{
                               		test_angle = atan((markx[mm1]/xsize-0.500)/(0.500-marky[mm1]/ysize))*180.000/pi;
                               		if(ABSV(test_angle-core_angle)<=variation)
						{
						test_center=pow(((markx[mm1]/xsize-0.500)*(markx[mm1]/xsize-0.500)+(marky[mm1]/ysize-0.500)*(marky[mm1]/ysize-0.500)),0.500);
                                                test_center2=pow(((markx[mm1]/xsize-0.500)*(markx[mm1]/xsize-0.500)+(marky[mm1]/ysize-0.500)*(marky[mm1]/ysize-0.500)),0.500);
						/**/
                       				max_radius_diapir=MAXV(test_center,max_radius_diapir);
						min_radius_diapir=MINV(test_center2,min_radius_diapir);
						}
					}
				/* Angle 90 <= alpha < 180 degrees */
                       		else if((core_angle>=90.000 && core_angle<180.000) && (markx[mm1]/xsize>0.500) && (marky[mm1]/ysize>=0.500))
					{
                               		test_angle = 90.000+atan((marky[mm1]/ysize-0.500)/(markx[mm1]/xsize-0.500))*180.000/pi;
                               		if(ABSV(test_angle-core_angle)<=variation)
						{
						test_center=pow(((markx[mm1]/xsize-0.500)*(markx[mm1]/xsize-0.500)+(marky[mm1]/ysize-0.500)*(marky[mm1]/ysize-0.500)),0.500);
                                                test_center2=pow(((markx[mm1]/xsize-0.500)*(markx[mm1]/xsize-0.500)+(marky[mm1]/ysize-0.500)*(marky[mm1]/ysize-0.500)),0.500);
						/**/
                       				max_radius_diapir=MAXV(test_center,max_radius_diapir);
						min_radius_diapir=MINV(test_center2,min_radius_diapir);
						}
					}
				/* Angle 180 <= alpha < 270 degrees */
                       		else if((core_angle>=180.000 && core_angle<270.000) && (markx[mm1]/xsize<=0.500) && (marky[mm1]/ysize>0.500))
					{
                               		test_angle = 180.000+atan((0.500-markx[mm1]/xsize)/(marky[mm1]/ysize-0.500))*180.000/pi;
                               		if(ABSV(test_angle-core_angle)<=variation)
						{
						test_center=pow(((markx[mm1]/xsize-0.500)*(markx[mm1]/xsize-0.500)+(marky[mm1]/ysize-0.500)*(marky[mm1]/ysize-0.500)),0.500);
                                                test_center2=pow(((markx[mm1]/xsize-0.500)*(markx[mm1]/xsize-0.500)+(marky[mm1]/ysize-0.500)*(marky[mm1]/ysize-0.500)),0.500);
						/**/
                       				max_radius_diapir=MAXV(test_center,max_radius_diapir);
						min_radius_diapir=MINV(test_center2,min_radius_diapir);
						}
					}
				/* Angle 270 <= alpha < 360 degrees */
                        	else if((core_angle>=270.000 && core_angle<360.000) && (markx[mm1]/xsize<0.500) && (marky[mm1]/ysize<=0.500))
					{
                                	test_angle = 270.000+atan((0.500-marky[mm1]/ysize)/(0.500-markx[mm1]/xsize))*180.000/pi;
                                	if(ABSV(test_angle-core_angle)<=variation)
						{
						test_center=pow(((markx[mm1]/xsize-0.500)*(markx[mm1]/xsize-0.500)+(marky[mm1]/ysize-0.500)*(marky[mm1]/ysize-0.500)),0.500);
                                                test_center2=pow(((markx[mm1]/xsize-0.500)*(markx[mm1]/xsize-0.500)+(marky[mm1]/ysize-0.500)*(marky[mm1]/ysize-0.500)),0.500);
						/**/
                       				max_radius_diapir=MAXV(test_center,max_radius_diapir);
						min_radius_diapir=MINV(test_center2,min_radius_diapir);
						}
					}
				}
                        }
                }
	/**/
        /* III.) Find maximum planetary radius [non-dim.] */
        /**/
	/* Set counter and other parameters initially to zero */
	max_radius_planet = 0.000000;      /* Used to compute maximum radius of planetesimal [non-dim.] */
        /**/
	/* Make sure not to count the zero position twice */
	for(core_angle=0.000000;core_angle<=(360.000-step);core_angle=core_angle+step)
		{
		/**/
                /* Reset distance test value for each sector */
                test_center = 0.000000;
                /**/
                /* Go through all markers */
		for(mm1=0;mm1<marknum;mm1++)
			{
			/**/
                        /* Reset test_angle always to zero */
                        test_angle  = 0.000000;
			/**/
			/* Calculate TOTAL radius of planetesimal [non-dim.] */
			/* Check all marker types, except sticky air */
        		if(markt[mm1]==5 || markt[mm1]==6 || markt[mm1]==7 || markt[mm1]==8 || markt[mm1]==9 || markt[mm1]==10 || markt[mm1]==17 || markt[mm1]==18 || markt[mm1]==19 || markt[mm1]==25 || markt[mm1]==26)
				{
				/* Angle 0 <= alpha < 90 degrees */
                       		if((core_angle>=0.000 && core_angle<90.000) && (markx[mm1]/xsize>=0.500) && (marky[mm1]/ysize<0.500))
					{
                               		test_angle = atan((markx[mm1]/xsize-0.500)/(0.500-marky[mm1]/ysize))*180.000/pi;
                               		if(ABSV(test_angle-core_angle)<=variation)
						{
						test_center=pow(((markx[mm1]/xsize-0.500)*(markx[mm1]/xsize-0.500)+(marky[mm1]/ysize-0.500)*(marky[mm1]/ysize-0.500)),0.500);
						/**/
                       				max_radius_planet=MAXV(test_center,max_radius_planet);
						}
					}
				/* Angle 90 <= alpha < 180 degrees */
                       		else if((core_angle>=90.000 && core_angle<180.000) && (markx[mm1]/xsize>0.500) && (marky[mm1]/ysize>=0.500))
					{
                               		test_angle = 90.000+atan((marky[mm1]/ysize-0.500)/(markx[mm1]/xsize-0.500))*180.000/pi;
                               		if(ABSV(test_angle-core_angle)<=variation)
						{
						test_center=pow(((markx[mm1]/xsize-0.500)*(markx[mm1]/xsize-0.500)+(marky[mm1]/ysize-0.500)*(marky[mm1]/ysize-0.500)),0.500);
						/**/
                       				max_radius_planet=MAXV(test_center,max_radius_planet);
						}
					}
				/* Angle 180 <= alpha < 270 degrees */
                       		else if((core_angle>=180.000 && core_angle<270.000) && (markx[mm1]/xsize<=0.500) && (marky[mm1]/ysize>0.500))
					{
                               		test_angle = 180.000+atan((0.500-markx[mm1]/xsize)/(marky[mm1]/ysize-0.500))*180.000/pi;
                               		if(ABSV(test_angle-core_angle)<=variation)
						{
						test_center=pow(((markx[mm1]/xsize-0.500)*(markx[mm1]/xsize-0.500)+(marky[mm1]/ysize-0.500)*(marky[mm1]/ysize-0.500)),0.500);
						/**/
                       				max_radius_planet=MAXV(test_center,max_radius_planet);
						}
					}
				/* Angle 270 <= alpha < 360 degrees */
                        	else if((core_angle>=270.000 && core_angle<360.000) && (markx[mm1]/xsize<0.500) && (marky[mm1]/ysize<=0.500))
					{
                                	test_angle = 270.000+atan((0.500-marky[mm1]/ysize)/(0.500-markx[mm1]/xsize))*180.000/pi;
                                	if(ABSV(test_angle-core_angle)<=variation)
						{
						test_center=pow(((markx[mm1]/xsize-0.500)*(markx[mm1]/xsize-0.500)+(marky[mm1]/ysize-0.500)*(marky[mm1]/ysize-0.500)),0.500);
						/**/
                       				max_radius_planet=MAXV(test_center,max_radius_planet);
						}
					}
				}
			}
		}
	/**/
/*==============================================================================================================================================================================================*/
	/**/
	/* Merger status */
	/**/
        /** No iron diapir markers in the mantle or no impact happened yet */
        if((no_iron==0) || (no1==1))
		{
		merger = 0;
                iron_layer = 0.000;
		}
        /**/
	/* Iron diapirs are in the target mantle, but no iron diapir-preexistent iron core merger has taken place yet */
	if((no_iron>0) && (no1>1) && ((timesum+timestep)>=(impact_time[no1-1]*1.000e+6*365.25*24.000*3600.000)))
                {
                merger = 1;
		/**/
		/* Make sure that correct merger condition is applied even when two collisions happen before a merger happens */
		/* Case1: One collision has occurred, no1 is already set to look for next impact, but no iron core-iron diapirs merger has taken place yet */
		if(no1==(no_merger+2)) iron_layer=d_iron;
		/**/
		/* Case2: A second collision has occurred, no1 is already set to look for next impact, but no merger has taken place yet */
		if(no1==(no_merger+3)) iron_layer=d_iron_ol;
		/**/
		/* Case3: A third collision has occurred, but no merger has taken place yet */
        	if(no1==(no_merger+4)) iron_layer=d_iron_old;
		}
	/**/
	/* Backup case 1: Iron diapir markers were found in the mantle although no merger is expected to happen as number of impacts and number of mergers agrees */
        if((no1==(no_merger+1)) && (no1>1) && (no_iron>0)) 
		{
		merger = 0;
		}
	/**/
	/* Backup case 2: A pure-silicate impactor is the next merger candidate, although no iron diapir-preexistent iron core merger can happen in this case */
	if((((no1>1) && (no1==(no_merger+2)) && (no_iron==0)) || ((no1>2) && (no1==(no_merger+3)) && (no_iron==0)) || ((no1>3) && (no1==(no_merger+4)) && (no_iron==0))) && ((timesum+timestep)>=(impact_time[no1-1]*1.000e+6*365.25*24.000*3600.000)))
		{
		printf("pure silicate merger event! \n");
		merger = 2;
		}
	/**/
	/* Backup case 3: Somehow the merger counter has reached an equal or higher value than the impact counter, reset the merger counter to the value corresponding to its currently maximum allowed value */
	if(no1<=no_merger)
		{
		no_merger = no1-1;
		merger    = 0;
		}	
	/**/
	/* Criterion for iron diapir-preexistent iron core merger: Larger part of the iron diapir markers are deeper than the preexistent iron core markers */
        if((no1>1) && (no_iron>0) && (merger==1) && (iron_layer>0.000) && (min_radius_diapir<=(mean_radius_core-0.900*(iron_layer/xsize))))
                {
	        printf("merger event! \n");
                merger = 2;
                }
	/**/
/*==============================================================================================================================================================================================*/
	/**/
	/* Save dimensional core radius [m] only when the iron core shape is unperturbed */
	if(merger==0)
		{
		/* Mean radius of the preexistent iron core [m] */
		core_radius=mean_radius_core*xsize;
		}
	/**/
/*==============================================================================================================================================================================================*/
	/**/
	/* In case no diapirs are in the target mantle OR no merger has happened between the diapirs and the preexistent iron core, compute the current mean temperature of the core [K] */
	/* In case no iron diapirs are sinking through the target mantle, keep the core temperature homogeneous */
	if(merger<=1)
		{
		/* Compute the mean temperatures of the preexistent iron core and the iron diapirs adding to the core */
		/**/
		/* Initialize counters */
		no_core  = 0;       /* Count number of preexistent iron core markers */
		/**/
		t_core   = 0.0000;  /* Used to stack preexistent core temperature AFTER removing the adiabatic term */
                t_mean   = 0.0000;  /* Used to compute mean temperature of preexistent iron core before merging */
		/**/
		/* Go through all markers */
		for(mm1=0;mm1<marknum;mm1++)
			{
			x=(double)(markx[mm1])/xsize;  /* non-dim. x coordinate */
			y=(double)(marky[mm1])/ysize;  /* non-dim. y coordinate */
			x1=x-0.50000;                  /* non-dim. x coordinate relative to planetesimal center */
			y1=(y-0.50000)*ysize/xsize;    /* non-dim. y coordinate relative to planetesimal center */
			/**/
			t_ad_core = 0.00000;           /* Set the temperature WITHOUT adiabatic contribution for every iron core marker initially to zero [K] */
			/**/
			/* Upper boundary for the iron core is the maximum iron marker radius obtained earlier */
			if((x1*x1+y1*y1)<=(uncert*mean_radius_core*uncert*mean_radius_core))
				{
				/* Lower boundary for the core is the planetesimal center */
				if((x1*x1+y1*y1)>=(0.0000*0.0000*ysize*ysize/xsize/xsize))
					{
					/* Calculate angles involved */
					ival=pow(x1*x1+y1*y1,0.500); if(!ival) ival=1.000;
					ival=x1/ival; ival=asin(ABSV(ival))/pi*180.000;
					if(x1>=0 && y1>0) ival=180.000-ival;
					if(x1<0 && y1>0) ival+=180.000;
					if(x1<0 && y1<=0) ival=360.000-ival;
					/**/
					if(ival>=ABSV(0.0000) && ival<=ABSV(360.000))
						{
						/**/
						/* Preexistent iron core material */
						if(markt[mm1]==10)
							{
							/* Compute the temperature in the preexistent iron core for the specific iron marker WITHOUT adiabatic temperature rise [K] */
							t_ad_core=markk[mm1]/exp(2.000*pi*alpha_core*rho_core*G*((max_radius_planet*max_radius_planet-(x1*x1+y1*y1))*xsize*xsize)/(3.000*cP_core));
							/**/
							t_core+=t_ad_core;   /* Stack core temperatures after removing the adiabatic contribution [K] */
							no_core+=1;          /* Count core markers involved for mean temperature calculation */
							}
						}
					}
				}
			}
		/**/
		/* Mean temperature of preexistent core without adiabatic term [K] */
		if(no_core>0) t_mean = t_core/(double)(no_core);
		/**/
		if(merger==0)
			{
			/* Keep core temperature homogenenous INCLUDING adiabatic heating */
			for(mm1=0;mm1<marknum;mm1++)
				{
				x=(double)(markx[mm1])/xsize;  /* non-dim. x coordinate */
				y=(double)(marky[mm1])/ysize;  /* non-dim. y coordinate */
				x1=x-0.50000;                  /* non-dim. x coordinate relative to planetesimal center */            
				y1=(y-0.50000)*ysize/xsize;    /* non-dim. y coordinate relative to planetesimal center */
				/**/
				/* Upper boundary for the core is the new iron core radius obtained earlier */
				if((x1*x1+y1*y1)<=(uncert*mean_radius_core*uncert*mean_radius_core))
					{
					/* Lower boundary for the core is the planetesimal center */
					if((x1*x1+y1*y1)>=(0.0000*0.0000*ysize*ysize/xsize/xsize))
						{
						/* Calculate angles involved */
						ival=pow(x1*x1+y1*y1,0.500); if(!ival) ival=1.000;
						ival=x1/ival; ival=asin(ABSV(ival))/pi*180.000;
						if(x1>=0 && y1>0) ival=180.000-ival;
						if(x1<0 && y1>0) ival+=180.000;
						if(x1<0 && y1<=0) ival=360.000-ival;
						/**/
						if(ival>=0.0000 && ival<=360.000)
							{
							/**/
							/* Only preexistent iron core markers and stray iron which already entered the core region are taken into account, all other materials are ignored */
							if(markt[mm1]==10 || markt[mm1]==17 || markt[mm1]==18 || markt[mm1]==19)
								{
                                                        	/* First reset marker temperature to zero */
                                                        	markk[mm1]=0.0000;
								/* Now put temperature of new iron core material INCLUDING adiabatic heating term on the specific marker [K] */
								markk[mm1]=t_mean*exp(2.000*pi*alpha_core*rho_core*G*((max_radius_planet*max_radius_planet-(x1*x1+y1*y1))*xsize*xsize)/(3.000*cP_core));  
								}
							}
						}
					}
				}
			}
		/* End of (merger<=1) if loop */
		}
	/**/
/*==============================================================================================================================================================================================*/
	/**/
	/* In case of diapir-preexistent core merger: Compute new mean core temperature and reset iron diapir markers */
	else if(merger==2)
		{
		/* Compute the mean temperatures of the preexistent iron core and the iron diapirs adding to the core */
		/**/
		/* Initialize counters */
		no_core       = 0;       /* Count number of preexistent iron core markers */
		no_diapir     = 0;       /* Reset to zero for reuse to count number of iron diapir markers */
		/**/
		t_core        = 0.0000;  /* Used to stack preexistent core temperature AFTER removing the adiabatic term */
		t_diapir      = 0.0000;  /* Used to stack iron diapir temperature AFTER removing the adiabatic term */
		t_mean_core   = 0.0000;  /* Mean temperature of preexistent iron core [K] */
		t_mean_diapir = 0.0000;  /* Mean temperature of iron diapirs [K] */
		t_mean        = 0.0000;  /* New mean temperature of merged iron core [K] */
		/**/
		/* Go through all markers */
		for(mm1=0;mm1<marknum;mm1++)
			{
			x=(double)(markx[mm1])/xsize;  /* non-dim. x coordinate */
			y=(double)(marky[mm1])/ysize;  /* non-dim. y coordinate */
			x1=x-0.50000;                  /* non-dim. x coordinate relative to planetesimal center */
			y1=(y-0.50000)*ysize/xsize;    /* non-dim. y coordinate relative to planetesimal center */
			/**/
			t_ad_core   = 0.00000;         /* Set the temperature WITHOUT adiabatic contribution for every iron core marker initially to zero [K] */
			t_ad_diapir = 0.00000;         /* Set the temperature WITHOUT adiabatic contribution for every iron diapir marker initially to zero [K] */
			/**/
			/* Upper boundary for the iron core is the maximum iron marker radius obtained earlier */
			if((x1*x1+y1*y1)<=(uncert*mean_radius_core*uncert*mean_radius_core))
				{
				/* Lower boundary for the core is the planetesimal center */
				if((x1*x1+y1*y1)>=(0.0000*0.0000*ysize*ysize/xsize/xsize))
					{
					/* Calculate angles involved */
					ival=pow(x1*x1+y1*y1,0.500); if(!ival) ival=1.000;
					ival=x1/ival; ival=asin(ABSV(ival))/pi*180.000;
					if(x1>=0 && y1>0) ival=180.000-ival;
					if(x1<0 && y1>0) ival+=180.000;
					if(x1<0 && y1<=0) ival=360.000-ival;
					/**/
					if(ival>=ABSV(0.0000) && ival<=ABSV(360.000))
						{
						/**/
						/* Newly added iron diapir material */
						if(markt[mm1]==(int)(iron_type[no_merger+1]))
							{
							/* Compute the temperature in the iron diapirs for the specific iron marker WITHOUT adiabatic temperature rise [K] */
							t_ad_diapir = markk[mm1]/exp(2.000*pi*alpha_core*rho_core*G*((max_radius_planet*max_radius_planet-(x1*x1+y1*y1))*xsize*xsize)/(3.000*cP_core));
							/**/
							t_diapir+=t_ad_diapir;    /* Stack diapir temperatures after removing the adiabatic contribution [K] */
							no_diapir+=1;             /* Count diapir markers involved for mean calculation */
							}
						/**/
						/* Preexistent iron core material */
						if(markt[mm1]==10)
							{
							/* Compute the temperature in the preexistent iron core for the specific iron marker WITHOUT adiabatic temperature rise [K] */
							t_ad_core = markk[mm1]/exp(2.000*pi*alpha_core*rho_core*G*((max_radius_planet*max_radius_planet-(x1*x1+y1*y1))*xsize*xsize)/(3.000*cP_core));
							/**/
							t_core+=t_ad_core;        /* Stack core temperatures after removing the adiabatic contribution [K] */
							no_core+=1;               /* Count core markers involved for mean temperature calculation */
							}
						}
					}
				}
			}
		/**/
		/* Mean temperature of preexistent core without adiabatic term [K] */
		if(no_core>0) t_mean_core = t_core/(double)(no_core);
		/**/
		/* Mean temperature of newly added iron diapir material [K] */
		if(no_diapir>0) t_mean_diapir = t_diapir/(double)(no_diapir);
		/**/
                /* Figure out which iron diapir mass has to be used for current computation */
                /* Case1: Merger of iron core with iron diapirs from latest impactor body */
		if(no1==(no_merger+2)) M_iron = M_diapir;
		/**/
		/* Case2: Merger of iron core with iron diapirs from previous impactor body */
                if(no1==(no_merger+3)) M_iron = M_diapir_ol;
		/**/
		/* Case3: Merger of iron core with iron diapirs from second to last impactor body */
		if(no1==(no_merger+4)) M_iron = M_diapir_old;
		/**/
		/* Compute mean temperatures of newly combined iron material [K], assume that both diapirs contribute in the same time to core growth */
                t_mean = (M_core*t_mean_core+M_iron*t_mean_diapir)/(M_core+M_iron);
		/**/
		/* Set all iron diapir markers inside the new mean core radius to preexistent iron core marker value */
		for(mm1=0;mm1<marknum;mm1++)
			{
			x=(double)(markx[mm1])/xsize;  /* non-dim. x coordinate */
			y=(double)(marky[mm1])/ysize;  /* non-dim. y coordinate */
			x1=x-0.50000;                  /* non-dim. x coordinate relative to planetesimal center */
			y1=(y-0.50000)*ysize/xsize;    /* non-dim. y coordinate relative to planetesimal center */
			/**/
			/* Upper boundary for the new iron core */
			if((x1*x1+y1*y1)<=((uncert*mean_radius_core)*(uncert*mean_radius_core)))
				{
				/* Lower boundary for the new iron core */
				if((x1*x1+y1*y1)>=(0.0000*0.0000*ysize*ysize/xsize/xsize))
					{
					/* Calculate angles involved */
					ival=pow(x1*x1+y1*y1,0.500); if(!ival) ival=1.000;
					ival=x1/ival; ival=asin(ABSV(ival))/pi*180.000;
					if(x1>=0 && y1>0) ival=180.000-ival;
					if(x1<0 && y1>0) ival+=180.000;
					if(x1<0 && y1<=0) ival=360.000-ival;
					/**/
					if(ival>=0.000 && ival<=360.000)
						{
						/**/
						/* All iron markers, which NEWLY contribute to the preexistent iron core are now set to be preexistent iron core markers */
						if(markt[mm1]==(int)(iron_type[no_merger+1]) || markt[mm1]==17 || markt[mm1]==18 || markt[mm1]==19)
							{
							markt[mm1]=10;
							}
						}
					}
				}
			}
		/**/
		/* All iron markers outside of the NEWLY merged iron core are assumed to be stray iron trapped in the mantle */
		for(mm1=0;mm1<marknum;mm1++)
			{
			x=(double)(markx[mm1])/xsize;  /* non-dim. x coordinate */
			y=(double)(marky[mm1])/ysize;  /* non-dim. y coordinate */
			x1=x-0.50000;                  /* non-dim. x coordinate relative to planetesimal center */
			y1=(y-0.50000)*ysize/xsize;    /* non-dim. y coordinate relative to planetesimal center */
			/**/
			/* Upper boundary for search for stray iron markers */
			/* Make sure to transform also material possibly entrained in sticky air by use of prefactor */
			if((x1*x1+y1*y1)<=((1.200*max_radius_planet)*(1.200*max_radius_planet)))
				{
				/* Lower boundary for search for stray iron markers */
				if((x1*x1+y1*y1)>(uncert*mean_radius_core*uncert*mean_radius_core*ysize*ysize/xsize/xsize))
					{
					/* Calculate angles involved */
					ival=pow(x1*x1+y1*y1,0.500); if(!ival) ival=1.000;
					ival=x1/ival; ival=asin(ABSV(ival))/pi*180.000;
					if(x1>=0 && y1>0) ival=180.000-ival;
					if(x1<0 && y1>0) ival+=180.000;
					if(x1<0 && y1<=0) ival=360.000-ival;
					/**/
					if(ival>=0.000 && ival<=360.000)
						{
						/* Iron diapir markers, which are still in the mantle of the target body */
						if(markt[mm1]==(int)(iron_type[no_merger+1]))
							{
							/* Set these markers to be stray iron markers in the mantle (physical properties identical to iron in core and standard iron diapirs) */
							markt[mm1]=(int)((int)(iron_type[no_merger+1])+10);
							}
						}
					}
				}
			}
		/**/
		/* Re-homogenize temperature of the newly merged iron core INCLUDING adiabatic heating */
		for(mm1=0;mm1<marknum;mm1++)
			{
			x=(double)(markx[mm1])/xsize;  /* non-dim. x coordinate */
			y=(double)(marky[mm1])/ysize;  /* non-dim. y coordinate */
			x1=x-0.50000;                  /* non-dim. x coordinate relative to planetesimal center */            
			y1=(y-0.50000)*ysize/xsize;    /* non-dim. y coordinate relative to planetesimal center */
			/**/
			/* Upper boundary for the core is the new iron core radius obtained earlier */
			if((x1*x1+y1*y1)<=(uncert*mean_radius_core*uncert*mean_radius_core))
				{
				/* Lower boundary for the core is the planetesimal center */
				if((x1*x1+y1*y1)>=(0.0000*0.0000*ysize*ysize/xsize/xsize))
					{
					/* Calculate angles involved */
					ival=pow(x1*x1+y1*y1,0.500); if(!ival) ival=1.000;
					ival=x1/ival; ival=asin(ABSV(ival))/pi*180.000;
					if(x1>=0 && y1>0) ival=180.000-ival;
					if(x1<0 && y1>0) ival+=180.000;
					if(x1<0 && y1<=0) ival=360.000-ival;
					/**/
					if(ival>=0.0000 && ival<=360.000)
						{
						/**/
						/* ONLY preexistent iron core markers are taken into account, all other materials are ignored */
						if(markt[mm1]==10)
							{
                                                        /* First reset marker temperature to zero */
                                                        markk[mm1]=0.0000;
							/* Now put temperature of new iron core material INCLUDING adiabatic heating term on the specific marker [K] */
							markk[mm1]=t_mean*exp(2.000*pi*alpha_core*rho_core*G*((max_radius_planet*max_radius_planet-(x1*x1+y1*y1))*xsize*xsize)/(3.000*cP_core));  
							}
						}
					}
				}
			}
                /**/
                /* The new mean radius of the preexistent iron core has to be computed [non-dim.] */
                /**/
	        /* Set counter and other parameters initially to zero */
	        ic       = 0;          /* Counter */
		ic_max   = 0;
	        g_core   = 0.000000;   /* Used to measure gravity acceleration of the upper bound of preexistent iron core [m/s^2] */
                stack    = 0.000000;   /* Used to stack maximum core radii from the different sectors [non-dim.] */
	        /**/
	        /* Make sure not to count the zero position twice */
	        for(core_angle=0.000000;core_angle<=(360.000-step);core_angle=core_angle+step)
			{
			/* Reset counter */
			ic = ic+1;
			/**/
			/* Set maximum iron core radius for each sector initially to zero */
			max_radius_core[ic] = 0.000000;
			/**/
			/* Reset distance test value for each sector */
			test_center = 0.000000;
			/**/
			/* Go through all markers */
			for(mm1=0;mm1<marknum;mm1++)
				{
				/* Reset test_angle always to zero */
				test_angle  = 0.000000;
				/**/
				/* Check only preexistent iron core markers, all other materials are excluded */
				/* Estimate radius of preexistent iron core [non-dim.] */
				if(markt[mm1]==10)
					{
					/* Angle 0 <= alpha < 90 degrees */
					if((core_angle>=0.000 && core_angle<90.000) && (markx[mm1]/xsize>=0.500) && (marky[mm1]/ysize<0.500))
						{
						test_angle = atan((markx[mm1]/xsize-0.500)/(0.500-marky[mm1]/ysize))*180.000/pi;
						if(ABSV(test_angle-core_angle)<=variation)
							{
							test_center=pow(((markx[mm1]/xsize-0.500)*(markx[mm1]/xsize-0.500)+(marky[mm1]/ysize-0.500)*(marky[mm1]/ysize-0.500)),0.500);
							/**/
							max_radius_core[ic]=MAXV(test_center,max_radius_core[ic]);
							}
						}
					/* Angle 90 <= alpha < 180 degrees */
					else if((core_angle>=90.000 && core_angle<180.000) && (markx[mm1]/xsize>0.500) && (marky[mm1]/ysize>=0.500))
						{
						test_angle = 90.000+atan((marky[mm1]/ysize-0.500)/(markx[mm1]/xsize-0.500))*180.000/pi;
						if(ABSV(test_angle-core_angle)<=variation)
							{
							test_center=pow(((markx[mm1]/xsize-0.500)*(markx[mm1]/xsize-0.500)+(marky[mm1]/ysize-0.500)*(marky[mm1]/ysize-0.500)),0.500);
							/**/
							max_radius_core[ic]=MAXV(test_center,max_radius_core[ic]);
							}
						}
					/* Angle 180 <= alpha < 270 degrees */
					else if((core_angle>=180.000 && core_angle<270.000) && (markx[mm1]/xsize<=0.500) && (marky[mm1]/ysize>0.500))
						{
						test_angle = 180.000+atan((0.500-markx[mm1]/xsize)/(marky[mm1]/ysize-0.500))*180.000/pi;
						if(ABSV(test_angle-core_angle)<=variation)
							{
							test_center=pow(((markx[mm1]/xsize-0.500)*(markx[mm1]/xsize-0.500)+(marky[mm1]/ysize-0.500)*(marky[mm1]/ysize-0.500)),0.500);
							/**/
							max_radius_core[ic]=MAXV(test_center,max_radius_core[ic]);
							}
						}
					/* Angle 270 <= alpha < 360 degrees */
					else if((core_angle>=270.000 && core_angle<360.000) && (markx[mm1]/xsize<0.500) && (marky[mm1]/ysize<=0.500))
						{
						test_angle = 270.000+atan((0.500-marky[mm1]/ysize)/(0.500-markx[mm1]/xsize))*180.000/pi;
							if(ABSV(test_angle-core_angle)<=variation)
							{
							test_center=pow(((markx[mm1]/xsize-0.500)*(markx[mm1]/xsize-0.500)+(marky[mm1]/ysize-0.500)*(marky[mm1]/ysize-0.500)),0.500);
							/**/
							max_radius_core[ic]=MAXV(test_center,max_radius_core[ic]);
							}
						}
					}
				}
			/**/
			/* Maximum number of steps */
			ic_max = ic;
			/**/
			/* Stack maximum core radii from the different sectors [non-dim.] */
			stack+=max_radius_core[ic];
			}
		/**/
		/* Compute current mean iron core radius [non-dim.] */
		mean_radius_core = stack/(double)(ic_max);
		/**/
		/* Reset mean dimensional radius of the preexistent iron core [m] */
		core_radius=mean_radius_core*xsize;
		/**/
		/* Recompute approximate gravity acceleration at CMB [m/s^2] */
		g_core = (4.000/3.000)*pi*rho_core*G*core_radius;
		/**/
                /* Reset the mass of the preexistent iron core [kg] */
                M_core+=M_iron;
		/**/
		/* Finally reset the merger counter [non-dim.] */
                no_merger+=1;
		/**/
		/* End of (merger==2) if loop */
		}
	/**/
/*==============================================================================================================================================================================================*/
	/**/
	/* Find silicate markers at the mantle side of the CMB */
	/**/
	/* Set counter and other parameters initially to zero */
	ic = 0;             /* Reset to zero for reuse */
	/**/
	/* Make sure not to count the zero position twice */
	for(core_angle=0.000000;core_angle<=(360.000-step);core_angle=core_angle+step)
		{
		ic                          = ic+1;      /* Reset counter */
		/**/
		k_mantle[ic]                = 0.000000;  /* Thermal conductivity of specific mantle marker in sector [W/(m*K)] */
		ml[ic]                      = 0;         /* Number of marker representing mantle side of CMB in sector */
		min_radius_mantle_local[ic] = 1.000000;  /* Used to check for lower bound of preexistent mantle in each sector, initially set to box size [non-dim.]  */
		t_cmb[ic]                   = 0.000000;  /* Used to measure T in each sector of the upper bound of preexistent iron core [K] */
		test_center                 = 1.000000;  /* Use for marker-planetesimal center distance measurement [non-dim.] */
		/**/
		/* Go through all markers */
		for(mm1=0;mm1<marknum;mm1++)
			{
			/* Reset the test_angle always to zero */
                        test_angle  = 0.000000;
                        /**/
			/* Check only silicate markers */
			/* Find lowermost silicate markers in planetesimal [non-dim.] */
        		if(markt[mm1]==5 || markt[mm1]==6 || markt[mm1]==25 || markt[mm1]==26)
				{
				/* Angle 0 <= alpha < 90 degrees */
                       		if((core_angle>=0.000 && core_angle<90.000) && (markx[mm1]/xsize>=0.500) && (marky[mm1]/ysize<0.500))
					{
                               		test_angle = atan((markx[mm1]/xsize-0.500)/(0.500-marky[mm1]/ysize))*180.000/pi;
                               		if(ABSV(test_angle-core_angle)<=variation)
						{
						test_center=pow(((markx[mm1]/xsize-0.500)*(markx[mm1]/xsize-0.500)+(marky[mm1]/ysize-0.500)*(marky[mm1]/ysize-0.500)),0.500);
						/**/
						/* Silicate marker potentially at mantle side of CMB sector */
						if((test_center<min_radius_mantle_local[ic]) && (test_center>max_radius_core[ic]))
							{
							if(markt[mm1]==5 || markt[mm1]==6)
								{
								alpha_mantle = markbb[5];    /* Thermal expansivity of solid silicates [1/K] */
								k_mantle[ic] = markkt[5];    /* Thermal conductivity of solid silicates [W/(m*K)] */
								}
							else if(markt[mm1]==25 || markt[mm1]==26)
								{
								alpha_mantle = markbb[25];   /* Thermal expansivity of molten silicates [1/K] */
								k_mantle[ic] = markkt[25];   /* Thermal conductivity of molten silicates [W/(m*K)] */
								}
							/**/
							/* Remove adiabatic contribution from silicate marker temperature and save temperature [K] */
							t_cmb[ic]=markk[mm1]/exp(2.000*pi*alpha_mantle*rho_mantle*G*((max_radius_planet*max_radius_planet-test_center*test_center)*xsize*xsize)/(3.000*cP_mantle));
							/**/
							/* Save marker number */
							ml[ic]    = mm1;
							}
						/**/
                       				min_radius_mantle_local[ic]=MINV(test_center,min_radius_mantle_local[ic]);
						}
					}
				/* Angle 90 <= alpha < 180 degrees */
                       		else if((core_angle>=90.000 && core_angle<180.000) && (markx[mm1]/xsize>0.500) && (marky[mm1]/ysize>=0.500))
					{
                               		test_angle = 90.000+atan((marky[mm1]/ysize-0.500)/(markx[mm1]/xsize-0.500))*180.000/pi;
                               		if(ABSV(test_angle-core_angle)<=variation)
						{
						test_center=pow(((markx[mm1]/xsize-0.500)*(markx[mm1]/xsize-0.500)+(marky[mm1]/ysize-0.500)*(marky[mm1]/ysize-0.500)),0.500);
						/**/
						/* Silicate marker potentially at mantle side of CMB sector */
						if((test_center<min_radius_mantle_local[ic]) && (test_center>max_radius_core[ic]))
							{
							if(markt[mm1]==5 || markt[mm1]==6)
								{
								alpha_mantle = markbb[5];    /* Thermal expansivity of solid silicates [1/K] */
								k_mantle[ic] = markkt[5];    /* Thermal conductivity of solid silicates [W/(m*K)] */
								}
							else if(markt[mm1]==25 || markt[mm1]==26)
								{
								alpha_mantle = markbb[25];   /* Thermal expansivity of molten silicates [1/K] */
								k_mantle[ic] = markkt[25];   /* Thermal conductivity of molten silicates [W/(m*K)] */
								}
							/**/
							/* Remove adiabatic contribution from silicate marker temperature and save temperature [K] */
							t_cmb[ic]=markk[mm1]/exp(2.000*pi*alpha_mantle*rho_mantle*G*((max_radius_planet*max_radius_planet-test_center*test_center)*xsize*xsize)/(3.000*cP_mantle));
							/**/
							/* Save marker number */
							ml[ic]    = mm1;
							}
						/**/
                       				min_radius_mantle_local[ic]=MINV(test_center,min_radius_mantle_local[ic]);
						}
					}
				/* Angle 180 <= alpha < 270 degrees */
                       		else if((core_angle>=180.000 && core_angle<270.000) && (markx[mm1]/xsize<=0.500) && (marky[mm1]/ysize>0.500))
					{
                               		test_angle = 180.000+atan((0.500-markx[mm1]/xsize)/(marky[mm1]/ysize-0.500))*180.000/pi;
                               		if(ABSV(test_angle-core_angle)<=variation)
						{
						test_center=pow(((markx[mm1]/xsize-0.500)*(markx[mm1]/xsize-0.500)+(marky[mm1]/ysize-0.500)*(marky[mm1]/ysize-0.500)),0.500);
						/**/
						/* Silicate marker potentially at mantle side of CMB sector */
						if((test_center<min_radius_mantle_local[ic]) && (test_center>max_radius_core[ic]))
							{
							if(markt[mm1]==5 || markt[mm1]==6)
								{
								alpha_mantle = markbb[5];    /* Thermal expansivity of solid silicates [1/K] */
								k_mantle[ic] = markkt[5];    /* Thermal conductivity of solid silicates [W/(m*K)] */
								}
							else if(markt[mm1]==25 || markt[mm1]==26)
								{
								alpha_mantle = markbb[25];   /* Thermal expansivity of molten silicates [1/K] */
								k_mantle[ic] = markkt[25];   /* Thermal conductivity of molten silicates [W/(m*K)] */
								}
							/**/
							/* Remove adiabatic contribution from silicate marker temperature and save temperature [K] */
							t_cmb[ic]=markk[mm1]/exp(2.000*pi*alpha_mantle*rho_mantle*G*((max_radius_planet*max_radius_planet-test_center*test_center)*xsize*xsize)/(3.000*cP_mantle));
							/**/
							/* Save marker number */
							ml[ic]   = mm1;
							}
						/**/
                       				min_radius_mantle_local[ic]=MINV(test_center,min_radius_mantle_local[ic]);
						}
					}
				/* Angle 270 <= alpha < 360 degrees */
                        	else if((core_angle>=270.000 && core_angle<360.000) && (markx[mm1]/xsize<0.500) && (marky[mm1]/ysize<=0.500))
					{
                                	test_angle = 270.000+atan((0.500-marky[mm1]/ysize)/(0.500-markx[mm1]/xsize))*180.000/pi;
                                	if(ABSV(test_angle-core_angle)<=variation)
						{
						test_center=pow(((markx[mm1]/xsize-0.500)*(markx[mm1]/xsize-0.500)+(marky[mm1]/ysize-0.500)*(marky[mm1]/ysize-0.500)),0.500);
						/**/
						/* Silicate marker potentially at mantle side of CMB sector */
						if((test_center<min_radius_mantle_local[ic]) && (test_center>max_radius_core[ic]))
							{
							if(markt[mm1]==5 || markt[mm1]==6)
								{
								alpha_mantle = markbb[5];    /* Thermal expansivity of solid silicates [1/K] */
								k_mantle[ic] = markkt[5];    /* Thermal conductivity of solid silicates [W/(m*K)] */
								}
							else if(markt[mm1]==25 || markt[mm1]==26)
								{
								alpha_mantle = markbb[25];   /* Thermal expansivity of molten silicates [1/K] */
								k_mantle[ic] = markkt[25];   /* Thermal conductivity of molten silicates [W/(m*K)] */
								}
							/**/
							/* Remove adiabatic contribution from silicate marker temperature and save temperature [K] */
							t_cmb[ic]=markk[mm1]/exp(2.000*pi*alpha_mantle*rho_mantle*G*((max_radius_planet*max_radius_planet-test_center*test_center)*xsize*xsize)/(3.000*cP_mantle));
							/**/
							/* Save marker number */
							ml[ic]   = mm1;
							}
						/**/
                       				min_radius_mantle_local[ic]=MINV(test_center,min_radius_mantle_local[ic]);
						}
					}
				}
			/* End marker loop */
			}
		/* End angle loop */
		}
	/**/
/*==============================================================================================================================================================================================*/
	/**/
	/* Compute now the temperature in the silicate layer just above the lowermost mantle [K] */
	/**/
	/* Set counter and other parameters initially to zero */
	ic = 0;                /* Reset to zero for reuse */
	/**/
	/**/
	/* Make sure not to count the zero position twice */
	for(core_angle=0.000000;core_angle<=(360.000-step);core_angle=core_angle+step)
		{
		ic            = ic+1;      /* Reset counter */
		/**/
		mll[ic]       = 0;         /* Number of silicate marker closest to silicate marker at mantle site of CMB */
		distance[ic]  = 1.000000;  /* Distance between CMB silicate marker and current silicate marker [non-dim.] */
		t_above[ic]   = 0.000000;  /* Used to measure T in each sector above the lowermost silicate layer [K] */
		test_center   = 0.000000;  /* For checking distance between planetesimal center and current marker [non-dim.] */
		test_distance = 1.000000;  /* For checking distance between CMB mantle marker and current marker [non-dim.] */
		/**/
		/* Go through all markers */
		for(mm1=0;mm1<marknum;mm1++)
			{
			/**/
			/* Check only silicate markers */
			/* Find layer ABOVE lowermost silicate markers in planetesimal [non-dim.] */
        		if(markt[mm1]==5 || markt[mm1]==6 || markt[mm1]==25 || markt[mm1]==26)
				{
				/* Angle 0 <= alpha < 90 degrees */
                       		if((core_angle>=0.000 && core_angle<90.000) && (markx[mm1]/xsize>=0.500) && (marky[mm1]/ysize<0.500))
					{
                               		test_angle = atan((markx[mm1]/xsize-0.500)/(0.500-marky[mm1]/ysize))*180.000/pi;
                               		if(ABSV(test_angle-core_angle)<=variation)
						{
						test_center=pow(((markx[mm1]/xsize-0.500)*(markx[mm1]/xsize-0.500)+(marky[mm1]/ysize-0.500)*(marky[mm1]/ysize-0.500)),0.500);
						/**/
						test_distance=pow(((markx[mm1]-markx[ml[ic]])/xsize*(markx[mm1]-markx[ml[ic]])/xsize+(marky[mm1]-marky[ml[ic]])/ysize*(marky[mm1]-marky[ml[ic]])/ysize),0.500);
						/**/
						/* Silicate marker potentially at mantle side of CMB sector, above the previous marker and not identical to it */
						if((mm1!=ml[ic]) && (test_center>min_radius_mantle_local[ic]) && (test_distance<distance[ic])) 
							{
							if(markt[mm1]==5 || markt[mm1]==6)
								{
								alpha_mantle = markbb[5];    /* Thermal expansivity of solid silicates [1/K] */
								}
							else if(markt[mm1]==25 || markt[mm1]==26)
								{
								alpha_mantle = markbb[25];   /* Thermal expansivity of molten silicates [1/K] */
								}
							/**/
							/* Remove adiabatic contribution from silicate marker temperature and save temperature [K] */
							t_above[ic]=markk[mm1]/exp(2.000*pi*alpha_mantle*rho_mantle*G*((max_radius_planet*max_radius_planet-test_center*test_center)*xsize*xsize)/(3.000*cP_mantle));
							/**/
							/* Measure distance between CMB silicate marker and current marker [non-dim.] */
                       					distance[ic]=MINV(distance[ic],test_distance);
							/**/
							/* Save marker number */
							mll[ic]      = mm1;
							}
						}
					}
				/* Angle 90 <= alpha < 180 degrees */
                       		else if((core_angle>=90.000 && core_angle<180.000) && (markx[mm1]/xsize>0.500) && (marky[mm1]/ysize>=0.500))
					{
                               		test_angle = 90.000+atan((marky[mm1]/ysize-0.500)/(markx[mm1]/xsize-0.500))*180.000/pi;
                               		if(ABSV(test_angle-core_angle)<=variation)
						{
						test_center=pow(((markx[mm1]/xsize-0.500)*(markx[mm1]/xsize-0.500)+(marky[mm1]/ysize-0.500)*(marky[mm1]/ysize-0.500)),0.500);
						/**/
						test_distance=pow(((markx[mm1]-markx[ml[ic]])/xsize*(markx[mm1]-markx[ml[ic]])/xsize+(marky[mm1]-marky[ml[ic]])/ysize*(marky[mm1]-marky[ml[ic]])/ysize),0.500);
						/**/
						/* Silicate marker potentially at mantle side of CMB sector, above the previous marker and not identical to it */
						if((mm1!=ml[ic]) && (test_center>min_radius_mantle_local[ic]) && (test_distance<distance[ic])) 
							{
							if(markt[mm1]==5 || markt[mm1]==6)
								{
								alpha_mantle = markbb[5];    /* Thermal expansivity of solid silicates [1/K] */
								}
							else if(markt[mm1]==25 || markt[mm1]==26)
								{
								alpha_mantle = markbb[25];   /* Thermal expansivity of molten silicates [1/K] */
								}
							/**/
							/* Remove adiabatic contribution from silicate marker temperature and save temperature [K] */
							t_above[ic]=markk[mm1]/exp(2.000*pi*alpha_mantle*rho_mantle*G*((max_radius_planet*max_radius_planet-test_center*test_center)*xsize*xsize)/(3.000*cP_mantle));
							/**/
							/* Measure distance between CMB silicate marker and current marker [non-dim.] */
                       					distance[ic]=MINV(distance[ic],test_distance);
							/**/
							/* Save marker number */
							mll[ic]     = mm1;
							}
						}
					}
				/* Angle 180 <= alpha < 270 degrees */
                       		else if((core_angle>=180.000 && core_angle<270.000) && (markx[mm1]/xsize<=0.500) && (marky[mm1]/ysize>0.500))
					{
                               		test_angle = 180.000+atan((0.500-markx[mm1]/xsize)/(marky[mm1]/ysize-0.500))*180.000/pi;
                               		if(ABSV(test_angle-core_angle)<=variation)
						{
						test_center=pow(((markx[mm1]/xsize-0.500)*(markx[mm1]/xsize-0.500)+(marky[mm1]/ysize-0.500)*(marky[mm1]/ysize-0.500)),0.500);
						/**/
						test_distance=pow(((markx[mm1]-markx[ml[ic]])/xsize*(markx[mm1]-markx[ml[ic]])/xsize+(marky[mm1]-marky[ml[ic]])/ysize*(marky[mm1]-marky[ml[ic]])/ysize),0.500);
						/**/
						/* Silicate marker potentially at mantle side of CMB sector, above the previous marker and not identical to it */
						if((mm1!=ml[ic]) && (test_center>min_radius_mantle_local[ic]) && (test_distance<distance[ic])) 
							{
							if(markt[mm1]==5 || markt[mm1]==6)
								{
								alpha_mantle = markbb[5];    /* Thermal expansivity of solid silicates [1/K] */
								}
							else if(markt[mm1]==25 || markt[mm1]==26)
								{
								alpha_mantle = markbb[25];   /* Thermal expansivity of molten silicates [1/K] */
								}
							/**/
							/* Remove adiabatic contribution from silicate marker temperature and save temperature [K] */
							t_above[ic]=markk[mm1]/exp(2.000*pi*alpha_mantle*rho_mantle*G*((max_radius_planet*max_radius_planet-test_center*test_center)*xsize*xsize)/(3.000*cP_mantle));
							/**/
							/* Measure distance between CMB silicate marker and current marker [non-dim.] */
                       					distance[ic]=MINV(distance[ic],test_distance);
							/**/
							/* Save marker number */
							mll[ic]     = mm1;
							}
						}
					}
				/* Angle 270 <= alpha < 360 degrees */
                        	else if((core_angle>=270.000 && core_angle<360.000) && (markx[mm1]/xsize<0.500) && (marky[mm1]/ysize<=0.500))
					{
                                	test_angle = 270.000+atan((0.500-marky[mm1]/ysize)/(0.500-markx[mm1]/xsize))*180.000/pi;
                                	if(ABSV(test_angle-core_angle)<=variation)
						{
						test_center=pow(((markx[mm1]/xsize-0.500)*(markx[mm1]/xsize-0.500)+(marky[mm1]/ysize-0.500)*(marky[mm1]/ysize-0.500)),0.500);
						/**/
						test_distance=pow(((markx[mm1]-markx[ml[ic]])/xsize*(markx[mm1]-markx[ml[ic]])/xsize+(marky[mm1]-marky[ml[ic]])/ysize*(marky[mm1]-marky[ml[ic]])/ysize),0.500);
						/**/
						/* Silicate marker potentially at mantle side of CMB sector, above the previous marker and not identical to it */
						if((mm1!=ml[ic]) && (test_center>min_radius_mantle_local[ic]) && (test_distance<distance[ic])) 
							{
							if(markt[mm1]==5 || markt[mm1]==6)
								{
								alpha_mantle = markbb[5];    /* Thermal expansivity of solid silicates [1/K] */
								}
							else if(markt[mm1]==25 || markt[mm1]==26)
								{
								alpha_mantle = markbb[25];   /* Thermal expansivity of molten silicates [1/K] */
								}
							/**/
							/* Remove adiabatic contribution from silicate marker temperature and save temperature [K] */
							t_above[ic]=markk[mm1]/exp(2.000*pi*alpha_mantle*rho_mantle*G*((max_radius_planet*max_radius_planet-test_center*test_center)*xsize*xsize)/(3.000*cP_mantle));
							/**/
							/* Measure distance between CMB silicate marker and current marker [non-dim.] */
                       					distance[ic]=MINV(distance[ic],test_distance);
							/**/
							/* Save marker number */
							mll[ic]     = mm1;
							}
						}
					}
				}
			/* End marker loop */
			}
		/* End angle loop */
		}
	/**/
/*==============================================================================================================================================================================================*/
	/**/
	/* Compute heat flux through lowermost mantle [W/m^2] */
	/**/
	/* Reset all general parameters to zero before each new use */
	q_adv_stack    = 0.00000;
	q_cond_stack   = 0.00000;
	q_c            = 0.00000;
	q_a            = 0.00000;
	t_cmb_stack    = 0.00000;
	t_cmb_mean     = 0.00000;
	/**/
	/**/
	/* Compute now the individual conductive heat flux through the lowermost silicate mantle [W/m^2] in each sector after adiabatic effects are removed */
	for(ic=1;ic<=ic_max;ic++)
		{
		/* Reset all local parameters to zero before each new use */
		x_dist     = 0.00000;
		y_dist     = 0.00000;
		x_nor      = 0.00000;
		y_nor      = 0.00000;
		l_nor      = 0.00000;
		x_nor1[ic] = 0.00000;
		y_nor1[ic] = 0.00000;
		q_cond     = 0.00000;
		zeta       = 0.00000;
		/**/
		/**/
		/* Compute x and y component of vector normal to iron core surface [non-dim.] */
		x_nor=((double)(markx[ml[ic]])/xsize)-0.5000;
		y_nor=((double)(marky[ml[ic]])/ysize)-0.5000;
		/**/
		/* Compute length of normal vector [non-dim.] */
		l_nor = pow((x_nor*x_nor+y_nor*y_nor),0.5000);
		/**/
		/* Now compute components for unity vector standing normal to the iron core surface [non-dim.] */
		x_nor1[ic] = x_nor/l_nor;
		y_nor1[ic] = y_nor/l_nor;
		/**/
		/* For completeness compute exact length of the unity vector [non-dim.] */
		l_nor1 = pow((x_nor1[ic]*x_nor1[ic]+y_nor1[ic]*y_nor1[ic]),0.5000);
		/**/
		/**/
		/* Compute vector between lowermost and overlying silicate marker [non-dim.] */
		x_dist=((double)(markx[mll[ic]])-(double)(markx[ml[ic]]))/xsize;
		y_dist=((double)(marky[mll[ic]])-(double)(marky[ml[ic]]))/ysize;
		/**/
		/* Compute angle between both vectors [degrees] */
		zeta = acos((x_nor1[ic]*x_dist+y_nor1[ic]*y_dist)/(l_nor1*distance[ic]))*(180.0000/pi);
		/**/
		/* Compute the conductive heat flux [W/m^2], consider only the radial component of the distance vector! */
		q_cond = k_mantle[ic]*(t_cmb[ic]-t_above[ic])/(distance[ic]*cos(zeta*pi/180.0000)*xsize);
		/**/
		/* Stack the conductive and advective heat flux from the lowermost silicate mantle [W/m^2] */
		q_cond_stack+=q_cond;
		t_cmb_stack+=t_cmb[ic];
		}
	/**/
	/* Current mean temperature of mantle side of CMB [K] */
	t_cmb_mean = t_cmb_stack/(double)(ic_max);
	/**/
	/**/
	/* Compute advective heat flux [W/m^2] in each sector after adiabatic effects are removed */
	for(ic=1;ic<=ic_max;ic++)
		{
		/* Reset all local parameters to zero before each new use */
		vx_marker = 0.00000;
		vy_marker = 0.00000;
		vx_rad    = 0.00000;
		vy_rad    = 0.00000;
		v_rad     = 0.00000;
		q_adv     = 0.00000;
		/**/
		/**/
		/* Now interpolate current x and y velocity components [m/s] of the specific silicate marker */
		allinterv((double)(markx[ml[ic]]),(double)(marky[ml[ic]]));
		/**/
		/* Save the x and y velocity components of the specific silicate marker [m/s] */
		vx_marker=eps[11];
		vy_marker=eps[12];
		/**/
		/* Now compute components of the velocity vector projected on the unity vector standing normal to the iron core surface [m/s] */
		/* based on equation: v_projected = (v_vec*e_vec)*e_vec */
		/* Use components of normal unity vector computed in the previous loop */
		vx_rad = (vx_marker*x_nor1[ic]+vy_marker*y_nor1[ic])*x_nor1[ic];
		vy_rad = (vx_marker*x_nor1[ic]+vy_marker*y_nor1[ic])*y_nor1[ic];	
		/**/
		/* Now compute the corresponding total velocity [m/s] */
		v_rad  = pow((vx_rad*vx_rad+vy_rad*vy_rad),0.5000);
		/**/
		/* Finally compute advective heat flux through the lowermost silicate mantle [W/m^2] after adiabatic effects and the mean CMB silicate temperature are removed */
		q_adv  = rho_mantle*cP_mantle*v_rad*(t_cmb[ic]-t_cmb_mean);
		/**/
		/* Stack the conductive and advective heat flux from the lowermost silicate mantle [W/m^2] */
		q_adv_stack+=q_adv;
		}
	/**/
	/* Now compute current mean conductive and advective heat flux through the bottom of the silicate mantle [W/m^2] */
	q_c = q_cond_stack/(double)(ic_max);
	q_a = q_adv_stack/(double)(ic_max);
	/**/
	/* Compute current combined mean heat flux through the lowermost silicate mantle [W/m^2] */
	q_mantle = q_c+q_a;
	/**/
	/**/
	/* Now compute the minimum conductive heat flux across the CMB needed to drive thermal convection in the iron core [W/m^2] */
	/* (D.J. Stevenson, EPSL, 208, 1-11, 2003a) */
	q_core = k_core*alpha_core*g_core*t_mean/cP_core;
	/**/
        /**/
/*==============================================================================================================================================================================================*/
	/**/
	/* Is the central region of the preexistent iron core molten? */
	pr_center      = 0.00000;
	t_melting      = 0.00000;
	/**/
	/* Node number of center of planetary body */
        m3=(((xnumx-1)/2)+1)*ynumy+(((ynumy-1)/2)+1);
	/**/
	/* Current pressure at the center of the planetary embryo [bar] */
	pr_center=pr[m3]/1.000e5;
        /**/
        /* Compute melting temperature of pure iron for center of planetary body [K] */
	if(fe_melting == 0)
		{
		/* Pure iron melting curve parametrized after (Boehler, Nature, 363, 534-536, 1993) */
		/**/
		if(pr_center>=0.000e6 && pr_center<0.100e6)            /*  0 GPa <= P < 10 GPa */
			{
			t_melting=1761.000+3.100e-3*pr_center;
			}
		if(pr_center>=0.100e6 && pr_center<0.200e6)            /* 10 GPa <= P < 20 GPa */
			{
			t_melting=1863.000+2.080e-3*pr_center;
			}
		if(pr_center>=0.200e6 && pr_center<0.600e6)            /* 20 GPa <= P < 60 GPa */
			{
			t_melting=2071.800+1.035e-3*pr_center;
			}
		if(pr_center>=0.600e6 && pr_center<1.000e6)            /* 60 GPa <= P < 100 GPa */
			{
			t_melting=2382.800+5.170e-4*pr_center;
			}
		if(pr_center>=1.000e6)                                 /* P >= 100 GPa */
			{
			t_melting=2900.000+(1.000/1.000e3)*pr_center;
			}
		}
	/**/
        /* Compute melting temperature of eutectic iron-sulfide for center of planetary body [K] */
	else if(fe_melting == 1)
		{
		/* Eutectic Fe-FeS melting curve using data compilation by (Chudinovskikh & Boehler, EPSL, 257, 97-103, 2007) */
		/**/
		if(pr_center<=0.4180e6)        /* P < (0.418e6 bar = 4.18e5 bar = 4.18e10 Pa =) 41.8 GPa */
			{
			t_melting=1260.1000+3.3171e-4*pr_center-3.395e-9*pr_center*pr_center+2.660e-14*pr_center*pr_center*pr_center-3.7688e-20*pow(pr_center,4.000);
			}
		else
		/* Linear extrapolation to higher pressures */
			{
			t_melting=1597.7000+4.2635e-4*(pr_center-0.4180e6);
			}
		}
	/**/
/*==============================================================================================================================================================================================*/
	/**/
	/* Is magnetic dynamo activity possible? */
	/**/
	/* - heat flux through lowermost mantle exceeds heat flux out of the iron core that would be carried by conduction along the core adiabat */
	/* - mean temperature of core is higher than the mean temperature of the mantle side of the CMB, thus the core is cooling (radiogenic heating could blanket iron core) */
        /* - temperature at center of planetary body is above melting temperature of pure iron or iron-sulfide eutectic */
	/**/
	if((q_mantle>q_core) && (t_mean>t_cmb_mean) && (t_mean>t_melting))
		{
		/* Compute magnetic Rayleigh number [non-dim.] (Aubert et al., GJI, 179, 1414-1428, 2009, eq. (12)) */
		/* Multiply mean heat flux difference [W/m^2] with surface area of core [m^2] and term (alpha_core/cP_core) to obtain proper unit [kg/s] */
		ra_q     = g_core*(q_mantle-q_core)*4.000*pi*pow(core_radius,2.000)*(alpha_core/cP_core)/(4.000*pi*rho_core*pow((2.000*pi*spin_rate[no1]/year_to_sec),3.000)*pow(core_radius,4.000));
		/**/
		/* Now compute power per unit volume [non-dim.] (Aubert et al., GJI, 179, 1414-1428, 2009, eq. (17)+(18)) */
		/* factor gamma computed using following assumption that no solid inner core exists: f_i=0, r_i=0, xi=0 */
		ppu      = 0.6000*ra_q;
		/**/
		/* Compute Ekman number [non-dim.] (Aubert et al., GJI, 179, 1414-1428, 2009, eq. (7)) */
		ek       = eta_fe/(rho_core*(2.000*pi*spin_rate[no1]/year_to_sec)*pow(core_radius,2.000));
		/**/
		/* Compute magnetic Prandtl number [non-dim.] (Aubert et al., GJI, 179, 1414-1428, 2009, eq. (8)) */
		prm      = eta_fe/(rho_core*magn_diff);
		/**/
		/* Magnetic Reynolds number (based on Aubert et al., GJI, 179, 1414-1428, 2009, eq. (19) and line above eq. (21)) */
		reynolds = 1.310*pow(ppu,0.420)*prm/ek;
		}
	/**/
	/**/
	/* Use following criteria to check whether the iron core is convecting and a magnetic field can be generated: */
	/* - heat flux through lowermost mantle exceeds heat flux out of the iron core that would be carried by conduction along the core adiabat */
	/* - mean temperature of core is higher than the mean temperature of the mantle side of the CMB, thus the core is cooling (radiogenic heating could blanket iron core) */
        /* - temperature at center of planetary body is above melting temperature of pure iron or iron-sulfide eutectic */
	/* - magnetic Reynolds number is larger than 40 (Christensen & Aubert, GJI, 117, 97-114, 2006) */
	/**/
	if((q_mantle>q_core) && (t_mean>t_cmb_mean) && (t_mean>t_melting) && (reynolds>40.000))
                {
                dynamo = 1;
                }
	/**/
	/* Compute virtual dipole moment [A*m^2] when planetary dynamo is active using spin rate data from N-body simulation (based on Aubert et al., GJI, 179, 1414-1428, 2009) */
	if(dynamo==1)
		{
		/* Magnetic field amplitude inside the shell [T = N/(A*m)] (Aubert et al., GJI, 179, 1414-1428, 2009, eq. (22)) */
		/* c1 taken as mean value [non-dim.] (Aubert et al., GJI, 179, 1414-1428, 2009, eq. (21)) */
		b_rms    = 1.170*pow(ppu,0.340)*pow((rho_core*magn_perm),0.500)*(2.000*pi*spin_rate[no1]/year_to_sec)*core_radius;
		/**/
		/* Ratio of mean strength inside the shell to the dipole strength on the outer boundary [non-dim.] (Aubert et al., GJI, 179, 1414-1428, 2009, eq. (28)) */
		b_dip    = 7.300;
		/**/
		/* Compute magnetic dipole moment [A*m^2] (Aubert et al., GJI, 179, 1414-1428, 2009, eq. (38)) */
		dipl_mom = 4.000*pi*pow(core_radius,3.000)*b_rms/(pow(2.000,0.500)*magn_perm*b_dip);
		/**/
		/**/
		/* Now compute parameters needed to estimate possibility of magnetic reversals */
		/**/
		/* Compute Prandtl number [non-dim.] (Aubert et al., GJI, 179, 1414-1428, 2009, eq. (9)) */
		prandtl  = eta_fe/(rho_core*diff_fe);
		/**/
		/* Compute local Rossby number (dimensionless rms flow velocity in iron core) (Aubert et al., GJI, 179, 1414-1428, 2009, eq. (31)) */
		/* Conservative estimation: Reversals when rossby > 0.10 (Aubert et al., GJI, 179, 1414-1428, 2009, eq. (31)) */
		rossby   = 0.540*pow(ppu,0.480)*pow(ek,(-0.320))*pow(prandtl,0.190)*pow(prm,(-0.190));
		}
	/**/
	/* End of time>0.000000 loop */
	}
	/**/
	return 0;
	/**/
}
/* End subroutine core */
