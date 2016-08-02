/* Compute pure olivine grain growth [Karato, 1989] and growth of olivine in partially/fully molten FeS in pallasites [Solferino et al., Geochim. Cosmochim. Acta, 162, 259-275, 2015] */
/* written by Gregor J. Golabek (last updated: 01/08/2016) */
int grain()
{
	long int mm1;
	/**/
	/* Old olivine grain growth parameters [Solferino et al., Geochim. Cosmochim. Acta, 2015] */
	/* double n_grain=2.415; */          /* Olivine grain growth exponent [non-dim.] */
	/* double E_grain=2.88683e+5; */     /* Activation energy [J/mol] */
	/* double k0_grain=9.434e+6; */      /* Olivine grain growth constant [mum^n/s] */
	/**/
	/* New olivine in partially molten FeS grain growth parameters */
	double n_grain=3.70;           /* Olivine in partially molten FeS grain growth exponent [non-dim.] */
	double E_grain=1.0098e+5;      /* Activation energy [J/mol] */
	double k0_grain=3.2033;        /* Olivine in partially molten FeS grain growth constant [mum^n/s] */
	/**/
	double n_oliv=2.250;           /* Pure olivine grain growth exponent [non-dim.] */
	double E_oliv=5.71674e+5;      /* Activation energy [J/mol] */
	double k0_oliv=1.06593e+16;    /* Pure olivine grain growth exponent [mum^n/s] */
	/**/
	double pival=3.141592654;      /* Define Pi */
	/**/
	double radius_core=1.000e+5;   /* Core radius [m] */
	double radius_planet=2.000e+5; /* Planetesimal radius [m] */
	double ref_angle=45.000;       /* Angle of region of interest [degrees] */
	double dev_angle=10.000;       /* Width of region of interest [degrees] */
	/**/
	double curr_angle,distancee,dist_x,dist_y;
	double oliv_gr1,oliv_gr2,oliv_gr3,pall_gr1,pall_gr2,pall_gr3;
	/**/
	/* Check all markers */
	for (mm1=0;mm1<marknum;mm1++)
		{
		/* Reset marker type in sector of circle after first impact event */
        	if((timesum+timestep)>=(impact_time[1]*1.000e+6*365.250*24.000*3600.000) && (timesum+timestep)<((impact_time[1]+0.014)*1.000e+6*365.250*24.000*3600.000))
			{
			/* Relative coordinates [m] */
			dist_x    = markx[mm1]-xsize/2.000;
			dist_y    = marky[mm1]-ysize/2.000;
			/**/
			/* Compute distance from planetesimal center [m] */
			distancee = pow((dist_x*dist_x+dist_y*dist_y),0.500);
			/**/
			/* Only markers inside the mantle will be reset */
			if(distancee>=radius_core && distancee<=radius_planet)
				{
				/* Compute angle of marker relative to center of planetesimal [degrees] */
				curr_angle = asin(dist_x/distancee)/pival*180.000;
				/**/
				/* Reset marker type in sector of interest */
				if(curr_angle>=ref_angle-dev_angle && curr_angle<=ref_angle+dev_angle && dist_x>0.000 && dist_y<0.000)
					{
					if(markt[mm1]==6)
						{
						markt[mm1]=5;
						}
					if(markt[mm1]==26)
						{
						markt[mm1]=25;
						}
					}
				}
			}
		/**/
		/* Olivine grain growth only activated after impact event */
		if((timesum+timestep)>=(impact_time[1]*1.000e+6*365.250*24.000*3600.000))
			{
			/* Olivine grains can only grow efficiently as long as (i) temperature is below olivine liquidus and (ii) Fe-Ni-S is still liquid [Solferino et al., Geochim. Cosmochim. Acta, 2015] */
			/* Olivine is solid, but Fe-Ni-S melting temperature is reached [Clark and Kullerud, 1963; Waldner and Pelton, 2004] */
			if(markt[mm1]==5 && markk[mm1]>=1073.000)
				{
				pall_gr1    = pow(((double)(markgr[mm1])),n_grain);
				pall_gr2    = k0_grain*exp(-E_grain/(8.3145*markk[mm1]))*timestep;
				pall_gr3    = pow((pall_gr1+pall_gr2),(1.000/n_grain));
				markgr[mm1] = (float)(pall_gr3);
				}
			/* Partial silicate melt, but still beneath olivine liquidus */
			if(markt[mm1]==25 && markk[mm1]<2103.000)
				{
				pall_gr1    = pow(((double)(markgr[mm1])),n_grain);
				pall_gr2    = k0_grain*exp(-E_grain/(8.3145*markk[mm1]))*timestep;
				pall_gr3    = pow((pall_gr1+pall_gr2),(1.000/n_grain));
				markgr[mm1] = (float)(pall_gr3);				
				}
			/* Growth of pure solid olivine */
			if(markt[mm1]==6)
				{
				oliv_gr1    = pow(((double)(markgr[mm1])),n_oliv);
				oliv_gr2    = k0_oliv*exp(-E_oliv/(8.3145*markk[mm1]))*timestep;
				oliv_gr3    = pow((oliv_gr1+oliv_gr2),(1.000/n_oliv));
				markgr[mm1] = (float)(oliv_gr3);
				}
			/* Growth of pure olivine only possible as temperature is below olivine liquidus */
			if(markt[mm1]==26 && markk[mm1]<2103.000)
				{
				oliv_gr1    = pow(((double)(markgr[mm1])),n_oliv);
				oliv_gr2    = k0_oliv*exp(-E_oliv/(8.3145*markk[mm1]))*timestep;
				oliv_gr3    = pow((oliv_gr1+oliv_gr2),(1.000/n_oliv));
				markgr[mm1] = (float)(oliv_gr3);
				}
			/**/
			/* Grain size is reset to start value in case temperature rises above olivine liquidus */
			if((markt[mm1]==25 || markt[mm1]==26) && markk[mm1]>=2103.000)
				{
				markgr[mm1] = (float)(gr_init);	
				}
			}
		}
	/**/
	return 0;
/* End grain growth subroutine */
}
