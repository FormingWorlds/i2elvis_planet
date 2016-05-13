/* Compute olivine grain growth in pallasites based on [Solferino et al., Geochim. Cosmochim. Acta, 162, 259-275, 2015] */
/* written by Greg (last updated: 13/05/2016) */
int grain()
{
	long int mm1;
	/**/
	/* Old olivine grain growth parameters [Solferino et al., Geochim. Cosmochim. Acta, 2015] */
	/* double n_grain=2.415; */          /* Olivine grain growth exponent [non-dim.] */
	/* double E_grain=2.88683e+5; */     /* Activation energy [J/mol] */
	/* double k0_grain=9.434e+6; */      /* Olivine grain growth constant [mum^n/s] */
	/**/
	/* New olivine grain growth parameters */
	double n_grain=3.70;         /* Olivine grain growth exponent [non-dim.] */
	double E_grain=1.0098e+5;    /* Activation energy [J/mol] */
	double k0_grain=3.2033;      /* Olivine grain growth constant [mum^n/s] */
	/**/
	double pall_gr1,pall_gr2,pall_gr3;
	/**/
	/* Check all markers */
	for (mm1=0;mm1<marknum;mm1++)
		{
		/* Olivine grain growth only activated after impact event */
		if((timesum+timestep)>=(impact_time[1]*1.000e+6*365.250*24.000*3600.000))
			{
			/* Olivine grains can only grow efficiently as long as (i) temperature is below olivine liquidus and (ii) Fe-Ni-S is still liquid [Solferino et al., Geochim. Cosmochim. Acta, 2015] */
			/* Olivine is solid, but Fe-Ni-S melting temperature is reached [Clark and Kullerud, 1963; Waldner and Pelton, 2004] */
			if((markt[mm1]==5 || markt[mm1]==6) && markk[mm1]>=1073.000)
				{
				pall_gr1    = pow(((double)(markgr[mm1])),n_grain);
				pall_gr2    = k0_grain*exp(-E_grain/(8.3145*markk[mm1]))*timestep;
				pall_gr3    = pow((pall_gr1+pall_gr2),(1.000/n_grain));
				markgr[mm1] = (float)(pall_gr3);
				}
			/* Partial silicate melt, but still beneath olivine liquidus */
			if((markt[mm1]==25 || markt[mm1]==26) && markk[mm1]<2103.000)
				{
				pall_gr1    = pow(((double)(markgr[mm1])),n_grain);
				pall_gr2    = k0_grain*exp(-E_grain/(8.3145*markk[mm1]))*timestep;
				pall_gr3    = pow((pall_gr1+pall_gr2),(1.000/n_grain));
				markgr[mm1] = (float)(pall_gr3);				
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
