/* Load information from configuration file mode.t3c ============== */
int loadconf()
{
/* Counter */
int n1,n2,n3,fln3=0;
double ival;
/**/
/**/
/**/
/* Open File file.t3c */
fl = fopen("file.t3c","rt");
ffscanf(); fln3=atoi(sa)-1;
fclose(fl);
/**/
/**/
/**/
/* Open File mode.t3c */
fl = fopen("mode.t3c","rt");
/**/
/* Data File name */
ffscanf();
for (n1=0;n1<50;n1++) fl1in[n1]=sa[n1];
ffscanf(); if(sa[0] == 'b') fl1itp=1;
/**/
/* Load first Results File names */
ffscanf();
fl0num=0;
while(sa[0]!='~')
	{
	/* Check file Counter */
	if(fl0num>=MAXFLN) {printf("Space out in fl0out[]"); exit(0);}
	/**/
	/* Save results file name */
	for (n1=0;n1<50;n1++) fl0out[fl0num][n1]=sa[n1];
	/**/
	/* Load TYPE_cyc0max_maxxystep_maxtkstep_maxtmstep */
	ffscanf(); if(sa[0] == 'b') fl0otp[fl0num]=1;
	ffscanf();fl0cyc[fl0num]=atoi(sa);
	ffscanf();fl0stp[fl0num][0]=atof(sa);
	ffscanf();fl0stp[fl0num][1]=atof(sa);
	ffscanf();fl0stp[fl0num][2]=atof(sa)*3.15576e+7;
	/**/
	/* Incr File Counters */
	fl0num++;
	/**/
	/* Load Next Results File names */
	ffscanf();
	}
/**/
/* Service */
ffscanf();printmod=atoi(sa);
ffscanf();timedir=atof(sa);
ffscanf();movemod=atoi(sa);
ffscanf();tempmod=atoi(sa);
ffscanf();markmod=atoi(sa);
ffscanf();ratemod=atoi(sa);
ffscanf();gridmod=atoi(sa);
ffscanf();intermod=atoi(sa);
ffscanf();intermod1=atoi(sa);
ffscanf();outgrid=atoi(sa);
ffscanf();densimod=atoi(sa);
/**/
/* V */
ffscanf();DIVVMIN=atof(sa);
ffscanf();STOKSMIN=atof(sa);
ffscanf();stoksmod=atoi(sa);
ffscanf();presmod=atoi(sa);
ffscanf();stoksfd=atoi(sa);
ffscanf();nubeg=atof(sa);
ffscanf();nuend=atof(sa);
ffscanf();nucontr=atof(sa);
ffscanf();hidry=atof(sa);
ffscanf();hidrl=atof(sa);
ffscanf();strmin=atof(sa);
ffscanf();strmax=atof(sa);
ffscanf();viscmod=atoi(sa);
/**/
/**/
/**/
/* T */
ffscanf();HEATMIN=atof(sa);
ffscanf();heatmod=atoi(sa);
ffscanf();heatfd=atoi(sa);
ffscanf();heatdif=atof(sa);
ffscanf();frictyn=atof(sa)*3.15576e+7;
ffscanf();adiabyn=atof(sa)*3.15576e+7;
ffscanf();core_form_time=atof(sa)*3.15576e+7;
ffscanf();si_melting=atoi(sa);
ffscanf();fe_melting=atoi(sa);
/**/
/**/
/**/
/* Data File name change after number */
if(fln3>=0 && fln3<fl0num)
	{
	for (n1=0;n1<50;n1++) fl1in[n1]=fl0out[fln3][n1];
	fl1itp=fl0otp[fln3];
	}
else
	{
	fln3=-1;
	}
/**/
fclose(fl);
/* End Load information from configuration file mode.t3c */
/**/
/* stop.yn file creation */
fl = fopen("stop.yn","wt");
fprintf(fl,"n \n");
fclose(fl);
/**/
/**/
/**/
/* Load thermodynamic database */
/* Database P-T grid parameters */
tknum=350;
pbnum=350;
tkmin=500.005;
pbmin=1000.5;
tkstp=(5000.0-tkmin)/(double)(tknum-1);
pbstp=(500000.0-pbmin)/(double)(pbnum-1);
/*
printf("%d %d %e %e %e %e",tknum,pbnum,tkmin,pbmin,tkstp,pbstp);getchar();
*/
/**/
/* Mars Mantle */
/* RO - density */
fl = fopen("mars_mantle","rt");
/* Read heading line */
for (n1=0;n1<46;n1++)
	{
	fscanf(fl,"%s",sa);
/*
printf("A   %d  %s",n1,sa);getchar();
*/
	}
/* Read Database %
/**/
for (n1=0;n1<pbnum;n1++)
for (n2=0;n2<tknum;n2++)
	{
	for (n3=0;n3<42;n3++)
		{
		fscanf(fl,"%s",sa);
		wi[n3]=atof(sa);
/*
printf("%d %e ",n3,wi[n3]);getchar();
*/
		}
	/* Density, g/cm3 */
	td[n2][n1][0][0]=wi[11]/1000.0;
	/* Enthalpy, kCal/kg */
	td[n2][n1][0][1]=wi[3]/wi[18]/4.1837;
/*
td[n2][n1][0][1]=wi[3]/wi[18]*1e+3;
if (n2>0)
{
printf("%d %d %e %e %e ",n1,n2,td[n2-1][n1][0][1],td[n2][n1][0][1],(td[n2][n1][0][1]-td[n2-1][n1][0][1])/tkstp);getchar();
}
printf("%d %d %e %e %e ",n1,n2,td[n2][n1][0][0],td[n2][n1][0][1],td[n2][n1][0][2]);getchar();
	td[n2][n1][0][1]=wi[3]/wi[2]*1e+5/wi[11];
	td[n2][n1][0][2]=wi[3]/wi[18]*1e+3;
*/
	}
fclose(fl);
printf("mars_mantle OK \n");
/**/
/* Mars Crust */
/* RO - density */
fl = fopen("mars_crust","rt");
/* Read heading line */
for (n1=0;n1<48;n1++)
	{
	fscanf(fl,"%s",sa);
/*
printf("A   %d  %s",n1,sa);getchar();
*/
	}
/* Read Database %
/**/
for (n1=0;n1<pbnum;n1++)
for (n2=0;n2<tknum;n2++)
	{
	for (n3=0;n3<44;n3++)
		{
		fscanf(fl,"%s",sa);
		wi[n3]=atof(sa);
/*
printf("%d %e ",n3,wi[n3]);getchar();
*/
		}
	/* Density, g/cm3 */
	td[n2][n1][1][0]=wi[11]/1000.0;
	/* Enthalpy, kCal/kg */
	td[n2][n1][1][1]=wi[3]/wi[18]/4.1837;
/*
td[n2][n1][0][1]=wi[3]/wi[18]*1e+3;
if (n2>0)
{
printf("%d %d %e %e %e ",n1,n2,td[n2-1][n1][0][1],td[n2][n1][0][1],(td[n2][n1][0][1]-td[n2-1][n1][0][1])/tkstp);getchar();
}
printf("%d %d %e %e %e ",n1,n2,td[n2][n1][0][0],td[n2][n1][0][1],td[n2][n1][0][2]);getchar();
td[n2][n1][0][1]=wi[3]/wi[2]*1e+5/wi[11];
td[n2][n1][0][2]=wi[3]/wi[18]*1e+3;
*/
	}
fclose(fl);
printf("mars_crust OK \n");
/**/
return fln3;
}
/* Load information from configuration file mode.t3c ============== */



/* Load Information from data file ------------------------------- */
void loader()
/* bondv[] - bondary value */
/* bondm[] - bondary mode 0=Not, -1=Value, 1,2...=LinNum+1 */
/* m1,m2 - node X,Y number */
{
/* Counter */
int n1;
char nn1,nn2,nn3;
long int m1,m2,m3;
long int mm1;
char szint,szlong,szfloat,szdouble,szcur;
float ival0;
double ival1;
/**/
/**/
/**/
/* Load Past Results from data file-------------------------------- */
if (printmod) printf("Load Past results from %s ...",fl1in);
/**/
/**/
/**/
/* Load in Text Format ---------------------------- */
if(fl1itp==0)
	{
	fl = fopen(fl1in,"rt");
	/**/
	/* Grid Parameters */
	ffscanf();xnumx=atoi(sa);
	ffscanf();ynumy=atoi(sa);
	ffscanf();mnumx=atoi(sa);
	ffscanf();mnumy=atoi(sa);
	ffscanf();marknum=atoi(sa);
	ffscanf();xsize=atof(sa);
	ffscanf();ysize=atof(sa);
	ffscanf();gamma_eff=atof(sa);
	ffscanf();memory_fe=atof(sa);
	ffscanf();memory_si=atof(sa);
	ffscanf();por_init=atof(sa);
	ffscanf();growth_model=atoi(sa);
	ffscanf();gr_init=atof(sa);
	ffscanf();znumz=atoi(sa);
	ffscanf();corr2d3d=atoi(sa);
	ffscanf();pinit=atof(sa);
	ffscanf();pkf[0]=atof(sa);
	ffscanf();pkf[1]=atof(sa);
	ffscanf();pkf[2]=atof(sa);
	ffscanf();pkf[3]=atof(sa);
	ffscanf();GXKOEF=atof(sa);
	ffscanf();GYKOEF=atof(sa);
	ffscanf();tmp_ambient=atof(sa);
	ffscanf();delta_tmp=atof(sa);
	ffscanf();timeexit=atof(sa)*3.15576e+7;
	ffscanf();al2627_init=atof(sa)*1.0e-5;
        ffscanf();fe6056_init=atof(sa)*1.0e-8;
	ffscanf();rocknum=atoi(sa);
	ffscanf();bondnum=atoi(sa);
	ffscanf();
	ffscanf();timesum=atof(sa)*3.15576e+7;
	/**/
	/* Calc,Check Grid parameters */
	gridcheck();
	/**/
	/* Rock Types information */
	for (n1=0;n1<rocknum;n1++)
		{
		ffscanf();
		ffscanf();markim[n1]=atoi(sa);;
		ffscanf();markn0[n1]=atof(sa);;
		ffscanf();markn1[n1]=atof(sa);;
		ffscanf();marks0[n1]=atof(sa);;
		ffscanf();marks1[n1]=atof(sa);;
		ffscanf();marknu[n1]=atof(sa);;
		ffscanf();markdh[n1]=atof(sa);;
		ffscanf();markdv[n1]=atof(sa);;
		ffscanf();markss[n1]=atof(sa);;
		ffscanf();markmm[n1]=atof(sa);;
		ffscanf();markll[n1]=atof(sa);;
		ffscanf();marka0[n1]=atof(sa);;
		ffscanf();marka1[n1]=atof(sa);;
		ffscanf();markb0[n1]=atof(sa);;
		ffscanf();markb1[n1]=atof(sa);;
		ffscanf();marke0[n1]=atof(sa);;
		ffscanf();marke1[n1]=atof(sa);;
		ffscanf();markro[n1]=atof(sa);;
		ffscanf();markbb[n1]=atof(sa);;
		ffscanf();markaa[n1]=atof(sa);;
		ffscanf();markcp[n1]=atof(sa);;
		ffscanf();markkt[n1]=atof(sa);;
		ffscanf();markkf[n1]=atof(sa);;
		ffscanf();markkp[n1]=atof(sa);;
		ffscanf();markht[n1]=atof(sa);;
		}
	/**/
	/* Nodes information */
	/* Vx,Vy,bondm[],ro[],nu[],ep[],et[],tk[],cp[],kt[],ht[] */
	for (m1=0;m1<xnumx;m1++)
	for (m2=0;m2<ynumy;m2++)
		{
		m3=m1*ynumy+m2;
		ffscanf();
		ffscanf();
		ffscanf();pr[m3]=atof(sa);
		ffscanf();vx[m3]=atof(sa);
		ffscanf();vy[m3]=atof(sa);
		ffscanf();bondm[m3*3+0]=atoi(sa);
		ffscanf();bondm[m3*3+1]=atoi(sa);
		ffscanf();bondm[m3*3+2]=atoi(sa);
		ffscanf();exx[m3]=atof(sa);
		ffscanf();eyy[m3]=atof(sa);
		ffscanf();exy[m3]=atof(sa);
		ffscanf();sxx[m3]=atof(sa);
		ffscanf();syy[m3]=atof(sa);
		ffscanf();sxy[m3]=atof(sa);
		ffscanf();ro[m3]=atof(sa);
		ffscanf();nu[m3]=atof(sa);
		ffscanf();nd[m3]=atof(sa);
		ffscanf();po[m3]=atof(sa);
		ffscanf();mu[m3]=atof(sa);
		ffscanf();ep[m3]=atof(sa);
		ffscanf();et[m3]=atof(sa);
		ffscanf();tk[m3]=atof(sa);
		ffscanf();bondm[nodenum3+m3]=atoi(sa);
		ffscanf();cp[m3]=atof(sa);
		ffscanf();kt[m3]=atof(sa);
		ffscanf();ht[m3]=atof(sa);
		}
	/**/
	/* Gridlines positions */
	for (m1=0;m1<xnumx;m1++)
		{
		ffscanf();
		ffscanf();gx[m1]=atof(sa);
		}
	for (m2=0;m2<ynumy;m2++)
		{
		ffscanf();
		ffscanf();gy[m2]=atof(sa);
		}
	/**/
	/* Boundary Conditions Equations */
	for (m1=1;m1<bondnum;m1++)
		{
		ffscanf();
		ffscanf();bondv[m1][0]=atof(sa);
		ffscanf();bondv[m1][1]=atof(sa);
		ffscanf();bondv[m1][2]=atof(sa);
		ffscanf();bondv[m1][3]=atof(sa);
		ffscanf();bondn[m1][0]=atoi(sa);
		ffscanf();bondn[m1][1]=atoi(sa);
		ffscanf();bondn[m1][2]=atoi(sa);
		}
	/**/
	/* Markers X,Y,types */
	for (mm1=0;mm1<=marknum;mm1++)
		{
		ffscanf();markx[mm1]=atof(sa);
		ffscanf();marky[mm1]=atof(sa);
		ffscanf();markk[mm1]=atof(sa);
		ffscanf();markv[mm1]=atof(sa);
		ffscanf();markhi[mm1]=atof(sa);      /* Greg: impact history variable */
		ffscanf();markpor[mm1]=atof(sa);     /* Greg: marker porosity */
		ffscanf();markgr[mm1]=atof(sa);      /* Greg: marker grain size */
		ffscanf();marktmax[mm1]=atof(sa);    /* Tim: marker maximum temperature */
        ffscanf();markacc[mm1]=atof(sa);     /* Tim: marker accretion time */
		ffscanf();markmg_old[mm1]=atof(sa);  /* Greg: silicate magnetization variable */
		ffscanf();markmg_time[mm1]=atof(sa); /* Greg: silicate magnetization time */
		ffscanf();markt[mm1]=atoi(sa);
		}
	}
/* Load in Text Format ---------------------------- */
/**/
/**/
/**/
/* Load in Binary Format ---------------------------- */
else
	{
	fl = fopen(fl1in,"rb");
	/**/
	/* Sizes of var definition */
	szint=sizeof(n1);
	szlong=sizeof(m1);
	szfloat=sizeof(ival0);
	szdouble=sizeof(ival1);
	/* Check sizes of variables */
	fread(&szcur,1,1,fl);
	if (szcur!=szint) {printf("Current INT size <%d> is different from given in file <%d> \n",szint,szcur); exit(0);}
	fread(&szcur,1,1,fl);
	if (szcur!=szlong) {printf("Current LONG INT size <%d> is different from given in file <%d> \n",szlong,szcur); exit(0);}
	fread(&szcur,1,1,fl);
	if (szcur!=szfloat) {printf("Current FLOAT size <%d> is different from given in file <%d> \n",szfloat,szcur); exit(0);}
	fread(&szcur,1,1,fl);
	if (szcur!=szdouble) {printf("Current DOUBLE size <%d> is different from given in file <%d> \n",szdouble,szcur); exit(0);}
	/**/
	/* Grid Parameters */
	fread(&xnumx,szlong,1,fl);
	fread(&ynumy,szlong,1,fl);
	fread(&mnumx,szlong,1,fl);
	fread(&mnumy,szlong,1,fl);
	fread(&marknum,szlong,1,fl);
	fread(&xsize,szdouble,1,fl);
	fread(&ysize,szdouble,1,fl);
	fread(&gamma_eff,szdouble,1,fl);
	fread(&memory_fe,szdouble,1,fl);
	fread(&memory_si,szdouble,1,fl);
	fread(&por_init,szdouble,1,fl);
	fread(&growth_model,szlong,1,fl);
	fread(&gr_init,szdouble,1,fl);
	fread(&znumz,szlong,1,fl);
	fread(&corr2d3d,szlong,1,fl);
	fread(&pinit,szdouble,1,fl);
	fread(pkf,szdouble,4,fl);
	fread(&GXKOEF,szdouble,1,fl);
	fread(&GYKOEF,szdouble,1,fl);
	fread(&tmp_ambient,szdouble,1,fl);
	fread(&delta_tmp,szdouble,1,fl);
	fread(&timeexit,szdouble,1,fl);timeexit*=3.15576e+7;
	fread(&al2627_init,szdouble,1,fl);al2627_init*=1.0e-5;
	fread(&fe6056_init,szdouble,1,fl);fe6056_init*=1.0e-8;
	fread(&rocknum,szint,1,fl);
	fread(&bondnum,szlong,1,fl);
	fread(&n1,szint,1,fl);
	fread(&timesum,szdouble,1,fl);timesum*=3.15576e+7;
	/**/
	/* Calc,Check Grid parameters */
	gridcheck();
	/**/
	/* Rock Types information */
	fread(markim,szint,rocknum,fl);
	fread(markn0,szdouble,rocknum,fl);
	fread(markn1,szdouble,rocknum,fl);
	fread(marks0,szdouble,rocknum,fl);
	fread(marks1,szdouble,rocknum,fl);
	fread(marknu,szdouble,rocknum,fl);
	fread(markdh,szdouble,rocknum,fl);
	fread(markdv,szdouble,rocknum,fl);
	fread(markss,szdouble,rocknum,fl);
	fread(markmm,szdouble,rocknum,fl);
	fread(markll,szdouble,rocknum,fl);
	fread(marka0,szdouble,rocknum,fl);
	fread(marka1,szdouble,rocknum,fl);
	fread(markb0,szdouble,rocknum,fl);
	fread(markb1,szdouble,rocknum,fl);
	fread(marke0,szdouble,rocknum,fl);
	fread(marke1,szdouble,rocknum,fl);
	fread(markro,szdouble,rocknum,fl);
	fread(markbb,szdouble,rocknum,fl);
	fread(markaa,szdouble,rocknum,fl);
	fread(markcp,szdouble,rocknum,fl);
	fread(markkt,szdouble,rocknum,fl);
	fread(markkf,szdouble,rocknum,fl);
	fread(markkp,szdouble,rocknum,fl);
	fread(markht,szdouble,rocknum,fl);
	/**/
	/* Nodes information */
	/* Vx,Vy,bondm[],ro[],nu[],ep[],et[],tk[],cp[],kt[],ht[] */
	for (m1=0;m1<nodenum;m1++)
		{
		fread(&ival0,szfloat,1,fl);pr[m1]=(double)(ival0);
		fread(&ival0,szfloat,1,fl);vx[m1]=(double)(ival0);
		fread(&ival0,szfloat,1,fl);vy[m1]=(double)(ival0);
		fread(&m2,szlong,1,fl);bondm[m1*3+0]=m2;
		fread(&m2,szlong,1,fl);bondm[m1*3+1]=m2;
		fread(&m2,szlong,1,fl);bondm[m1*3+2]=m2;
		fread(&ival0,szfloat,1,fl);exx[m1]=(double)(ival0);
		fread(&ival0,szfloat,1,fl);eyy[m1]=(double)(ival0);
		fread(&ival0,szfloat,1,fl);exy[m1]=(double)(ival0);
		fread(&ival0,szfloat,1,fl);sxx[m1]=(double)(ival0);
		fread(&ival0,szfloat,1,fl);syy[m1]=(double)(ival0);
		fread(&ival0,szfloat,1,fl);sxy[m1]=(double)(ival0);
		fread(&ival0,szfloat,1,fl);ro[m1]=(double)(ival0);
		fread(&ival0,szfloat,1,fl);nu[m1]=(double)(ival0);
		fread(&ival0,szfloat,1,fl);nd[m1]=(double)(ival0);
		fread(&ival0,szfloat,1,fl);po[m1]=(double)(ival0);
		fread(&ival0,szfloat,1,fl);mu[m1]=(double)(ival0);
		fread(&ival0,szfloat,1,fl);ep[m1]=(double)(ival0);
		fread(&ival0,szfloat,1,fl);et[m1]=(double)(ival0);
		fread(&ival0,szfloat,1,fl);tk[m1]=(double)(ival0);
		fread(&m2,szlong,1,fl);bondm[nodenum3+m1]=m2;
		fread(&ival0,szfloat,1,fl);cp[m1]=(double)(ival0);
		fread(&ival0,szfloat,1,fl);kt[m1]=(double)(ival0);
		fread(&ival0,szfloat,1,fl);ht[m1]=(double)(ival0);
		}
	/**/
	/* Gridlines positions */
	for (m1=0;m1<xnumx;m1++)
		{
		fread(&ival0,szfloat,1,fl);gx[m1]=(double)(ival0);
		}
	for (m2=0;m2<ynumy;m2++)
		{
		fread(&ival0,szfloat,1,fl);gy[m2]=(double)(ival0);
		}
	/**/
	/* Bondary Conditions Equations */
	for (m1=1;m1<bondnum;m1++)
		{
		fread(&ival0,szfloat,1,fl);bondv[m1][0]=(double)(ival0);
		fread(&ival0,szfloat,1,fl);bondv[m1][1]=(double)(ival0);
		fread(&ival0,szfloat,1,fl);bondv[m1][2]=(double)(ival0);
		fread(&ival0,szfloat,1,fl);bondv[m1][3]=(double)(ival0);
		fread(&m2,szlong,1,fl);bondn[m1][0]=m2;
		fread(&m2,szlong,1,fl);bondn[m1][1]=m2;
		fread(&m2,szlong,1,fl);bondn[m1][2]=m2;
		}
	/**/
	/* Markers X,Y,types */
	for (mm1=0;mm1<=marknum;mm1++)
		{
		fread(&ival0,szfloat,1,fl);markx[mm1]=ival0;
		fread(&ival0,szfloat,1,fl);marky[mm1]=ival0;
		fread(&ival0,szfloat,1,fl);markk[mm1]=ival0;
		fread(&ival0,szfloat,1,fl);markv[mm1]=ival0;
		fread(&ival0,szfloat,1,fl);markhi[mm1]=ival0;   /* Greg: marker variable for impact history */
		fread(&ival0,szfloat,1,fl);markpor[mm1]=ival0;  /* Greg: marker porosity */
		fread(&ival0,szfloat,1,fl);markgr[mm1]=ival0;   /* Greg: marker grain size */
		fread(&ival0,szfloat,1,fl);marktmax[mm1]=ival0; /* Tim: marker maximum temperature */
        fread(&ival0,szfloat,1,fl);markacc[mm1]=ival0;  /* Tim: marker accretion time */
		fread(&nn3,1,1,fl);markmg_old[mm1]=ival0;       /* Greg: marker variable for silicate magnetization */
		fread(&nn2,1,1,fl);markmg_time[mm1]=ival0;      /* Greg: marker variable for silicate magnetization time */
		fread(&nn1,1,1,fl);markt[mm1]=nn1;
/*
if(mm1>850) {printf("%ld %e %e %e %d",mm1,markx[mm1],marky[mm1],markk[mm1],markt[mm1]);getchar();}
*/
		}
	}
/* Load in Binary Format ---------------------------- */
/**/
/**/
/**/
fclose(fl);
if (printmod) printf("OK!\n");
}
/* Load Information from data file ------------------------------- */



/* Print Results to data file ----------------------------------- */
void saver(int f0, int n0)
/* n0 - circle number */
/* f0 - file number */
{
/* Counters */
int n1,mm2;
char nn1,nn2,nn3;
long int m1,m2,m3;
long int mm1;
/* Buffers */
char szint,szlong,szfloat,szdouble;
float ival0;
double ival1;
/**/
/**/
/**/
if (printmod) printf("Print %d circle results to %s...",n0+1,fl1out);
/**/
/**/
/**/
/* Save data in text format ---------------------------- */
if (fl1otp==0)
	{
	fl = fopen(fl1out,"wt");
	/**/
	/* Grid Parameters */
	fprintf(fl,"%ld-xnumx\n",xnumx);
	fprintf(fl,"%ld-ynumy\n",ynumy);
	fprintf(fl,"%ld-mnumx\n",mnumx);
	fprintf(fl,"%ld-mnumy\n",mnumy);
	fprintf(fl,"%ld-marknum\n",marknum);
	fprintf(fl,"% 9.8e-xsize\n",xsize);
	fprintf(fl,"% 9.8e-ysize\n",ysize);
	fprintf(fl,"% 9.8e-gamma_eff\n",gamma_eff);
	fprintf(fl,"% 9.8e-memory_fe\n",memory_fe);
	fprintf(fl,"% 9.8e-memory_si\n",memory_si);
	fprintf(fl,"% 9.8e-por_init\n",por_init);
	fprintf(fl,"%ld-growth_model\n",growth_model);
	fprintf(fl,"% 9.8e-gr_init\n",gr_init);
	fprintf(fl,"%ld-znumz\n",znumz);
	fprintf(fl,"%ld-corr2d3d\n",corr2d3d);
	fprintf(fl,"% 9.8e-pinit\n",pinit);
	fprintf(fl,"% 9.8e-p0\n",pkf[0]);
	fprintf(fl,"% 9.8e-px\n",pkf[1]);
	fprintf(fl,"% 9.8e-py\n",pkf[2]);
	fprintf(fl,"% 9.8e-pxy\n",pkf[3]);
	fprintf(fl,"% 9.8e-GXKOEF\n",GXKOEF);
	fprintf(fl,"% 9.8e-GYKOEF\n",GYKOEF);
	fprintf(fl,"% 9.8e-tmp_ambient\n",tmp_ambient);
	fprintf(fl,"% 9.8e-delta_tmp\n",delta_tmp);
	fprintf(fl,"% 9.8e-timeexit",timeexit/3.15576e+7);
	fprintf(fl,"% 9.8e-al2627_init\n",al2627_init*1.0e+5);
	fprintf(fl,"% 9.8e-fe6056_init\n",fe6056_init*1.0e+8);
	fprintf(fl,"%d-rocknum\n",rocknum);
	fprintf(fl,"%ld-bondnum\n",bondnum);
	fprintf(fl,"%d-n0cycle\n",n0+1);
	fprintf(fl,"% 9.8e-timesum",timesum/3.15576e+7);
	fprintf(fl,"\n\n\n");
	/**/
	/* Rock Types information */
	for (n1=0;n1<rocknum;n1++)
		{
		fprintf(fl,"% 3d % 1d % 9.8e % 9.8e % 9.8e % 9.8e % 9.8e % 9.8e % 9.8e % 9.8e % 9.8e % 9.8e % 9.8e % 9.8e % 9.8e % 9.8e % 9.8e % 9.8e % 9.8e % 9.8e % 9.8e % 9.8e % 9.8e % 9.8e % 9.8e % 9.8e\n",n1,markim[n1],markn0[n1],markn1[n1],marks0[n1],marks1[n1],marknu[n1],markdh[n1],markdv[n1],markss[n1],markmm[n1],markll[n1],marka0[n1],marka1[n1],markb0[n1],markb1[n1],marke0[n1],marke1[n1],markro[n1],markbb[n1],markaa[n1],markcp[n1],markkt[n1],markkf[n1],markkp[n1],markht[n1]);
		}
	fprintf(fl,"\n\n\n");
	/**/
	/* Nodes information */
	/* Vx,Vy,bondm[],ro[],nu[],ep[],et[],tk[],cp[],kt[],ht[] */
	for (m1=0;m1<xnumx;m1++)
		{
		for (m2=0;m2<ynumy;m2++)
			{
			m3=m1*ynumy+m2;
			/* X, Y */
			fprintf(fl,"% 5.4e % 5.4e ",gx[m1],gy[m2]);
			/* P,Vx,Vy, Bond P,Vx,Vy */
			fprintf(fl,"% 9.8e % 9.8e % 9.8e %ld %ld %ld ",pr[m3],vx[m3],vy[m3],bondm[m3*3+0],bondm[m3*3+1],bondm[m3*3+2]);
			/* EPS, SIG */
			fprintf(fl,"% 9.8e % 9.8e % 9.8e % 9.8e % 9.8e % 9.8e ",exx[m3],eyy[m3],exy[m3],sxx[m3],syy[m3],sxy[m3]);
			/* ro, Nu koef */
			fprintf(fl,"% 9.8e % 9.8e % 9.8e % 9.8e % 9.8e % 9.8e % 9.8e ",ro[m3],nu[m3],nd[m3],po[m3],mu[m3],ep[m3],et[m3]);
			/* T, Bond T, Cp, Kt, Ht */
			fprintf(fl,"% 9.8e %ld % 9.8e % 9.8e % 9.8e\n",tk[m3],bondm[nodenum3+m3],cp[m3],kt[m3],ht[m3]);
			}
		fprintf(fl,"\n");
		}
	fprintf(fl,"\n\n\n");
	/**/
	/* Gridlines positions */
	for (m1=0;m1<xnumx;m1++)
		{
		fprintf(fl,"%ld % 9.8e\n",m1,gx[m1]);
		}
	fprintf(fl,"\n");
	for (m2=0;m2<ynumy;m2++)
		{
		fprintf(fl,"%ld % 9.8e\n",m2,gy[m2]);
		}
	fprintf(fl,"\n\n\n");
	/**/
	/* Boundary Conditions Equations */
	for (m1=1;m1<bondnum;m1++)
		{
		fprintf(fl,"%6ld % 9.8e % 9.8e % 9.8e % 9.8e %ld %ld %ld \n",m1,bondv[m1][0],bondv[m1][1],bondv[m1][2],bondv[m1][3],bondn[m1][0],bondn[m1][1],bondn[m1][2]);
		}
	fprintf(fl,"\n\n\n");
	/**/
	/* Markers X,Y,Type */
	for (m1=0;m1<=marknum;m1++)
		{
		mm2=markt[m1];
		fprintf(fl,"% 9.8e % 9.8e % 9.8e % 9.8e % 9.8e % 9.8e % 9.8e % 9.8e % 9.8e %d % 9.5e %d \n",markx[m1],marky[m1],markk[m1],markv[m1],markhi[m1],markpor[m1],markgr[m1],marktmax[m1],markacc[m1],markmg_old[mm1],markmg_time[mm1],mm2);
		}
	}
/* Save data in text format ---------------------------- */
/**/
/**/
/**/
/* Save data in binary format ---------------------------- */
else
	{
	fl = fopen(fl1out,"wb");
	/**/
	/* Sizes of var definition */
	szint=sizeof(n1);
	szlong=sizeof(m1);
	szfloat=sizeof(ival0);
	szdouble=sizeof(ival1);
	fwrite(&szint,1,1,fl);
	fwrite(&szlong,1,1,fl);
	fwrite(&szfloat,1,1,fl);
	fwrite(&szdouble,1,1,fl);
	/**/
	/* Grid Parameters */
	fwrite(&xnumx,szlong,1,fl);
	fwrite(&ynumy,szlong,1,fl);
	fwrite(&mnumx,szlong,1,fl);
	fwrite(&mnumy,szlong,1,fl);
	fwrite(&marknum,szlong,1,fl);
	fwrite(&xsize,szdouble,1,fl);
	fwrite(&ysize,szdouble,1,fl);
	fwrite(&gamma_eff,szdouble,1,fl);
	fwrite(&memory_fe,szdouble,1,fl);
	fwrite(&memory_si,szdouble,1,fl);
	fwrite(&por_init,szdouble,1,fl);
	fwrite(&growth_model,szlong,1,fl);
	fwrite(&gr_init,szdouble,1,fl);
	fwrite(&znumz,szlong,1,fl);
	fwrite(&corr2d3d,szlong,1,fl);
	fwrite(&pinit,szdouble,1,fl);
	fwrite(pkf,szdouble,4,fl);
	fwrite(&GXKOEF,szdouble,1,fl);
	fwrite(&GYKOEF,szdouble,1,fl);
	fwrite(&tmp_ambient,szdouble,1,fl);
	fwrite(&delta_tmp,szdouble,1,fl);
	ival1=timeexit/3.15576e+7;fwrite(&ival1,szdouble,1,fl);
	ival1=al2627_init*1.0e+5;fwrite(&ival1,szdouble,1,fl);
	ival1=fe6056_init*1.0e+8;fwrite(&ival1,szdouble,1,fl);
	fwrite(&rocknum,szint,1,fl);
	fwrite(&bondnum,szlong,1,fl);
	fwrite(&n0,szint,1,fl);
	ival1=timesum/3.15576e+7;fwrite(&ival1,szdouble,1,fl);
	/**/
	/* Rock Types information */
	fwrite(markim,szint,rocknum,fl);
	fwrite(markn0,szdouble,rocknum,fl);
	fwrite(markn1,szdouble,rocknum,fl);
	fwrite(marks0,szdouble,rocknum,fl);
	fwrite(marks1,szdouble,rocknum,fl);
	fwrite(marknu,szdouble,rocknum,fl);
	fwrite(markdh,szdouble,rocknum,fl);
	fwrite(markdv,szdouble,rocknum,fl);
	fwrite(markss,szdouble,rocknum,fl);
	fwrite(markmm,szdouble,rocknum,fl);
	fwrite(markll,szdouble,rocknum,fl);
	fwrite(marka0,szdouble,rocknum,fl);
	fwrite(marka1,szdouble,rocknum,fl);
	fwrite(markb0,szdouble,rocknum,fl);
	fwrite(markb1,szdouble,rocknum,fl);
	fwrite(marke0,szdouble,rocknum,fl);
	fwrite(marke1,szdouble,rocknum,fl);
	fwrite(markro,szdouble,rocknum,fl);
	fwrite(markbb,szdouble,rocknum,fl);
	fwrite(markaa,szdouble,rocknum,fl);
	fwrite(markcp,szdouble,rocknum,fl);
	fwrite(markkt,szdouble,rocknum,fl);
	fwrite(markkf,szdouble,rocknum,fl);
	fwrite(markkp,szdouble,rocknum,fl);
	fwrite(markht,szdouble,rocknum,fl);
	/**/
	/* Nodes information */
	/* Vx,Vy,bondm[],ro[],nu[],ep[],et[],tk[],cp[],kt[],ht[] */
	for (m1=0;m1<nodenum;m1++)
		{
		ival0=(float)(pr[m1]);fwrite(&ival0,szfloat,1,fl);
		ival0=(float)(vx[m1]);fwrite(&ival0,szfloat,1,fl);
		ival0=(float)(vy[m1]);fwrite(&ival0,szfloat,1,fl);
		m2=bondm[m1*3+0];fwrite(&m2,szlong,1,fl);
		m2=bondm[m1*3+1];fwrite(&m2,szlong,1,fl);
		m2=bondm[m1*3+2];fwrite(&m2,szlong,1,fl);
		ival0=(float)(exx[m1]);fwrite(&ival0,szfloat,1,fl);
		ival0=(float)(eyy[m1]);fwrite(&ival0,szfloat,1,fl);
		ival0=(float)(exy[m1]);fwrite(&ival0,szfloat,1,fl);
		ival0=(float)(sxx[m1]);fwrite(&ival0,szfloat,1,fl);
		ival0=(float)(syy[m1]);fwrite(&ival0,szfloat,1,fl);
		ival0=(float)(sxy[m1]);fwrite(&ival0,szfloat,1,fl);
		ival0=(float)(ro[m1]);fwrite(&ival0,szfloat,1,fl);
		ival0=(float)(nu[m1]);fwrite(&ival0,szfloat,1,fl);
		ival0=(float)(nd[m1]);fwrite(&ival0,szfloat,1,fl);
		ival0=(float)(po[m1]);fwrite(&ival0,szfloat,1,fl);
		ival0=(float)(mu[m1]);fwrite(&ival0,szfloat,1,fl);
		ival0=(float)(ep[m1]);fwrite(&ival0,szfloat,1,fl);
		ival0=(float)(et[m1]);fwrite(&ival0,szfloat,1,fl);
		ival0=(float)(tk[m1]);fwrite(&ival0,szfloat,1,fl);
		m2=bondm[nodenum3+m1];fwrite(&m2,szlong,1,fl);
		ival0=(float)(cp[m1]);fwrite(&ival0,szfloat,1,fl);
		ival0=(float)(kt[m1]);fwrite(&ival0,szfloat,1,fl);
		ival0=(float)(ht[m1]);fwrite(&ival0,szfloat,1,fl);
		}
	/**/
	/* Gridlines positions */
	for (m1=0;m1<xnumx;m1++)
		{
		ival0=(float)(gx[m1]);fwrite(&ival0,szfloat,1,fl);
		}
	for (m2=0;m2<ynumy;m2++)
		{
		ival0=(float)(gy[m2]);fwrite(&ival0,szfloat,1,fl);
		}
	/**/
	/* Boundary Conditions Equations */
	for (m1=1;m1<bondnum;m1++)
		{
		ival0=(float)(bondv[m1][0]);fwrite(&ival0,szfloat,1,fl);
		ival0=(float)(bondv[m1][1]);fwrite(&ival0,szfloat,1,fl);
		ival0=(float)(bondv[m1][2]);fwrite(&ival0,szfloat,1,fl);
		ival0=(float)(bondv[m1][3]);fwrite(&ival0,szfloat,1,fl);
		m2=bondn[m1][0];fwrite(&m2,szlong,1,fl);
		m2=bondn[m1][1];fwrite(&m2,szlong,1,fl);
		m2=bondn[m1][2];fwrite(&m2,szlong,1,fl);
		}
	/**/
	/* Markers X,Y,types */
	for (mm1=0;mm1<=marknum;mm1++)
		{
		ival0=markx[mm1];fwrite(&ival0,szfloat,1,fl);
		ival0=marky[mm1];fwrite(&ival0,szfloat,1,fl);
		ival0=markk[mm1];fwrite(&ival0,szfloat,1,fl);
		ival0=markv[mm1];fwrite(&ival0,szfloat,1,fl);
		ival0=markhi[mm1];fwrite(&ival0,szfloat,1,fl);   /* Greg: impact history variable */
		ival0=markpor[mm1];fwrite(&ival0,szfloat,1,fl);  /* Greg: marker porosity */
        ival0=markgr[mm1];fwrite(&ival0,szfloat,1,fl);   /* Greg: marker grain size */
		ival0=marktmax[mm1];fwrite(&ival0,szfloat,1,fl); /* Tim: marker maximum temperature */
        ival0=markacc[mm1];fwrite(&ival0,szfloat,1,fl);  /* Tim: marker accretion time */
		nn3=markmg_old[mm1];fwrite(&nn3,1,1,fl);         /* Greg: silicate magnetization variable */
		nn2=markmg_time[mm1];fwrite(&nn2,1,1,fl);        /* Greg: silicate magnetization time */
		nn1=markt[mm1];fwrite(&nn1,1,1,fl);
		}
	}
/* Save data in binary format ---------------------------- */
/**/
/**/
/**/
fclose(fl);
if (printmod) printf("OK!\n");
/**/
/* file.t3c file creation */
fl = fopen("file.t3c","wt");
fprintf(fl,"%d \n",f0);
fclose(fl);
/**/
/* stop.yn file information read */
fl = fopen("stop.yn","rt");
/* Read String */
ffscanf();
/**/
/* Stop Y/N */
if (sa[0]=='y' || sa[0]=='Y')
	{
	fclose(fl);
	printf("PROGRAM TERMINATED FROM stop.yn \n");
	exit(0);
	}
/**/
/* Change printmod */
if (sa[0]>='0' && sa[0]<='9')
	{
	printmod=atoi(sa);
	}
fclose(fl);
}
/* Print Results to data file ----------------------------------- */



/* LOAD WITHOUT EMPTY LINES from fl =================================== */
void ffscanf()
{
/* Counter */
int n1;
/**/
/* Read cycle */
do
	{
	/* Input string */
	n1=fscanf(fl,"%s",sa);
	/* Check end of file */
	if (n1<1)
		{
		printf("\n Unexpected end of file\n");
		fclose(fl);
		exit(0);
		}
	/* Delete last symbol <32 */
	for(n1=strlen(sa)-1;n1>=0;n1--)
	if (*(sa+n1)<=32)
	*(sa+n1)=0;
	else
	break;
	}
while (strlen(sa)<1 || sa[0]=='/');
}
/* End LOAD WITHOUT EMPTY LINES from fl =================================== */



/* LOAD WITHOUT EMPTY LINES from fl1 =================================== */
void ffscanf1()
{
/* Counter */
int n1;
/**/
/* Read cycle */
do
	{
	/* Input string */
	n1=fscanf(fl1,"%s",sa);
	/* Check end of file */
	if (n1<1)
		{
		printf("\n Unexpected end of file\n");
		fclose(fl1);
		exit(0);
		}
	/* Delete last symbol <32 */
	for(n1=strlen(sa)-1;n1>=0;n1--)
	if (*(sa+n1)<=32)
	*(sa+n1)=0;
	else
	break;
	}
while (strlen(sa)<1 || sa[0]=='/');
}
/* End LOAD WITHOUT EMPTY LINES from fl1 =================================== */


/* Calc,Check Parameters of Grid */
void gridcheck()
{
/* Gridlines NUM */
if(xnumx>MAXXLN) {printf("Space out in gx[] %ld",xnumx); exit(0);}
if(ynumy>MAXYLN) {printf("Space out in gy[] %ld",ynumy); exit(0);}
/**/
/* Nodes Num */
nodenum=xnumx*ynumy;
if(nodenum>MAXNOD) {printf("Space out in vx[],vy[] %ld",nodenum); exit(0);}
/**/
/* Cells Num */
cellnum=(xnumx-1)*(ynumy-1);
if(cellnum>MAXCEL) {printf("Space out in pr[]"); exit(0);}
/**/
/* Mark num */
if(marknum>MAXMRK+1) {printf("Space out in markx[]"); exit(0);}
/**/
/* Rock types Num */
if (rocknum>MAXTMR){printf("Space out in marknu[]"); exit(0);}
/**/
/* Bondary condit Equations Num */
if (bondnum>MAXBON){printf("Space out in bondv[]"); exit(0);}
/**/
/* Koef for processing */
xstpx=xsize/(double)(xnumx-1);
ystpy=ysize/(double)(ynumy-1);
kfx=1.0/xstpx;
kfy=1.0/ystpy;
kfxx=kfx*kfx;
kfyy=kfy*kfy;
kfxy=kfx*kfy;
/* Marker size */
mardx=xstpx/(double)(mnumx);
mardy=ystpy/(double)(mnumy);
/* Step for numerical differentiation */
numdx=5e-1*mardx;
numdy=5e-1*mardy;
/**/
/* Spec counters */
nodenum2=nodenum*2;
nodenum3=nodenum*3;
xnumx1=xnumx-1;
ynumy1=ynumy-1;
/**/
}
/* Calc,Check Parameters of Grid */
