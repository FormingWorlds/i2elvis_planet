/* Start save marker data */
/* added by Gregor, last modified: 18/10/2011) */
int markersave()
{
long int m3;
int mmm3;
double mmm4;
/**/
/* Write out marker data, however amir.t3c should be set to different output field than composition! */
/**/
if(fl1otp!=5) /* Set amir.t3c other than composition (c) */
	{
/**/
/* Saves marker coordinates and silicate melt fraction */
	fl = fopen("marker_melt.txt","wt");
	for(m3=0;m3<=marknum;m3++)
		{
		fprintf(fl,"% 5.4e % 5.4e % 5.4e \n",markx[m3]/xsize,marky[m3]/ysize,melt_frac[m3]);
		}
	fclose(fl);
/**/
/* Saves marker coordinates and marker types */
	fl = fopen("marker_comp.txt","wt");
	for(m3=0;m3<=marknum;m3++)
		{
		mmm3=markt[m3];
		fprintf(fl,"% 5.4e % 5.4e %d \n",markx[m3]/xsize,marky[m3]/ysize,mmm3);
		}
	fclose(fl);
/**/
/* Saves marker coordinates and marker magnetization */
	fl = fopen("marker_magnet.txt","wt");
	for(m3=0;m3<=marknum;m3++)
		{
		mmm3=markmg_old[m3];
		fprintf(fl,"% 5.4e % 5.4e %d \n",markx[m3]/xsize,marky[m3]/ysize,mmm3);
		}
	fclose(fl);
/**/
/* Saves marker coordinates and marker magnetization */
	fl = fopen("marker_magnet_time.txt","wt");
	for(m3=0;m3<=marknum;m3++)
		{
		mmm4=markmg_time[m3];
		fprintf(fl,"% 5.4e % 5.4e % 5.4e \n",markx[m3]/xsize,marky[m3]/ysize,mmm4);
		}
	fclose(fl);
	}
/**/
return 0;
}
/* End marker saving */
/*====================================================================================================*/
/**/
/**/
/* Load information from configuration file amir.t3c ============== */
void loadconf()
{
/* Counter */
int n1,n2,n3;
/**/
/**/
/**/
/* Open File amir.t3c */
fl = fopen("amir.t3c","rt");
/**/
/* Load first Input File name */
ffscanf();
fl0num=0;
while(sa[0]!='~')
	{
	/* Check file Counter */
	if(fl0num>=MAXFLN) {printf("Space out in fl0out[]"); exit(0);}
	/* Save Input file name */
	for (n1=0;n1<50;n1++) fl0in[fl0num][n1]=sa[n1];
	/* Load TYPE of input - b=binary t=text */
	ffscanf(); if(sa[0] == 'b') fl0itp[fl0num]=1;
	/**/
	/* Load, Save Output File name */
	ffscanf();
	for (n1=0;n1<50;n1++) fl0out[fl0num][n1]=sa[n1];
	/* Load TYPE of output */
	/* 1 - T,t = temperature */
	/* 2 - V,v = Log viscosity */
	/* 3 -  */
	/* 4 - H,h = Log shear heating */
	/* 5 - C,c = chemical field (rock types) */
	/* 6 - D,d = density */
	/* 7 - K,k = heat conductivity */
	/* 8 - G,g = Log total heat generation */
	/* 9 - vp - vp-seismic velosity, km/s */
	/* 10 - vs - vs-seismic velosity, km/s */
	/* 11 - vp/vs - seismic velosity ratio */
	/* 12 - R,r - vorticity/shear ratio */
	/* 13 - O,o - topography */
	/* 14 - L,l - overriding plate length */
	/**/
	ffscanf(); 
	if(sa[0] == 'T' || sa[0] == 't') fl0otp[fl0num]=1;
	if(sa[0] == 'V' || sa[0] == 'v') fl0otp[fl0num]=2;
	if(sa[0] == 'H' || sa[0] == 'h') fl0otp[fl0num]=4;
	if(sa[0] == 'C' || sa[0] == 'c') fl0otp[fl0num]=5;
	if(sa[0] == 'D' || sa[0] == 'd') fl0otp[fl0num]=6;
	if(sa[0] == 'K' || sa[0] == 'k') fl0otp[fl0num]=7;
	if(sa[0] == 'G' || sa[0] == 'g') fl0otp[fl0num]=8;
	if(sa[1] == 'P' || sa[1] == 'p') fl0otp[fl0num]=9;
	if(sa[1] == 'S' || sa[1] == 's') fl0otp[fl0num]=10;
	if(sa[1] == 'R' || sa[1] == 'r') fl0otp[fl0num]=11;
	if(sa[0] == 'R' || sa[0] == 'r') fl0otp[fl0num]=12;
	if(sa[0] == 'O' || sa[0] == 'o') fl0otp[fl0num]=13;
	if(sa[0] == 'L' || sa[0] == 'l') fl0otp[fl0num]=14;
	/**/
	/* Resolution,position of the output */
	ffscanf();fl0stp[fl0num][0]=atof(sa);
	ffscanf();fl0stp[fl0num][1]=atof(sa);
	ffscanf();fl0stp[fl0num][2]=atof(sa);
	ffscanf();fl0stp[fl0num][3]=atof(sa);
	ffscanf();fl0stp[fl0num][4]=atof(sa);
	ffscanf();fl0stp[fl0num][5]=atof(sa);
	/**/
	/* Incr File Counters */
	fl0num++;
	/**/
	/* Load Next Input File names */
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
fclose(fl);
/* End Load information from configuration file amir.t3c */
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
}
/* Load information from configuration file mode.t3c ============== */



/* Load Information from data file ------------------------------- */
void loader()
/* bondv[] - boundary value */
/* bondm[] - boundary mode 0=Not, -1=Value, 1,2...=LinNum+1 */
/* m1,m2 - node X,Y number */
{
/* Counter */
int n1;
char nn1,nn2,nn3;
long int m1,m2,m3,nv[100];
long int mm1,mm2,erosmark=0;
char szint,szlong,szfloat,szdouble,szcur;
float ival0,iv[100];
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
		ffscanf();markhi[mm1]=atof(sa);      /* Greg: marker impact history variable */
		ffscanf();markpor[mm1]=atof(sa);     /* Greg: marker porosity */
		ffscanf();markgr[mm1]=atof(sa);      /* Greg: marker grain size */
		ffscanf();marktmax[mm1]=atof(sa);    /* Tim: marker maximum temperature */
		ffscanf();markacc[mm1]=atof(sa);     /* Tim: marker accretion time */
		ffscanf();markmg_old[mm1]=atof(sa);  /* Greg: marker magnetization variable */
		ffscanf();markmg_time[mm1]=atof(sa); /* Greg: marker magnetization time */
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
		fread(iv,szfloat,3,fl);
		pr[m1]=(double)(iv[0]);
		vx[m1]=(double)(iv[1]);
		vy[m1]=(double)(iv[2]);
		fread(nv,szlong,3,fl);
		bondm[m1*3+0]=nv[0];
		bondm[m1*3+1]=nv[1];
		bondm[m1*3+2]=nv[2];
		fread(iv,szfloat,14,fl);
		exx[m1]=(double)(iv[0]);
		eyy[m1]=(double)(iv[1]);
		exy[m1]=(double)(iv[2]);
		sxx[m1]=(double)(iv[3]);
		syy[m1]=(double)(iv[4]);
		sxy[m1]=(double)(iv[5]);
		ro[m1]=(double)(iv[6]);
		nu[m1]=(double)(iv[7]);
		nd[m1]=(double)(iv[8]);
		po[m1]=(double)(iv[9]);
		mu[m1]=(double)(iv[10]);
		ep[m1]=(double)(iv[11]);
		et[m1]=(double)(iv[12]);
		tk[m1]=(double)(iv[13]);
		fread(&m2,szlong,1,fl);bondm[nodenum3+m1]=m2;
		fread(iv,szfloat,3,fl);
		cp[m1]=(double)(iv[0]);
		kt[m1]=(double)(iv[1]);
		ht[m1]=(double)(iv[2]);
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
		fread(iv,szfloat,4,fl);
		bondv[m1][0]=(double)(iv[0]);
		bondv[m1][1]=(double)(iv[1]);
		bondv[m1][2]=(double)(iv[2]);
		bondv[m1][3]=(double)(iv[3]);
		fread(nv,szlong,3,fl);
		bondn[m1][0]=nv[0];
		bondn[m1][1]=nv[1];
		bondn[m1][2]=nv[2];
		}
	/**/
	/* Markers X,Y,types */
	for (mm1=0;mm1<=marknum;mm1++)
		{
		fread(iv,szfloat,7,fl);
		markx[mm1]=iv[0];
		marky[mm1]=iv[1];
		markk[mm1]=iv[2];
		markv[mm1]=iv[3];
		markhi[mm1]=iv[4];                       /* Greg: marker history variable */
		markpor[mm1]=iv[5];                      /* Greg: marker porosity */
        markgr[mm1]=iv[6];                       /* Greg: marker grain size */
		marktmax[mm1]=iv[7];                     /* Tim: marker maximum temperature  */
        markacc[mm1]=iv[8];                      /* Tim: marker accretion time */
		fread(&nn3,1,1,fl);markmg_old[mm1]=nn3;  /* Greg: marker magnetization variable */
		fread(&nn2,1,1,fl);markmg_time[mm1]=nn2; /* Greg: marker magnetization time */
		fread(&nn1,1,1,fl);markt[mm1]=nn1;
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
void saver(int n0)
/* n0 - file number */
{
/* Counters */
int n1,n2,mm2,mm3,ccur,cmin=0,cmax=0,yn;
char nn1,nn2,nn3;
long int m1,m2,m3,m4,m5,m6,m7;
long int mm1,erosmark=0;
/* Buffers for XY */
double x,y;
char szint,szlong,szfloat,szdouble;
float ival0;
double ival, ival1,xresx,yresy,xminx,xmaxx,yminy,ymaxy,ea,na,mpb,mtk,mro,mbb,mnu,mkt,mcp,mht,mvp,mvs,xcur,ycur,mwa;
long int xresol,yresol,nresol,markstp;
long int xresol0,xresolsum;
double xincr,xminx0,yminy0;
/* TD Database variables  */
double n,e;
/**/
/**/
/**/
if (printmod) printf("Print data to %s...",fl1out);
/**/
/**/
/**/
/* Open file */
fl = fopen(fl1out,"wt");
/* Read check resolution */
xincr=xresol=xresol0=(long int)(fl0stp[n0][2]);
yresol=(long int)(fl0stp[n0][5]);
nresol=xresol0*yresol;
/* Check coordinates, calc step */
if (fl0stp[n0][0]<0) fl0stp[n0][0]=0;
if (fl0stp[n0][1]>1.0) fl0stp[n0][1]=1.0;
if (fl0stp[n0][3]<0) fl0stp[n0][3]=0;
if (fl0stp[n0][4]>1.0) fl0stp[n0][4]=1.0;
xminx0=xminx=fl0stp[n0][0]*xsize;
xmaxx=fl0stp[n0][1]*xsize;
xresx=(xmaxx-xminx)/(double)(xresol-1);
yminy0=yminy=fl0stp[n0][3]*ysize;
ymaxy=fl0stp[n0][4]*ysize;
yresy=(ymaxy-yminy)/(double)(yresol-1);
yminy-=yresy; if(yminy<0) yminy=0;
ymaxy+=yresy; if(ymaxy>ysize) ymaxy=ysize;
if((nresol*2)>MAXMAT)
	{
	printf("\n Limited space in val0[] lin0[], spleating ...\n");
	xincr=(long int)((double)(MAXMAT)/(double)(yresol*2));
	}
/* X, Y Resolution save */
/*
printf("%ld %ld %ld   %e %e %e %e %e %e\n",xresol,yresol,nresol,xminx,xmaxx,xresx,yminy,ymaxy,yresy); getchar();
*/
fprintf(fl,"%e \n %ld %ld \n",timesum/3.15576e+7,xresol0,yresol);
/* Step for markers definition */
ival=0.1*(double)(marknum)/(double)(nresol);
markstp=(long int)(ival);
if (markstp<1) markstp=1;
/**/
/**/
/**/
/* Compute vorticity tensor */
for (m1=1;m1<xnumx-1;m1++)
for (m2=1;m2<ynumy-1;m2++)
	{
	/* Pos of in ol0[] */
	m3=m1*ynumy+m2;
	/**/
	/* Min,Max Vx definition */
	esp[m3]=spncalc(m1,m2);
/*
printf("%d %ld %ld %ld %e %e %e %e ",fl1otp,m1,m2,m3,exy[m3],eps[0],eps[1],esp[m3]); getchar();
*/
	}
/* Compute vorticity tensor */
/**/
/**/
/**/
/* Save data in text format ---------------------------- */
xresolsum=0;
do
	{
	/* Set x resolution */
	xresol=xincr;
	if(xresolsum+xresol>xresol0) xresol=xresol0-xresolsum;
	nresol=xresol*yresol;
	xminx=xminx0=fl0stp[n0][0]*xsize+xresx*((double)(xresolsum));
	xminx-=xresx; if(xminx<0) xminx=0;
	xmaxx=fl0stp[n0][0]*xsize+xresx*((double)(xresolsum+xresol));
	if(xmaxx>xsize) xmaxx=xsize;
	/**/
/*
printf("%ld %ld %ld %ld  %e %e %e %e %e %e\n",xresol,yresol,nresol,markstp,xminx,xmaxx,xresx,yminy,ymaxy,yresy); getchar();
*/
	/* Clear visual arrays */
	for (m1=0;m1<nresol;m1++)
		{
		val0[m1]=val0[nresol+m1]=0;
		lin0[m1]=lin0[nresol+m1]=0;
		}
	/**/
	/**/
	/**/
	/* CHEMICAL COMPONENT VISUALISATION ----------------------------*/
	/* Variations in markers types for chemical component visualisation */
	if(fl1otp==5)
		{
		cmin=1000,cmax=0;
		/* Markers cycle */
		for (m3=0;m3<=marknum;m3+=markstp) 
		/* Check markers out of visualisation area */
		if ((double)(markx[m3])>xminx && (double)(markx[m3])<xmaxx && (double)(marky[m3])>yminy && (double)(marky[m3])<ymaxy && markt[m3]<50)
			{
			if (cmin>markt[m3]) cmin=markt[m3];
			if (cmax<markt[m3]) cmax=markt[m3];
			}
		/**/
		/* Marker type cycle for chemical component Visualisation */
		for (ccur=cmin;ccur<=cmax;ccur++) 
			{
			/* Markers cycle */
			for (m3=0;m3<=marknum;m3+=markstp)
			/* Check markers out of visualisation area */
			if ((double)(markx[m3])>xminx && (double)(markx[m3])<xmaxx && (double)(marky[m3])>yminy && (double)(marky[m3])<ymaxy && markt[m3]==ccur)
				{
/*
printf("%d %d %d   %ld %d   %e %e %e \n",cmin,cmax,ccur,m3,markt[m3],markx[m3],marky[m3],markk[m3]); getchar();
*/
				/* Define relative position of the marker */
				ea=((double)(markx[m3])-xminx0)/xresx+1.0;
				m1=(long int)(ea);
				if(m1<0) m1=0;
				if(m1>xresol) m1=xresol;
				ea-=(double)(m1);
				na=((double)(marky[m3])-yminy0)/yresy+1.0;
				m2=(long int)(na);
				if(m2<0) m2=0;
				if(m2>yresol) m2=yresol;
				na-=(double)(m2);
				/* Add weights for four nodes surrounding the marker */
				m4=nresol+(m1-1)*yresol+(m2-1);
				if(m1>0 && m2>0) val0[m4]+=(float)((1.0-ea)*(1.0-na));
				if(m1>0 && m2<yresol) val0[m4+1]+=(float)((1.0-ea)*na);
				if(m1<xresol && m2>0) val0[m4+yresol]+=(float)(ea*(1.0-na));
				if(m1<xresol && m2<yresol) val0[m4+yresol+1]+=(float)(ea*na);
/*
if(xresolsum){printf("%ld %ld %ld %ld  %e %e %e %e \n",m3,m1,m2,m4,markx[m3],marky[m3],ea,na); getchar();}
*/
				}
			/**/
			/* Change types for nodes of visual arrays */
			for (m1=0;m1<nresol;m1++)
				{
				if(val0[nresol+m1]>val0[m1]) 
					{
					val0[m1]=val0[nresol+m1]; 
					lin0[m1]=ccur;
					/* Clear value */
					val0[nresol+m1]=0; 
					}
				}
			}
		/**/
		/* Marker type reload for chemical component Visualisation */
		yn=0;
		for (m1=0;m1<nresol;m1++)
			{
			if(!val0[m1]) 
				{
				lin0[m1]=-1;
				yn=1;
				}
			}
		/**/
		/* Marker type reinterpolate for empty nodes */
		if(gridmod && yn)
			{
			/* Node cycle */
			for (m1=0;m1<xresol;m1++)
			for (m2=0;m2<yresol;m2++)
				{
				m3=m1*yresol+m2;
				/* Interpolate marker types from other nodes */
				if(lin0[m3]==-1) 
					{
					/* Serch for surrounding non empty nodes */
					m4=1;
					yn=0;
					do
						{
						cmin=1000,cmax=0;
						for (m5=m1-m4;m5<=m1+m4;m5++)
						for (m6=m2-m4;m6<=m2+m4;m6++)
						if (m5>=0 && m5<xresol && m6>=0 && m6<yresol)
							{
							m7=m5*yresol+m6;
							if(lin0[m7]!=-1)
								{
								yn=1;
								if (cmin>lin0[m7]) cmin=lin0[m7];
								if (cmax<lin0[m7]) cmax=lin0[m7];
								}
							}
						m4++;
						}
					while(!yn && m4<=gridmod);
					/**/
					/* Recalc using non empty nodes */
					if(yn) 
						{
						/* Marker type cycle for chemical component Visualisation */
						for (ccur=cmin;ccur<=cmax;ccur++) 
							{
							for (m5=m1-m4+1;m5<m1+m4;m5++)
							for (m6=m2-m4+1;m6<m2+m4;m6++)
							if (m5>=0 && m5<xresol && m6>=0 && m6<yresol)
								{
								m7=m5*yresol+m6;
								/* Add weight */
								if(lin0[m7]==ccur)
									{
									ea=ABSV(((double)(m5-m1))/(double)(m4));
									na=ABSV(((double)(m6-m2))/(double)(m4));
									val0[nresol+m3]+=(float)((1.0-ea)*(1.0-na))*val0[m7];
/*
									val0[nresol+m3]+=(float)((1.0-ea)*(1.0-na));
{printf("%ld %ld %ld   %ld   %ld %ld %ld   %d \n",m1,m2,m3,m4,m5,m6,m7,ccur); getchar();}
*/
									}
								}
							/* Set new marker type */
							if(val0[nresol+m3]>val0[m3]) 
								{
								val0[m3]=val0[nresol+m3];
								lin0[nresol+m3]=ccur;
								/* Clear value */
								val0[nresol+m3]=0; 
								}
							}
						}
					}
				}
			/**/
			/* Marker type reload for chemical component Visualisation */
			for (m1=0;m1<nresol;m1++)
				{
				if(lin0[m1]==-1 && val0[m1]) 
					{
					lin0[m1]=lin0[nresol+m1];
					lin0[nresol+m1]=0;
					yn=1;
					}
				}
			}
		}
	/* End CHEMICAL COMPONENT VISUALISATION ----------------------------*/
	/**/
	/**/
	/**/
	/* TRANSPORT PROPERTIES VISUALISATION +++++++++++++++++++++++++++++ */
	if(fl1otp!=5)
		{
/**/
/**/
		/* Markers cycle (not for shear heating) */
		if(fl1otp!=4 && fl1otp!=12) 
		for (m3=0;m3<=marknum;m3+=markstp) 
		/* Check markers out of visualisation area */
		if ((double)(markx[m3])>xminx && (double)(markx[m3])<xmaxx && (double)(marky[m3])>yminy && (double)(marky[m3])<ymaxy && markt[m3]<50)
			{
/*
printf("%ld %d   %e %e %e \n",m3,markt[m3],markx[m3],marky[m3],markk[m3]); getchar();
*/
			/* P, T parameters calc */
			allinter1((double)(markx[m3]),(double)(marky[m3]));
			mpb=eps[10]*1e-5;
			/* Reset marker temperature for newly coming markers */
			if(markk[m3]<1.0) markk[m3]=(float)(eps[2]);
			mtk=(double)(markk[m3]);
			/* Set starting value for marker maximum temperature field (Tim, 28/12/2016) */
			if(markt[m3]!=0)
                		{
                        	if(marktmax[m3]<markk[m3]) marktmax[m3]=(float)(markk[m3]);
				}
			/**/
			/* Marker type */
			mm2=markt[m3];
			/* Melt fraction */
			meltpart1(mtk,mpb,mm2);
			melt_frac[m3]=eps[21];
			/**/
			/* Marker Properties */
			mnu=markv[m3];
			/* Min,Max NU limitation */
			if(mnu<nubeg) mnu=nubeg;
			if(mnu>nuend) mnu=nuend;
			mro=dencalc(mtk,mpb,(double)(markpor[m3]),(double)(markx[m3]),(double)(marky[m3]),mm2);
			mbb=eps[20];
			mcp=markcp[mm2];
			mkt=(markkt[mm2]+markkf[mm2]/(mtk+77.0))*exp(markkp[mm2]*mpb);
			/* Test Heat conductivity k=ko/(1+b*(T-To)/To) */
			if (markkt[mm2]<0) mkt=-markkt[mm2]/(1.0+markkf[mm2]*(mtk-markkp[mm2])/markkp[mm2]);
			mht=markht[mm2];
			/**/
			/**/
			/**/
			/* Thermodynamic database use for ro, Cp */
			if (densimod==2)
			if(mm2==5 || mm2==6 || mm2==11 || mm2==12)
				{
				/* Compute TD variables */
				eps[47]=mkt;
				tdbasecalc(mtk,mpb,mm2,mm1);
				mro=eps[41];
				mcp=eps[43];
				mbb=eps[44];
				mkt=eps[47];
				/**/
				}
			/**/
/*
printf("num=%ld type=%d  x=%e y=%e mpb=%e mtk=%e nu=%e ro=%e cp=%e kt=%e ht=%e",mm1,mm2,markx[mm1],marky[mm1],mpb,mtk,mnu,mro,mcp,mkt,mht);getchar();
*/
			/**/
			/* Define parameter for the visualisation */
			/* 1 - T,t = temperature */
			/* 2 - V,v = Log viscosity */
			/* 3 - S,s = bulk strain */
			/* 4 - H,h = Log shear heating */
			/* 5 - C,c = chemical field (rock types) */
			/* 6 - D,d = density */
			/* 7 - K,k = heat conductivity */
			/* 8 - G,g = Log total heat generation */
			/* 9 - vp - vp-seismic velosity, km/s */
			/* 10 - vs - vs-seismic velosity, km/s */
			/* 11 - vr - vp/vs */
			/* 12 - R,r vorticity/strainrate ratio */
			switch(fl1otp)
				{
				/* 1 - T,t = temperature */
				case 1: ival=mtk; break;
				/**/
				/* 2 - V,v = Log viscosity */
				case 2: ival=log(mnu)/log(10.0); break;
				/**/
				/* 3 -  */
				case 3: break;
				/**/
				/* 4 - H,h = Log shear heating */
				/* EPSxx*SIGxx  EPSyy*SIGyy  EPSxy*SIGxy interpolated values  */
				case 4: ival=2.0*eps[13]+eps[14]+eps[15]; 
				if(ival<1e-30) ival=1e-30; 
				ival=log(ival)/log(10.0); break;
				break;
				/**/
				/* 6 - D,d = density */
				case 6: ival=mro; break;
				/**/
				/* 7 - K,k = heat conductivity */
				case 7: ival=mkt; break;
				/**/
				/* 8 - G,g = Log total heat generation */
				case 8: ival=mht+2.0*eps[13]+eps[14]+eps[15]; break;
				if(ival<1e-30) ival=1e-30; 
				ival=log(ival)/log(10.0); 
				break;
				/**/
				/* 9 - vp - vp-seismic velosity, km/s */
				case 9: ival=mvp; break;
				/**/
				/* 10 - vs - vs-seismic velosity, km/s */
				case 10: ival=mvs; break;
				/**/
				/* 11 - vr -  vp/vs */
				case 11:
				ival=0;
				if(mvs>0) ival=mvp/mvs; 
				break;
				/**/
				/* 12 - R,r - vorticity/shear ratio */
				case 12: ival=pow(eps[4]*eps[4]+0.5*(eps[6]*eps[6]+eps[8]*eps[8]),0.5); 
/*
printf("%ld %e %e %e %e %e %e %e \n",m3,markx[m3],marky[m3],eps[4],eps[6],eps[8],eps[27],ival); getchar();
*/
				if(ival<1e-30) ival=1e-30;
				ival=eps[27]/ival;
				break;
				}
			/**/
			/* Define relative position of the marker */
			ea=((double)(markx[m3])-xminx0)/xresx+1.0;
			m1=(long int)(ea);
			if(m1<0) m1=0;
			if(m1>xresol) m1=xresol;
			ea-=(double)(m1);
			na=((double)(marky[m3])-yminy0)/yresy+1.0;
			m2=(long int)(na);
			if(m2<0) m2=0;
			if(m2>yresol) m2=yresol;
			na-=(double)(m2);
/*
if(m2>yresol-2){printf("%ld %ld %e %e %e %e %e ",m3,m1,marky[m3],na,yminy,ymaxy,yresy); getchar();}
if(m1>xresol-2){printf("%ld %ld %e %e %e %e %e ",m3,m1,markx[m3],ea,xminx,xmaxx,xresx); getchar();}
*/
			/* Add values and weights for four nodes surrounding the marker */
			m4=(m1-1)*yresol+(m2-1);
			if(m1>0 && m2>0) {val0[m4]+=(float)(ival*(1.0-ea)*(1.0-na)); val0[nresol+m4]+=(float)((1.0-ea)*(1.0-na));}
			if(m1>0 && m2<yresol) {val0[m4+1]+=(float)(ival*(1.0-ea)*na); val0[nresol+m4+1]+=(float)((1.0-ea)*na);}
			if(m1<xresol && m2>0) {val0[m4+yresol]+=(float)(ival*ea*(1.0-na)); val0[nresol+m4+yresol]+=(float)(ea*(1.0-na));}
			if(m1<xresol && m2<yresol) {val0[m4+yresol+1]+=(float)(ival*ea*na); val0[nresol+m4+yresol+1]+=(float)(ea*na);}
/*
if(xresolsum){printf("%ld %ld %ld %ld  %e %e %e %e %e \n",m3,m1,m2,m4,markx[m3],marky[m3],ea,na); getchar();}
*/
			}
		/* End Markers cycle (not for shear heating) */
		/**/
		/**/
		/**/
		/* Transport properties recalc */
		for (m1=0;m1<xresol;m1++)
		for (m2=0;m2<yresol;m2++)
			{
			m4=m1*yresol+m2;
			/* Interpolate parameters from markers */
			if(val0[nresol+m4]) 
				{
				val0[m4]/=val0[nresol+m4];
				}
			/* Interpolate other parameters from grid */
			else
				{
				xcur=xminx0+xresx*((double)(m1));
				ycur=yminy0+yresy*((double)(m2));
				/* Interpolated parameters calc */
				allinter1(xcur,ycur);
				/* Define parameter for the visualisation */
				/* 1 - T,t = temperature */
				/* 2 - V,v = Log viscosity */
				/* 3 - S,s = bulk strain */
				/* 4 - H,h = Log shear heating */
				/* 5 - C,c = chemical field (rock types) */
				/* 6 - D,d = density */
				/* 7 - K,k = heat conductivity */
				/* 8 - G,g = Log total heat generation */
				/* 9 - vp - vp-seismic velosity, km/s */
				/* 10 - vs - vs-seismic velosity, km/s */
				/* 11 - vr - vp/vs */
				/* 12 - R,r vorticity/strainrate ratio */
				switch(fl1otp)
					{
					/* 1 - T,t = temperature */
					case 1: val0[m4]=(float)(eps[2]); break;
					/**/
					/* 2 - V,v = Log viscosity */
					case 2: val0[m4]=(float)(log(eps[16])/log(10.0)); break;
					/**/
					/* 3 - S,s = bulk strain */
					/* Calc Second strain Tenzor invariant EEii=(1/2SUM(EPSik^2))^(1/2) */
					case 3: val0[m4]=0; break;
					/**/
					/* 4 - H,h = Log shear heating */
					/* EPSxx*SIGxx  EPSyy*SIGyy  EPSxy*SIGxy interpolated values  */
					case 4: val0[m4]=(float)(2.0*eps[13]+eps[14]+eps[15]); 
					if(val0[m4]<1e-30) val0[m4]=1e-30; 
					val0[m4]=(float)(log(val0[m4])/log(10.0)); break;
					break;
					/**/
					/* 6 - D,d = density */
					case 6: val0[m4]=(float)(eps[17]); break;
					/**/
					/* 7 - K,k = heat conductivity */
					case 7: val0[m4]=(float)(eps[18]); break;
					/**/
					/* 8 - G,g = Log total heat generation */
					case 8: val0[m4]=(float)(eps[19]+2.0*eps[13]+eps[14]+eps[15]); break;
					if(val0[m4]<1e-30) val0[m4]=1e-30; 
					val0[m4]=(float)(log(val0[m4])/log(10.0)); 
					break;
					/* 9 - vp - vp-seismic velosity, km/s */
					case 9: val0[m4]=0; break;
					/* 10 - vs - vs-seismic velosity, km/s */
					case 10: val0[m4]=0; break;
					/**/
					/* 11 - vr -  vp/vs */
					case 11: val0[m4]=0; break;
					/**/
					/* 12 - R,r - vorticity/shear ratio */
					case 12: ival=pow(eps[4]*eps[4]+0.5*(eps[6]*eps[6]+eps[8]*eps[8]),0.5); 
/*
printf("%ld %ld %e %e %e %e %e %e %e \n",m1,m2,xcur,ycur,eps[4],eps[6],eps[8],eps[27],ival); getchar();
*/
					if(ival<1e-30) ival=1e-30; 
					val0[m4]=eps[27]/ival; 
					break;
					}
				}
			}
		}
	/* End TRANSPORT PROPERTIES VISUALISATION +++++++++++++++++++++++++++++ */
	/**/
	/**/
	/**/
	/* Visualisation parameters output */
	for (m1=0;m1<xresol;m1++)
	for (m2=0;m2<yresol;m2++)
		{
		/* Visual node number */
		m4=m1*yresol+m2;
		/* Visual node coordinates */
		xcur=xminx0+xresx*((double)(m1));
		ycur=yminy0+yresy*((double)(m2));
		/* Value for Visualisation */
		if(fl1otp!=5) 
			{
			fprintf(fl," % 5.4e",val0[m4]);
			}
		else
			{
			/*Compression of rock type */
			m5=1;
			for (m3=m2+1;m3<yresol;m3++)
				{
				if(lin0[m4]==lin0[m4+m3-m2]) m5++; else break;
				}
			if(m5>3)
				{
				fprintf(fl," %d %ld %ld",-2,m5,lin0[m4]);
				m2=m3-1;
				}
			else
				{
				fprintf(fl," %ld",lin0[m4]);
				}
			}
/*
		fprintf(fl,"% 9.8e % 9.8e % 9.8e \n",xcur,ycur,val0[m4]);
	fprintf(fl,"\n\n\n");
*/
		}
	/**/
	/* Add x resolution */
	xresolsum+=xresol;
	}
while(xresolsum<xresol0);
/* Save data in text format ---------------------------- */
/**/
/**/
/**/
fclose(fl);
if (printmod) printf("OK!\n");
/**/
/**/
/* Save marker informations for post-processing */
if (printmod) printf("Read marker data\n");
markersave();
if (printmod) printf("Marker data saved!\n");
/**/
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



/* Calculation of Vx,Vy,EPSxx,EPSyy,EPSxy,SIGxx,SIGyy,SIGxy,T,T0,T1,T2,P for current location by Linear Interpolation */
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
void allinter1(double x, double y)
/* x,y - XY location of point for Vx,Vy calc */
{
/* Counters */
long int m1,m2,m3,m10,m20,m1min,m1max,m2min,m2max;
/* en-NormalisedDistance */
double e,n,ival;
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
/* Buffer clear */
for (m1=0;m1<=40;m1++) eps[m1]=0;
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
/* Interpolation weights calc  after Fornberg (1996) */
nodewt(m1min,m1max,m2min,m2max,x,y,0,0);
/**/
/* T Interpolate after interpolation weights */
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
	/* Material properties */
	eps[16]+=ival*nu[m3];
	eps[17]+=ival*ro[m3];
	eps[18]+=ival*kt[m3];
	eps[19]+=ival*ht[m3];
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
/* SIGxy,EPSxy, SIGxy*EPSxy, Esp interpolation ------------------------ */
/* Horizontal,Vertical limits for interpolation calc */
m1min=m10; if(m1min<1) m1min=1; if(m1min>xnumx-3) m1min=xnumx-3;
m1max=m1min+1+intermod; if(m1max>xnumx-2) m1max=xnumx-2;
m1min=m1min-intermod; if(m1min<1) m1min=1;
/**/
m2min=m20; if(m2min<1) m2min=1; if(m2min>ynumy-3) m2min=ynumy-3;
m2max=m2min+1+intermod; if(m2max>ynumy-2) m2max=ynumy-2;
m2min=m2min-intermod; if(m2min<1) m2min=1;
/**/
/* Interpolation weights calc  after Fornberg (1996) */
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
	eps[27]+=ival*esp[m3];
	}
/* End SIGxy,EPSxy interpolation ------------------------ */
/**/
/**/
/**/
/* SIGxx,EPSxx,SIGyy,EPSyy,P, SIGxx*EPSxx, SIGyy*EPSyy interpolation ------------------------ */
/* Horizontal,Vertical limits for interpolation calc */
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
/* Interpolation weights calc  after Fornberg (1996) */
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
/* Interpolation weights calc  after Fornberg (1996) */
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
/* Interpolation weights calc  after Fornberg (1996) */
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
/*
printf("eps %e %e ",m1,m2,e,n); getchar();
*/
}
/* Calculation of Vx,Vy,EPSxx,EPSyy,EPSxy,SIGxx,SIGyy,SIGxy,T,T0,T1,T2,P for current location by Linear Interpolation */



/* Spin tensor calc */ 
/* Spn=1/2(dVy/dX-dVx/dY) */
double spncalc(long int m1, long int m2)
/* m1,m2 - node X,Y number */
{
/* Exy position */
double xi,leftsxy=0,leftsxy1=0,leftsxy2=0;
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
/* Return Spn val ----------------------------*/
/* Exy=1/2(dVx/dY+dVy/dX)=0 */
/* 1/2dVx/dY */
/* Add Vx with koefficients */
for (m3=m2min;m3<=m2max;m3++)
	{
	leftsxy1+=vx[m1*ynumy+m3]*cn[m3-m2min][2]/2.0;
	}
/**/
/* 1/2dVy/dX */
/* Add Vy with koefficients */
for (m3=m1min;m3<=m1max;m3++)
	{
	leftsxy2+=vy[m3*ynumy+m2]*cn[m3-m1min][1]/2.0;
	}
/**/
/* Save Exy */
eps[0]=leftsxy=leftsxy1+leftsxy2;
/* Save Esp (rotation rate) */
eps[1]=leftsxy2-leftsxy1;
return leftsxy2-leftsxy1;
}
/* Spin tensor calc */ 


