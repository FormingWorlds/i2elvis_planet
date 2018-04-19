/* MACROS --------------------------------------------- */
/* Min, Max, Abs */
#define MINV(a,b) ((a)<=(b)? (a):(b))
#define MAXV(a,b) ((a)>=(b)? (a):(b))
#define ABSV(a)   ((a)>=0? (a):-(a))
/**/
/* Main Sizes */
#define MAXMAT  40000000 /* General Matrix Size */
#define MAXVAL  40000000 /* General Matrix Size */
#define MAXBUF  1000000 /* General Matrix Size */
#define MAXXLN  1201 /* Max vertical line num in Grid */
#define MAXYLN  1201 /* Max horizontal line num in Grid */
#define MAXXMR   20 /* Max marker in cell for X direction */
#define MAXYMR   20 /* Max marker in cell for Y direction */
#define MAXTMR   200 /* Max markers types */
#define MAXPOS   100 /* Max Pos on Line num for wi[], wn[] buffers */
#define MAXFLN  16000 /* Max Output file names Num */
#define MAXNOD    MAXXLN*MAXYLN  /* Max Total nodes num for all grid */
#define MAXCEL    (MAXXLN-1)*(MAXYLN-1)  /* Max Total nodes num for all grid */
#define MAXPAR    MAXXLN*MAXYLN*4   /* Max Total par num for all grid */
#define MAXMRK  64000000 /* Max total marker Num */
#define MAXBON    MAXPAR /* Max Total bondary condition Equations for all grid */
#define no_max   2000    /* Max number of impact events */
/* End MACROS --------------------------------------------- */



/* ARRAYS --------------------------------------------- */
/* Processing+Service Arrays for solutions */
/* val0[]        - matrix contents */
/* fre0[]        - free member for lines */
/* bufv[]        - buffer for matrix organisation */
/* lin0[]        - line numbers for matrix contents */
/* num0[]        - pos numbers for line in  val0[], pos0[] */
/* pos0[]        - first pos numbers for line in  val0[], pos0[] */
/* sol0[],sol1[] - solution buffers */
/* wn[],wi[]     - Line Num,Koef buffer */
double sol0[MAXPAR],sol1[MAXPAR],buf0[MAXBUF];
double val0[MAXMAT],fre0[MAXPAR],bufv[MAXPAR];
long int lin0[MAXMAT],num0[MAXPAR],cur0[MAXPAR];
long int pos0[MAXPAR];
long int bufn[MAXPOS];
long int wn[MAXPOS];
double wi[MAXPOS];
/**/
/**/
/**/
/* PARDISO SOLVER ARRAYS */
/* ia[] - first coloumn of each row */
/* ja[] - coloumn position of each koef in the respective row from left to right */
/* a[]  - value of each koef */
/* b[]  - Right hand side */
/* x[]  - solutions */
int ia[MAXPAR];
int ja[MAXVAL];
double b[MAXPAR];
double x[MAXPAR];
double a[MAXVAL];
/**/
/**/
/**/
/* Nodes information */
/* gx[], gy[]          - coordinates of gridlines */
/* nu[]                - Viscosity in Stokes equation, Pa*sek */
/* mu[]                - Standard viskosity for node */
/* ep[]                - Surface trace */
/* et[]                - Free array */
/* ro[]                - density, kg/m^3 */
/* tk[]                - temperature K */
/* cp[]                - Heat capacity, J/kg */
/* kt[]                - Thermal conductivity koef, Wt/m/K */
/* ht[]                - Heat Sources, W/m^3 */
/* vx[],vy[]           - Cur circle Vx,Vy,  m/sek */
/* tk0[]               - Last circle TK */
/* wa0[], wa1[]        - Water content */
/* td[]                - Thermodynamic data base */
/* bondm[]             - bondary position in bondn[],bondv[] for cur par (0) - no bondary */
/* bondn[]             - PAR1+1,PAR2+1,PAR3+1 num in bond equat: CURPAR=CONST+KOEF1*PAR1 */
/* bondv[]             - CONST,KOEF1,KOEF2,KOEF3 val in bond equat: CURPAR=CONST+KOEF1*PAR1 */
/* pr[]                - pressure, Pa */
/* pr0[]               - last cycle pressure, Pa */
/* pkf[]               - koef for pressure aproximation */
/* exx[], eyy[], exy[] - strain rates for cells (exx[], eyy[]) and nodes (exy[]) */
/* sxx[], syy[], sxy[] - stresses for cells (sxx[], syy[]) and nodes (sxy[]) */
/* gp[]                - gravity potential per unit mass */
/* ga[]                - gravity acceleration magnitude, m/s^2 */
double gx[MAXXLN],gy[MAXYLN];
double nu[MAXNOD],nd[MAXNOD],mu[MAXNOD],ep[MAXNOD],et[MAXNOD],ro[MAXNOD];
double nu0[MAXNOD],nd0[MAXNOD],ep0[MAXNOD],et0[MAXNOD],ro0[MAXNOD];
double tk[MAXNOD],cp[MAXNOD],kt[MAXNOD],ht[MAXNOD];
double cp0[MAXNOD],kt0[MAXNOD],ht0[MAXNOD];
double vx[MAXNOD],vy[MAXNOD];
double tk0[MAXNOD], tk1[MAXNOD], tk2[MAXNOD], tk3[MAXNOD];
double wa0[MAXNOD], wa1[MAXNOD];
double bondv[MAXBON][4];
double pr[MAXNOD],po[MAXNOD],po0[MAXNOD],pkf[20],gp[MAXNOD],ga[MAXNOD];
long int bondm[MAXPAR],bondn[MAXBON][3];
double exx[MAXNOD],eyy[MAXNOD],exy[MAXNOD];
double sxx[MAXNOD],syy[MAXNOD],sxy[MAXNOD];
/**/
/**/
/**/
/* Markers and Rock information */
/* Markers information */
/* markx[], marky[] - X,Y of markers */
/* markk[]          - temperature in K for marker */
/* markhi[]         - impact history variable */
/* markmg_time[]    - silicate magnetization time variable [Ma] */
/* markmg_old[]     - global silicate magnetization variable */
/* markeii          - second strain rate in variant in [1/s] for marker */
/* markt[]          - rock type of markers */
/* Information for different rock types */
/* markim[]         -  Immobility of marker type  Y(1)/No(0) */
/* marknu[], markdh[], markdv[], markss[] markmm[]                                - Koef in ductile rheology Eq */
/* markll[], marka0[], marka1[], markb0[] markb1[], marke0[], marke1[]            - Koef in brittle rheology Eq */
/* markn0[], markn1[], marks0[], marks1[]                                         - viscosity and stress limits for individual rock types */
/* markro[], markaa[], markbb[], markcp[], markkt[], markkf[], markkv[], markht[] - ro,aro,bro,cp,kt, Tkt, Pkt, ht */
/* markpor[]                                                                      - porosity [non-dim.] */
/* markgr[]									  - grain size [mum] */
/* marktmax[]                                                                     - maximum temperature [K] */
/* markacc[]                                                                      - accretion time [yr] */
float markx[MAXMRK],marky[MAXMRK];
float markk[MAXMRK],markv[MAXMRK],markp[MAXMRK],markhi[MAXMRK];
float markpor[MAXMRK],markgr[MAXMRK],marktmax[MAXMRK],markacc[MAXMRK];
char markt[MAXMRK],markmg_old[MAXMRK];
double markeii,marksii,markrii;
double markmg_time[MAXMRK];
double marknu[MAXTMR],markdh[MAXTMR],markdv[MAXTMR],markss[MAXTMR],markmm[MAXTMR];
double markll[MAXTMR],marka0[MAXTMR],marka1[MAXTMR],markb0[MAXTMR],markb1[MAXTMR],marke0[MAXTMR],marke1[MAXTMR];
double markro[MAXTMR],markaa[MAXTMR],markbb[MAXTMR],markcp[MAXTMR],markkt[MAXTMR],markht[MAXTMR],markkf[MAXTMR],markkp[MAXTMR];
double markn0[MAXTMR],markn1[MAXTMR],marks0[MAXTMR],marks1[MAXTMR];
int markim[MAXTMR];
/**/
/**/
/**/
/* Thermodynamic database parameters <loadjxxx.c> */
double td[350][350][2][5];
double tkmin,pbmin,tkstp,pbstp;
int tknum,pbnum;
/**/
/**/
/**/
/* Service, Connection Buffers */
/* vxy[]             - Cur Vx,Vy */
/* eps[]             - Cur Eps tensors */
/* errbuf[]          - Error buffer */
/* nunu[]            - NUik buffer */
/* sa[]              - input string buffer */
/* fl1in[]           - input data file name */
/* fl1itp            - input data file type */
/* fl1otp            - output data file type */
/* fl0out[]          - output data All file names */
/* fl0otp[]          - output data files types */
/* fl0stp[],fl0cyc[] - output data file X,T,time steps cyc0max for Out files */
/* fl1out[]          - output data Cur file name */
/* xn[], cn[]        - high order FD cordinates and coefficients */
/* *fl               - Streem for File Input/Output */
double vxy[20],errbuf[20],nunu[9][5],eps[500];
char sa[250];
char fl1in[50],fl1out[50];
int fl1itp,fl1otp;
char fl0out[MAXFLN][50];
int fl0otp[MAXFLN];
double fl0stp[MAXFLN][3];
int fl0cyc[MAXFLN];
double xn[1000],cn[1000][10];
FILE *fl,*fl1,*fl2,*fl3;
/**/
/* End ARRAYS --------------------------------------------- */


/* VARIABLES -------------------------------------------------- */
/**/
/* Grid Parameters */
/* xnumx,ynumy                        - num of lines in grid in X,Y directions */
/* znumz                              - num of lines in grid in Z directions */
/* xnumx1,ynumy1                      - num of cells in grid for X,Y directions */
/* mnumx,mnumy                        - num of markers in one cell for X,Y directions */
/* xsize,ysize                        - size of grid in X,Y directions, m */
/* GXKOEF,GYKOEF                      - Gravitation g for X,Y directions, m/sek^2 */
/* tmp_ambient                        - Ambient temperature of sticky air */
/* al2627_init,fe6056_init            - Initial abundance of 26Al (e-5) and 60Fe (e-8) relative to their stable isotopes */
/* pinit                              - pressure on the upper boundary, Pa */
/* gamma_eff                          - efficiency of impact heating [0-1] [Monteux et al., GRL (2007) use value of 0.3] */
/* memory_fe,memory_si                - memory of impactor's iron/silicate material + ejecta of thermal history of impactor or target body [0-1] */
/* xstpx,ystpy                        - X,Y steps for grid */
/* kfx,kfy,kfxx,kfxy,kfyy,numdx,numdy - Koef for numeric differentiation */
/* cellnum                            - total Cell Num */
/* nodenum                            - total Node Num */
/* marknum                            - total Marker Num */
/* rocknum                            - rock types num */
/* bondnum                            - bondary condition equation num */
/* nodenum3                           - nodenum*3 */
/* corr2d3d                           - 0: Accretion treated as 2D cylindrical problem, 1: Problem treated as 3D spherical */
long int xnumx,ynumy,znumz,mnumx,mnumy;
long int xnumx1,ynumy1,corr2d3d;
double xsize,ysize;
double GXKOEF,GYKOEF;
double tmp_ambient;
double al2627_init,fe6056_init;
double pinit;
double gamma_eff,memory_fe,memory_si;
double xstpx,ystpy;
double kfx,kfy,kfxx,kfxy,kfyy,numdx,numdy,mardx,mardy;
long int cellnum,nodenum;
long int marknum,bondnum;
int rocknum=0;
long int nodenum2,nodenum3;
long int leftnum,rightnum;
/**/
/**/
/**/
/* Service parameters */
/* printmod  - print information on the monitor Y(1)/N(0) */
/* intermod  - order of interpolation from nodes to markers 0-linear >=1-Non linear */
/* intermod1 - order of interpolation from markers to nodes 0-linear >=1-Non linear */
/* densimod  - mode of  density calculation: 0-constant, 1-PT-dependent */
/* outgrid   - marker move out of grid Y(0)/N(1) Orthogonal Only (2) */
/* fl0num    - number of otput file Names */
/* pos0cur   - Cur pos Counter in val0[] lin0[] */
long int pos0cur=0,printmod;
int fl0num,intermod,intermod1,outgrid=0,densimod;
/**/
/**/
/**/
/* Motion parameters */
/* cyc0max   - Num of circles for cur calculation */
/* maxxystep - max Distance change for one time step, m */
/* maxtkstep - max Temperature for one time step, K */
/* maxtmstep - max Time for one time step, sek */
/* timestep  - time step, sek */
/* timesum   - time from start, sek */
/* timeexit  - time to end simulation, sek */
/* timedir   - direction of motion in time:1 - forward, -1 - backward */
/* movemod   - solve Stokes+Continuity equations Y(1)/N(0) */
/* tempmod   - solve Heat Transfer equation Y(1)/N(0) */
/* markmod   - move markers Y(1)/N(0) */
/* ratemod   - reset velosity and pressure Y(1)/N(0) */
/* gridmod   - recalc grid parameters Y(1)/N(0) */
int cyc0max,movemod,tempmod,markmod,ratemod,gridmod;
double maxxystep,maxtkstep,maxtmstep,timestep,timesum,timeexit,timedir;
/**/
/**/
/**/
/* General iteration parameters in <gaus.c> */
/* ckoef - Cur Koef for Seidel method of iteration */
double ckoef;
/**/
/**/
/**/
/* V parameters in <move.c> */
/* DIVVMIN,STOKSMIN    - Min Valid absolut Err val for Contin,Stokes Eq */
/* stoksmod            - dNu/dT, dNu/dP Y(0)/N(1) in Stokes calc */
/* presmod             - Pressure as a function of depth Y(1)/N(0) for allinter() */
/* stoksfd             - Order of FD for EPSxx,EPSyy,EPSxy calculation */
/* nukoef              - No function */
/* nubeg,nuend,nucontr - Min, Max, Max/Min limits of viscosity */
/* hidry               - max depth with hidrostatic pressure of pore fluid */
/* hidrl               - brittle weackening factor for hidrostatic pressure of pore fluid */
/* strmin, strmax      - min max values of stress Pa allowed */
double DIVVMIN,STOKSMIN;
int stoksmod,presmod,stoksfd,viscmod;
double nubeg,nuend,nucontr,hidry,hidrl,strmin,strmax;
/**/
/**/
/**/
/* T parameters in <heat.c> */
/* HEATMIN    - Min Valid absolut Err val for Heat Equation */
/* heatmod    - dK/dX account Y(1)/N(0) */
/* heatfd     - Order of FD for Qx, Qy calculation */
/* heatdif    - Numerical diffusion coefficient */
/* frictyn    - Viscous friction heat Y(1)/N(0) */
/* adiabyn    - Adiabatic heat calculation: N(0)/Y(1) */
/* si_melting - Peridotite solidus T: 0: Hirschmann [2000], 1: Herzberg et al. [2000] */
/* fe_melting - Iron melting T: 0: Pure iron after Boehler [1993], 1: Fe-FeS eutectic after Chudinovskikh & Boehler [2007] */
double HEATMIN;
double heatdif;
int fe_melting,heatmod,heatfd,si_melting;
double frictyn,adiabyn;
/**/
/* Standard Volume in <mark2d.c> */
double rostand=0;
/* End VARIABLES -------------------------------------------------- */



/* FUNCTIONS PROTOTYPES ------------------------------------------- */
/**/
/* <load.c> FILE */
/* ffscanf()   - Load single word from file fl without empty lines */
/* ffscanf1()  - Load single word from file fl1 without empty lines */
/* loadconf()  - Load configuration from mode.t3c */
/* loader()    - Load information from data file */
/* saver()     - Save information to data file */
/* gridcheck() - Calc,Check parameters of Grid */
void ffscanf();
void ffscanf1();
int loadconf();
void loader();
void saver(int,int);
void gridcheck();
/**/
/* <gaus.c> FILE */
/* gausmat2() - Solve system of linear equation by economic frontal Gauss method */
int gausmat2(int, long int, long int);
int gausmat3(int, long int, long int);
int gausmat4(int, long int, long int);
/**/
/* <move.c> FILE */
/* viterate()               - General soubr for Vx,Vy,P Calc by Iterativ method */
/* xstokserr()              - Right part or Err in X Stokes Equat for cur node calc */
/* ystokserr()              - Right part or Err in Y Stokes Equat for cur node calc */
/* conterr()                - Right part or Err in Contin Equat for cur node calc */
/* xbonderr()               - Right part or Err in Boundary vX Equat for cur node calc */
/* ybonderr()               - Right part or Err in Boundary vY Equat for cur node calc */
/* pbonderr()               - Right part or Err in Boundary P Equat for cur cell calc */
/* maxvelstep()             - Max time step for markers definition */
/* sxxcalc(),syycalc() etc. - Value or add EPS and SIG equations */
/* fdweight()               - weight for FD calculation after Fornberg [1996] */
void viterate(int);
double xstokserr(long int, long int, int);
double ystokserr(long int, long int, int);
double conterr(long int, long int, int);
double xbonderr(long int, int);
double ybonderr(long int, int);
double pbonderr(long int, int);
void maxvelstep();
double sxxcalc(long int, long int, double);
double syycalc(long int, long int, double);
double sxycalc(long int, long int, double);
void fdweight(int, int, double);
/**/
/* <heat.c> FILE */
/* titerate()        - Temperature recalc after time step */
/* heatserr()        - Right+Right part or  Err in Heat Equat for cur node calc */
/* tbonderr()        - Right part or  Err in Boundary T Equat for cur node calc */
/* tkrecalc()        - Calc average temperature after new marker position */
/* qxcalc(),qycalc() - Coefficients or value for Qx,Qy  Equations */
void titerate(int);
double heaterr(long int, long int, int);
double tbonderr(long int, int);
void tkrecalc();
double qxcalc(long int, long int, double);
double qycalc(long int, long int, double);
/**/
/* <mark.c> FILE */
/* movemark()             - move markers by Runge-Kutta method */
/* erosmark()             - Erosion/Sedimentation Function for markers */
/* erosion()              - Erosion Surface Tracing */
/* allinter()             - Vx,Vy,P,T,EPS,SIG calc for marker by interpolation */
/* allinterv()            - Vx,Vy calc for marker by interpolation */
/* allintert()            - T calc for marker by interpolation */
/* allintere()            - EPS calc for marker by interpolation */
/* allinters()            - Vx,Vy,EPS*SIG calc for marker by interpolation */
/* allinterp()            - P calc for marker by interpolation */
/* depthp()               - P calc for marker as function of depth below the surface */
/* nodewt()               - Weights for horisontal and vertical nodes calculation for marker interpolation */ 
/* ronurecalc()           - recalc ro[],nu[] etc after new marker position */
/* dencalc()              - Calc density for given P,T and porosity after ro equation */
/* m1serch(), m2serch()   - search of nearest upper-left node of cell for  current location */
/* m1serch1(), m2serch1() - search of nearest node for current location */
/* tdbasecalc()           - compute TD properties at given P,T,Composition */
/* meltpart()             - melt fraction computing */
void movemark();
void ronurecalc();
void allinter(double, double);
void allinterv(double, double);
void allintert(double, double);
void allintere(double, double);
void allinters(double, double);
void allinterp(double, double);
void depthp(double, double);
void nodewt(long int, long int, long int, long int, double, double, int, int);
double dencalc(double, double, double, double, double, int);
double viscalc(double, double, double, double, long int, int);
long int m1serch(double);
long int m2serch(double);
void tdbasecalc(double, double, int, long int);
void meltpart1(double, double, int);
/**/
/* impact.c file */
/* impact() - Gregor's impact routine */
int impact();
int impactread();                                    /* Gregor addition */
int impactsave();                                    /* Gregor addition */
double core_radius;                       /* Gregor addition */
double M_diapir,M_diapir_ol,M_diapir_old; /* Mass of the current iron diapir pair [kg] */
double d_iron,d_iron_ol,d_iron_old;       /* Resulting new and old iron pond thicknesses [m] */
double k_cutoff       = 1.000e6;          /* Numerical cut-off value for effective thermal conductivity [W/(m*K)] */
double M_acc,M_init;                      /* Gregor addition */
double relax_time     = 0.010;            /* Time after which silicate material regains strength [Ma] (experience value) */
double start_time;                        /* Start time of N-body model after CAI formation [Ma] */
double core_form_time;                    /* Start time of thermochemical model after core formation in target body [a] */
double fe_frac;                           /* Linear iron fraction of the impactor bodies [non-dim.] */
double por_init;                          /* Initially prescribed porosity of solid silicates [non-dim.] */
double gr_init;                           /* Initial grain size for solid silicates [mum] */
int growth_model;                         /* Grain growth model: (0) No olivine grain growth, (1) Olivine grain growth only after first impact, (2) grain growth after model start */
/**/
/* core.c file */
/* core() - Gregor's core routine */
int core();
int dynamo;                               /* Is planetary dynamo likely active or not? [non-dim. switch] */
int no_merger;                            /* How many iron core-iron diapir mergers have already happened? [non-dim. counter] */
double dipl_mom;                          /* Estimated dipole moment of the planetary dynamo [A*m^2] */
double M_core;                            /* Initial mass of the preexistent iron core [kg] */
double reynolds,rossby;                   /* Reynolds and local Rossby number, can be used to assess possibility of magnetic field and magnetic reversal [non-dim.] */
double q_core,q_mantle;                   /* Current and minimum heat flux from iron core needed to drive planetary dynamo [W/m^2] */
double t_mean,t_cmb_mean;                 /* Current mean temperature of the iron core and the mean temperature of the mantle site of the CMB [K] */
/**/
/* Marker post-processing */
int markersave();                         /* Gregor addition */
/**/
/* Pebble accretion routine */
void pebbleaccr();			  /* Pebble accretion function */
double pebble_time[10000];		  /* Accretion time of pebbles from history file [Myr] */
double pebble_mass[10000];		  /* Pebble mass accreted at accretion time [mass of Ceres = 9.47e20 kg] */
/**/
/* End FUNCTIONS PROTOTYPES ------------------------------------------- */
