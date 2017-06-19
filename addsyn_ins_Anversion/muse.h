#ifndef MUSE_H
#define MUSE_H

/******************************************************************************
 *
 *    header file containing definitions needed by M4C unit and function
 *    generators and by user's instrument and orchestra files.
 *   MUSE.H   -  Collection of global variables, structure
 *    For        templates, and macros available for user
 *    M4C        Orchestra development.
 *
 *   file organized by R. Maher  16-Jul-86 
 *   11/21/87  talmage	changed _afcard and added defn of maptypF
 *    9/15/88  maher    modified to reflect new unit gens, structures
 *    7/23/90  jwb	added NLPTABLEN, NLPTAB2 definitions
 *    1/29/91	jwb	For NeXT install: remove char *calloc()
 *    2/11/91	jwb	For Next:  remove abs() macro define
 *    6/  /91  kelley   new ADeSR structure
 *    7/20/91	jwb	Add complete set of nonvoid function defs
 *    7/21/91   jwb     install NHMAX, NHM2 defines
 *    9/8/91    jwb	Install structure and function def for comb and
 *			allpass ugens
 *    9/25/91	jwb	Install structure and function def for reverb ugen
 *    10/91     jwb     muse.h replaces common_m4c.h
 *    11/91   camilleg  Add defs for InitDone() and NoQuad()
 *    1/19/92 camilleg  Unify different versions of muse.h 
 *    7/3/92  camilleg	vfiltset() renamed to varbpfset()
 *    3/28/92  jwb      add Dfunc function table pointers (double **)
 *                      add definition of TERM, used to terminate arrays
 *    4/12/92  jwb      add def of NABMAX for Alvbr table for nlf synth
 *                      and declaration of alvbr() function
 *    4/06/92 camilleg  strings as scorecard params
 *    4/20/92 camilleg  random() renamed to rand01()
 *    5/28/92 camilleg  "endnote" added to maptyp macro;
 *			TacetNote(), OrchNil() declarations added;
 *			"#ifndef MUSE_H / #define MUSE_H / #endif" wrapper
 *    01/24/93 camilleg #ifdef NeXT ANSI-style function definitions
 *    01/23/94 jwb	install opsndfile(), getsamp(), FILDES declarations
 *    11/08/95 jwb	Change some function args back from double to float
 *			to avoid compilation problem.
 *    11/12/95 jwb	Insert PI #define.
 *			ANSI-style declarations for mono(), stereo(), etc.
 *			Remove #ifdef NeXT for ANSI-style unit generator
 *			function defs and comment out traditional function defs.
 *    12/15/96 jwb      Add filename to FILDES structure definition.
 *    04/24/00 jwb      Replaced old-c with OLD_C for linux compatibility.
 *    08/08/02 mertbay  Added fildes.filendian type to store the byte order 
 *			of the file and added byteswap2(),byteswap4(),
 *			byte_order() to handle byteorder problem. 
 *    09/08/06 jwb      Add include <sys/types.h> and <unistd.h> for Mac OS X
 *                    	installation.
 *****************************************************************************/

#include <sys/types.h>					    /* jwb 09/08/06 */
#include <unistd.h>					    /* jwb 09/08/06 */
#include <math.h>				


/******************************************************************************
 * 			Defines and Macros
 *****************************************************************************/

#define		MONO		1		    /* One channel of sound */
#define		STEREO		2		    /* 	       Two channels */
#define		QUAD		4		    /*        Four channels */
#define		R_IBASE		4		    /* P_array base - Icard */
#define		R_FBASE		3		    /* P_array base - Fcard */
#define		NLPTABL2	256.	/* nonlinear table center index     */
#define		NLPTABLEN	513	/* nonlinear table length           */
#define		NHMAX		100     /* no. of harmonics for nlf instru */
#define		NHM2		(NHMAX/2)
#define		NABMAX		50	/* size of Alvbr table */
#define		TERM		9999.	/* termination of float no. list */
#define 	DUR 		(card.I.dur)        /* I card duration  */
#define 	STIME 		(card.I.stime)      /* note start time) */
#define 	P 		printf		    /* jwb 10/19/93 */
#define 	PI   3.14159265358979323846	    /* jwb 11/12/95 */

/* Map P-array onto a structure definition */
#define Istencil(sname,base)  ((sname *)(P_buf + R_IBASE + base))
#define Fstencil(sname,base)  ((sname *)(P_buf + R_FBASE + base))

/* Declare a voice to m4c */
#define maptyp(nam,max,es,set,or,endnote)  \
	((es *)_maptyp((nam),(max),sizeof(es),(set),(or),(endnote)))

/* Declare a function to m4c				Talmage 11/21/87 */
#define maptypF(nam,es,set)  \
	((es *)_maptypF((nam),sizeof(es),(set)))

/* get memory for tables, etc */
#define getmem(type,num) ( (type *)calloc( (num),sizeof(type)) )

/* free memory for tables, etc. */
#define freemem(fptr) ( cfree(fptr) )

/* The lesser of two numbers */
#define min(x,y)  (((x) < (y)) ? (x) : (y))

/* The greater of two numbers */
#define max(x,y) (((x) > (y)) ? (x) : (y))

/*  A floating-point modulo macro */
#define amod(x,y) ((x) - ((int)((x)/(y)))*(y))

/* Absolute value of a number */
/*  #define abs(x)  (((x) > 0) ? (x) : (-(x)))  */

/* Call mono(), stereo(), or exit() as appropriate from a sample function */
/* Camille Goudeseune, 11/91 */
#define NoQuad(x, lprop) \
	switch (Nchans) { \
	case MONO:	mono(x); break; \
	case STEREO:	stereo(x, lprop); break; \
	default:	printf("m4c: instrument error: Nchans != 1 or 2\n"); \
			exit(-1); \
	}

/* for .init section of an instrument template file */
/* Camille Goudeseune, 11/91 */
#define InitDone() \
	printf("Complete init of %d instances of %s.\n", maxv, name);\
	fflush(stdout)

/******************************************************************************
 * 			Structure Typedefs
 *****************************************************************************/

typedef  struct _aicard		I_SCORE;    	/* Instrument cards */
typedef  struct _afcard		F_SCORE;    	/*   Function cards */
typedef  struct _aocard		O_SCORE;    	/*           -Icard */
typedef  union  _score		SCORE;	    	/*  Score structure */
typedef  struct _partials       PARTIALS;	/* formnt: partials (Maher)  */
typedef  struct _formnt         FORMNT;      	/* formnt: state    (Maher)  */
typedef  struct _linens		LINENS;	    	/* State for linens */
typedef  struct _adsr           ADSR;        	/* ADSR envelope    (Maher)  */
typedef  struct _env3           ENV3;        	/* ENV3 envelope    (Maher)  */
typedef  struct _rstate       RSTATE;        	/* randh,randi structure     */
typedef  struct _bstate       BSTATE;        	/*     balance structure     */
typedef  struct _tone           TONE;        	/*     tone structure        */
typedef  struct _fstate       FSTATE;        	/*     filter structure      */
typedef  struct _fparms       FPARMS;        	/*     filter structure      */
typedef  struct _envp           ENVP;        	/*     envelop structure     */
typedef struct 				     	/*  structure used for       */
{						/* opsndfile() and getsamp() */
  int unit, ftype, headsize, nobytes, 		/* sound file input programs */
      sammax, sr, nchans, filendiantype;
  char* filename;					    /* jwb 12/15/96 */
} FILDES;

/******************************************************************************
 * This structure holds the state of 
 * LINENS, the envelope generator
 *****************************************************************************/
struct _linens
{
    int segment;
    int att_samp;
    int dur_samp;
    int dec_samp;
    float att_incr;
    float dec_incr;
    float output;
};

/******************************************************************************
*  specification of partial characteristics for formnt() (Maher)
******************************************************************************/
struct _partials
{
   float freq;
   float coef;
   float phase;
};

/******************************************************************************
*  specification of state for formnt()                  (Maher)
******************************************************************************/
struct _formnt
{
   float phase;
   float ss;
   float scale;
};

/******************************************************************************
*  specification of ADSR structure  for adsr()          (Maher)
******************************************************************************/
struct _adsr
{
  int seg;        int tss;
  float a_incr;   float d_incr;   float r_fac;
  float sus;      float out;
};

/******************************************************************************
 *	State structure for ADeSR envelope for adesr()   (jwb, kelley 6/91)
 *****************************************************************************/

typedef struct
{
  float a_incr, d_incr, s_fac, r_fac, out;
  int seg, tas, tds, tss;
}
 ADeSR;

/******************************************************************************
 	State structure for for comb() and allpass() ugens  (jwb)
******************************************************************************/

typedef struct
{
     float Delay, Rvt, Atten, Atten2;
     int l, Maxnum;
     /* below is a pointer to an array */
     float * Delbuf;
} COMBSTATE;

/******************************************************************************
   State structure for reverb() unit generator   (jwb)
******************************************************************************/
typedef struct
{
  int ncombs;
  int naps;
  float *delay;
  COMBSTATE **cstates;
} REVSTATE;

/******************************************************************************
    structure for three segment envelope generator  (Beauchamp 7/10/86)
    modified by Maher  10/22/86
******************************************************************************/
struct _env3 
{
  int seg;  double ephase, incr[3], incr1[3], incr2[3], end[3];
};

/******************************************************************************
*  specification of RSTATE structure for randh, randi   (Maher)
******************************************************************************/
struct _rstate
{
	   float rphase;
	   float starta;
	   float stopa;
	   float delta;
	  };
/******************************************************************************
*  specification of TONE structure for tone()      (Maher)
******************************************************************************/
struct _tone
{
    float tone[3];  
  };
/******************************************************************************
*  specification of BSTATE structure for balance()      (Maher)
******************************************************************************/
struct _bstate
{
    TONE state1;
    TONE state2;
  };
/******************************************************************************
*  specification of FSTATE structure for filters      (Maher)
******************************************************************************/
struct _fstate
{
    float fstate[4];  
  };
/******************************************************************************
*  specification of FPARMS structure for filters      (Maher)
******************************************************************************/
struct _fparms
{
    float fparms[10];  
  };
/******************************************************************************
*  specification of ENVP structure for envelop() (Maher)
******************************************************************************/
struct _envp
{
    float envp[6];  
  };

/******************************************************************************
 * 		Structures returned by getscore()
 *****************************************************************************/
struct _aicard
{
	int type;		    /*     Instrument type-number assigned */
	int which;		    /* Instrument Instance number assigned */
	float stime;		    /* 		 Note start-time (seconds) */
	float dur;		    /* 		   Note duration (seconds) */
	float *p;		    /* 	         Floating-point parameters */
	int parc;	            /* 		      Number of parameters */
};

struct _afcard						/* spt 11/21/87 */
{
	int fnum;		    /* 	 	  Function number assigned */
	float stime;		    /* 		   	   Generation time */
	float *p;		    /* 	 	 Floating-point parameters */
	int parc;	            /* 		      Number of parameters */
};

struct _aocard
{
	int type;		    /*     Instrument type-number assigned */
	int which;		    /* Instrument Instance number assigned */
	float etime;		    /* 		 Note start-time (seconds) */
};

union _score
{
	I_SCORE I;			     /* Use this if  Icard */
	F_SCORE F;			     /* Use this if  Fcard */
	O_SCORE O;			     /* Use this if -Icard */
};


/******************************************************************************
    		Unit Generator Function Declarations (ANSI-style)
******************************************************************************/
float adesr(float, ADeSR *);
int adesr_set(float, float, float, float, float, float, float, ADeSR *);
float adsr(float, ADSR *);
int   adsr_set(float, float, float, float, float, float, ADSR *);
float allpass(double, double, COMBSTATE *);
float *alvbr(void);
float balance(double, double, BSTATE *);
int   balset(double, double, BSTATE *);
float buzz(float, float, int, float *, float *);
float bpf2(double, FSTATE *, FPARMS *);
float bsf2(double, FSTATE *, FPARMS *);
void byteswap2(short *);                           /*mb 08/08/02*/
void byteswap4(int *);                             /*mb 08/08/02*/
int byte_order();                                  /*mb 08/08/02*/
float comb(double, double, COMBSTATE *);
void  combset(double, int, COMBSTATE *);
float env3(float, float *, register ENV3 *);
float envelop(float, float *, float *, ENVP *);
int   evpset(float, float, float, float, ENVP *);
double expon(double *, double *);
int   expset(double *, double *, float);
double filhp(double, double, double);
float formnt(float, float, float *, PARTIALS *, float *, float *, FORMNT *);
int   fourfun(float **, float *, int, int);
double *getacoeff(void);
float getmaxnlf(float *, float *, float *, double, double *, int, double, 
		double, double);
float *getsamp(FILDES *, int, int);
float hpf1(double, FSTATE *, FPARMS *);
float hpf2(double, FSTATE *, FPARMS *);
int   krand(void);
float linens(float, LINENS *);
int   linsegfun(float **, float *, int, int);
int   linset(float, float, float, LINENS *);
float lpf1(double, FSTATE *, FPARMS *);
float lpf2(double, FSTATE *, FPARMS *);
int   mono(float);
float nlf(float, float, float, float *, float *, float *, FSTATE *, FPARMS *);
void  nlpfun(float **, int, double, double, double, double, float *);
float noise(float, float *);
void  opsndfile(char *, FILDES *);
float oscil(float, float, float *, float *);
float oscil1(float, float, float *, float *);
float oscili(float, float, float *, float *);
int   output(float, float, float, float);
float phasmod(float, float, float*, float*, float, float, float*, float*);
float rand01(float);
float randh(float, float, RSTATE *);
float randi(float, float, RSTATE *);
void  resonset(double, double, double, int, FSTATE *, FPARMS *);
void  reverbinit(int, int, float *, REVSTATE *);
void  reverbset(int, REVSTATE *);
float reverb(float, float, REVSTATE *);
int   setenv3(float, float, float, float *, ENV3 *);
float *sintab(int);
void  skrand(unsigned int);
float slope(float *, float *);
int   slopeset(float *, float *, float);
int   splinsegfun(float **, float *, int, int);
int   stereo(float, float);
int   stereoManual(float, float);
float tone(double, TONE *);
int   tonset(double, double, TONE *);
float varbpf(double, double, double, double, FSTATE *, FPARMS *);
int   varbpfset(double, double, double, double, double, double, FPARMS *);
float vfmult(double, double, float *);
float voctave(float, float, float, float *, float *, float *, float *);
float vresnorm(double, double);
float vreson(double, double, double, float *, FSTATE *, FPARMS *);
void  vresonset(double, int, FSTATE *, FPARMS *);
float reson(double, FSTATE *, FPARMS *);
/* conversion functions in U_convert.c: */
float cpsoct(float);
float cpspitch(float);
double frfacpitch(float);
float octcps(float);
float octpitch(float);
float pitchcps(float);
float pitchoct(float);
float sicps(float);
float sioct(float);
float sipitch(float);
float siper(float);
/* special functions: */
/* defined in pass3.c: */
void TacetNote();	/* PLIST isn't exported, so no function prototype */
int OrchNil(void *);
/* defined in pass1.c for handling scorecard strings: */
char *GetScorecardString(double);

/******************************************************************************
    		Unit Generator Function Declarations (traditional)
******************************************************************************/
#ifdef OLD_C					    /* jwb 04/24/00 */

float adesr();				/*	adesr envelope generator */
int adesr_set();                        /*    setup for adesr gen       */
float adsr();				/*      adsr envelope generator */
int   adsr_set();                       /*           setup for adsr gen */
float allpass();			/*	all pass reverberator   */
float balance();                        /*    balance processed val     */
int   balset();                         /*        setup for balance gen */
float buzz();                           /*    band-limited pulse source */
float bpf2();				/*    2nd order bandpass filter */
float bsf2();				/*    2nd order bandstop filter */
float comb();				/*   comb filter reverberator   */
void combset();				/* initialize comb or allpass reverb */
float env3();                           /*       3 segment envelope gen */
float envelop();			/*    table lookup envelope gen */
int   evpset();                         /*        setup for envelop gen */
double expon();				/*     exponential envelope gen */
int   expset();                         /*          setup for expon gen */
float formnt();				/*   formant additive synthesis */ 
float hpf1();				/*    1st order highpass filter */
float hpf2();				/*    2nd order highpass filter */
double filhp();				/* calculate hpf2 response      */
float linens();				/*    linear envelope generator */
int   linset();                         /*         setup for linens gen */
float lpf1();				/*     1st order lowpass filter */
float lpf2();				/*     2nd order lowpass filter */
int mono();				/*  monaural sample output gen  */
float noise();				/*       random noise generator */
float oscil();                          /*      non-interpolating oscil */
float oscil1();                         /*            single-pass oscil */ 
float oscili();				/*     Interpolating Oscillator */
int output();				/* 4-channel sample output gen  */
float phasmod();			/*   phase modulation generator */
float rand01();                         /*    0 -> 1. random number gen */
float randh();                          /*         random staircase gen */
float randi();                          /*     random linear interp gen */
void  resonset();                       /*         setup of reson gen   */
void  reverbset();			/* setup for reverb()           */
float reverb();				/*  reverberator unit gen       */
int   setenv3();                        /*           setup for env3 gen */
float slope();				/*     linear segment generator */
int   slopeset();                       /*          setup for slope gen */
int stereo();				/* stereo sample output gen     */
int stereoManual();			/* stereo sample output gen     */
float tone();				/* 1st order treble/bass filter */
int   tonset();                         /*           setup for tone gen */
float varbpf();                         /*     variable bandpass filter */
int   varbpfset();                      /*         setup for varbpf gen */
float vfmult();                         /*            amp*func[residue] */
float voctave();                        /*       vibrato/glissando unit */
float vresnorm();                       /*       bandwidth compensation */
float vreson();                         /*     variable bandpass filter */
void  vresonset();                      /*         setup for vreson gen */
float reson();				/*     fixed bandpass filter    */
int    linsegfun();                     /* allocate, generate lseg tabl */
int    splinsegfun();                   /* allocate, generate lseg tabl */
float nlf();                            /*     nlf processor            */
int nlpfun();				/* alloc., gen. nonlin. poly table */
double *getacoeff();		        /* get ptr to nlf acoef array   */
float getmaxnlf();			/* predicter of max nlf amplitude  */
float *alvbr();				/* get ptr to nlf alvbr array */
int    fourfun();			/* alloc., gen. harm sum table  */
float *sintab();			/* allocate, generate sine table*/ 
float cpsoct();				/*            octaves ==> Hertz */
float cpspitch();			/*             oct.pc ==> Hertz */
double frfacpitch();			/* delta oct.pc ==> freq factor */
float octcps();				/*            Hertz ==> octaves */
float octpitch();			/*           oct.pc ==> octaves */
float pitchcps();			/*             Hertz ==> oct.pc */
float pitchoct();			/*           octaves ==> oct.pc */
float sicps();				/*   Hertz ==> sample increment */
float sioct();				/* octaves ==> sample increment */
float sipitch();			/*  oct.pc ==> sample increment */
float siper();				/* seconds ==> sample increment */
void skrand();				/* set seed for krand() */
int krand();				/* random number in [0,32767] */
void TacetNote();			/* to silence the rest of a note */
int OrchNil();				/* do-nothing .sample function */
char *GetScorecardString();		/* convert scorecard float to string */
void opsndfile();			/* open sound file */
float *getsamp();			/* get samples from sound file */

#endif

/******************************************************************************
 * 			Global Variables
 *****************************************************************************/

extern SCORE card;			   /* Structure from score data */
extern float *P_buf;			   /*          Parameter-buffer */
extern float **Func;                       /*  Function table pointers  */
extern double **DFunc;                     /*  Function table pointers  */
extern float SR;			   /*         The Sampling Rate */
extern int Nchans;			   /*  How many channels to use */
extern int FUN_LENGTH;			   /* Length of table used by ugens */
extern float FFUN_LEN;			   /* Length of table used by ugens */
extern float *Sine;			   /* pointer to a sine table  */

#endif 
