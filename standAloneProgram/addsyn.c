/******************************************************************************
 *		addsyn.c	latest edit: 06/19/08
 *
 *	Time-variant additive synthesis from analysis file data
 *	using linear interpolation of original data points.
 *
 *	For synthesis: fundamental freq., no. of harmonics, and
 *	               duration are read from analysis file header.
 *
 *	execution:  addsyn an_file syn_file.SUF
 *		or  addsyn an_file syn_file.SUF >& list.file &
 *
 *	        where SUF = .snd, .au, .fp., .sh, or .wav.
 *
 *	compile:  cc -O -o addsyn addsyn.o anread.o header.o gettime.o \
				byteorder.o sndhdr.o wavhdr.o -lm
 *
 *      09/11/87	Programmed by  James W. Beauchamp, Univ. of Ill. U-C
 *	05/17/91 jwb	NeXT installation: replace TABSIZ. by TABSIZF
 *	05/17/91 jwb	Remove calloc() declaration
 *	05/17/91 jwb	Install gettime() to give beg. & end real ex. times
 *	05/22/91 jwb	Modify to output next sound files; install wrNxtHdr()
 *	11/02/91 jwb	Complete modification to handle .snd, .sh, and .fp
 *			  files; install getfiltype(); for .snd files use
 *			  fs = 22050, and tailor nhar accordingly.
 *	10/18/93 jwb	Replace fmod by amod to speed up add synth loop
 *	04/23/94 jwb	Remove anread() and link with separate anread.o
 *	08/17/94 jwb	Port for DEC Alpha requires byte swapping
 *	12/05/94 jwb	Install time counter printout
 *	09/17-22/96 jwb	Adjust the starting phase of the additive synthesis.
 *	11/29/97 zheng	Change #ifdef alpha to if (byte_reverse)
 *      01/29/98 jjjm   'byte_order()' now returns a value (endianness)
 *      02/05/98 jjjm   Moved wrNxtHdr() to sndhdr.c as writeSndHdr()
 *      02/10/98 jjjm   Modify to handle WAVE files.  'byte_reverse' is now
 *                      set in anread() for input and passed to addsyn() as
 *                      a parameter for output samples so that input and output
 *                      byte orderings are handled separately.
 *	03/09/98 jwb	Add ".au" suffix option for output sound files.
 *	03/28/98 jjjm   Add a call to 'fixSndHdr()' after samples written.
 *	04/02/98 jwb	Revert AU filetype to SND filetype.
 *	04/15/98 jwb	Print message that .au is converted to .snd.
 *	04/18/98 jwb	Change NEXT to SND.
 *	04/19/98 jwb    addsyn(): Rearrange some things. Install function
 *			getout().
 *	04/25/98 jwb	addsyn(): Introduce clipping to +-32767 before output.
 *	04/30/98 jwb	addsyn(): Print warning when clipping occurs.
 *	12/06/99 jwb	addsyn(): Add 2 frames at beginning. Fix dfr index.
 *	02/02/00 jwb	addsyn(): Bug fix: narg0++
 *	03/01/00 jwb	addsyn(): Remove fade-in from first two frames.
 *				  Install output file clobber protection.
 *	03/02/00 jwb	addsyn(): Add two frames at the end.
 *	06/19/08 jwb	addsyn(): Install trap for "not a number" sample output
 *****************************************************************************/
/*    includes and defines   */
#include <stdlib.h>					    /* jwb 11/02/91 */
#include <math.h>
#include <stdio.h>
#include <strings.h>					    /* jwb 06/19/08 */
#include "macro.h"
#include "header.h"
#include "byteorder.h"					    /* jjjm 01/29/98 */
#include "sndhdr.h"					    /* jjjm 02/05/98 */
#include "wavhdr.h"					    /* jjjm 02/10/98 */
#define TABSIZ 5000
#define TABSIZF 5000.
#define P printf
#define FLOAT 0						    /* jwb 11/02/91 */
#define SHORT 1						    /* jwb 11/02/91 */
#define SND   2						    /* jwb 04/18/98 */
#define WAVE 3						   /* jjjm 02/10/98 */
#define AU 4						    /* jwb 03/09/98 */
#define MAXFTYPES 5					    /* jwb 03/09 98 */
char *tail[MAXFTYPES] = {".fp",".sh","snd","wav",".au"};    /* jwb 03/09/98 */

/*  function declarations */
char *gettime(); int getfiltype();

/*    global variables:     */
HEADER header;
int nhar, nhar1, npts;
float *cmag, *dfr, *phase, dt, tl, smax, fs, fa, scalefac;

main(argc,argv) int argc; char **argv;
{
  char *anfil, *synfil;

  /* Assume that FLOAT and SHORT files are in the machine's native byte order
   * unless the '-x' command-line option is specified.
   * (For SND or WAVE files, 'byte_reverse' is set appropriately below.)*/
  int byte_reverse = 0;			       /* zheng 10/29/97, rev. jjjm */

  if(argc < 3) {P("Usage: %s [-x] anfile synfile\n%s%s\n",
		  argv[0],"       where '-x' means exchange byte order for ",
		  "'.fp' or '.sh' files."); exit(0);}  /* rev.jjjm 02/11/98 */

  /* Assume the input data has the same byte order as the host architecture.
     Exchange the byte order if the '-x' command-line flag is specified.
     (The input byte order may also be determined by the type of the input
     file -- see handling of NeXT '.snd' files and WAVE files below.)
                                                            jjjm 02/11/98 */
  if ( *argv[1] == '-' ) 	/* check for option flag */
  {
    if ( argv[1][1] == 'x' )	/* '-x' flag means exchange bytes */
      byte_reverse = 1;
    else
    {
      fprintf(stderr, "%s: unknown option flag", argv[0]);
      exit(-1);
    }
    /* shift argument pointer past option flag */
    argc--;
    argv++;
  }

/*     figure out input and output file names  */
  anfil = *(argv + 1);  synfil = *(argv + 2);

/*  read in analysis data for given file  */
  if(anread(anfil))
  {
    P("analysis file %s isn't available; try another\n",anfil);
    exit(-1);
  }
  P("Analysis data has been read in\n");

/*  resynthesize tone at sample rate fs */
  P("Begin synthesis\n");
  addsyn(synfil, byte_reverse);				    /* jwb 12/06/99 */
  P("\nSynthesis of file %s complete \n", 		    /* jwb 04/02/98 */
			synfil);		    /* jwb 04/02/98 */
} /* end main() */

addsyn(filnam,byte_reverse) char *filnam; int byte_reverse; /* jwb 12/06/99 */
{
  int uno,i,imax,j,L;
  int filetype, samptype;
  int samp;						    /* jwb 04/30/98 */
  register int k,narg0,narg1,narg2;  			    /* jwb 12/06/99 */
  short int isum; 					    /* jwb 11/02/91 */
  float fac,sifac,t,st,sum,si,fj,frac1,frac2,ckj,*costab;   /* jwb 09/17/96 */
  float tm, tmp=1.;                                         /* jwb 12/05/94 */
  float tfade, fade, dfade, tend;			    /* jwb 03/01/00 */
  register float phasek;

/*  open output file  */
  if((uno = open(filnam,0)) != -1)                         /* jwb 03/01/00 */
  { 							    /* jwb 03/01/00 */
    P("File %s already exists! Abort\n\n",filnam);	    /* jwb 03/01/00 */
    exit(-1);						    /* jwb 03/01/00 */
  }							    /* jwb 03/01/00 */
  uno = creat(filnam,0644);
  if(uno == -1)
  {
    P("addsyn output file %s cannot be created\n",filnam);
    exit(-1);
  }
/*  determine output file type  */
  filetype = getfiltype(filnam);			    /* jwb 11/02/91 */
  if(filetype == MAXFTYPES) getout();			    /* jwb 04/19/98 */
  fs = header.sr;					    /* jwb 04/19/98 */
  if(filetype == AU)                                        /* jwb 04/15/98 */
  {                                                         /* jwb 04/15/98 */
    P("\nThis program does not support the .au extension.\n");
    P("Extension .snd will be used instead\n");             /* jwb 04/15/98 */
    filetype == SND;                                        /* jwb 04/18/98 */
  }                                                         /* jwb 04/15/98 */
  if(filetype == SND)				     	    /* jwb 04/18/98 */
  {
/*  create next header */
    writeSndHdr(uno, (int)fs, 1, 0, SND_FORMAT_LINEAR_16);  /* jwb 11/02/91 */
				                            /* rev. jjjm 98 */
    samptype = SHORT;					    /* jwb 04/19/98 */

    /*
     * '.snd' data is big-endian, so if the host is little-endian,
     * byte-reverse the data. [This is not strictly necessary here,
     * since 'byte_reverse' is set correctly in main() above for reading
     * the analysis data, which has the same endianness as '.snd' data.
     * However, this is added in case this function is reused in another
     * context in which 'byte_reverse' was not set for bit-endian data.]
     */
    byte_reverse = (byte_order() != big_endian);           /* jjjm 02/10/98 */
  }
  else if (filetype == WAVE) 				   /* jjjm 02/10/98 */
  {
    writeWavHdr(uno, (int) fs, 1, 0, 16);

    /* The WAVE format doesn't support float samples. */
    samptype = SHORT;

    /* WAVE data is little-endian, so if the host is big-endian,
       byte-reverse the data. */
    byte_reverse = (byte_order() != little_endian);	   /* jjjm 02/10/98 */
  }
  else samptype = filetype;				    /* jwb 08/17/94 */

/* lower no. of harmonics if necessary to fit the sample rate  */
  nhar = min(.5*fs/fa,nhar);	 		    	    /* jwb 04/19/98 */
  P("Synthesis sample frequency = %.0f\n", fs);	    	    /* jwb 11/02/91 */
  P("Synthesis number of harmonics = %d\n", nhar);  	    /* jwb 11/02/91 */
  fflush(stdout);

/*  store cosine table  */
  fac = 8.*atan(1.)/(float)TABSIZ;			    /* jwb 09/17/96 */
  costab = (float *)calloc(TABSIZ+1,sizeof(float));
  for(L=0;L<=TABSIZ;L++) costab[L] = cos(fac*L);

/*  store initial phase values  */
  for(k=1; k <= nhar; k++) phase[k] = phase[k]/fac;	    /* jwb 04/23/94 */

  P("\nOriginal time in input file:  ");		    /* jwb 12/05/94 */
/*  count time from -2.*dt to tl  */
  sifac = TABSIZ/fs; st = 1./fs;			    /* jwb 09/22/96 */
  tfade = 2.*dt; fade = 0.; dfade = st/tfade;		    /* jwb 12/06/99 */
  tend = tl + tfade -st;				    /* jwb 03/01/00 */
  for(t=-tfade,samp=0; t < tend; t += st,samp++) 	    /* jwb 03/01/00 */
  {
    sum = 0.;						    /* jwb 12/06/99 */
   /* for 1st tfade*fs samples do synthesis based on first frame */
    if(t<0.)						    /* jwb 12/06/99 */
    {							    /* jwb 12/06/99 */
      for(k=1; k <= nhar; k++)				    /* jwb 12/06/99 */
      {							    /* jwb 12/06/99 */
        phasek = amod(phase[k], TABSIZ);	  	    /* jwb 12/06/99 */
        if(phasek < 0.) phasek += TABSIZ;		    /* jwb 12/06/99 */
        sum += cmag[k]*costab[(int)phasek];		    /* jwb 03/01/00 */
        phase[k] = phasek + sifac*k*fa;			    /* jwb 03/01/00 */
      }							    /* jwb 12/06/99 */
      fade += dfade;					    /* jwb 12/06/99 */
    } /* end if */					    /* jwb 12/06/99 */
    else if((t>=0.)&&(t<tl))				    /* jwb 03/01/00 */
    {							    /* jwb 12/06/99 */
	/*  find frame number j */
      fj = t/dt; j = fj;
      frac2 = fj - j;  frac1 = 1. - frac2;
      narg0 = (j-1)*nhar1;				    /* jwb 12/06/99 */
      narg1 = j*nhar1;  narg2 = (j+1)*nhar1;
      if((j+1) > (npts-1)) continue;			    /* jwb 03/01/00 */

      for(k=1; k <= nhar; k++)
      {
        narg0++; narg1++; narg2++;			    /* jwb 02/02/00 */
        phasek = amod(phase[k], TABSIZ);	  	    /* jwb 09/17/96 */
        if(phasek < 0.) phasek += TABSIZ;
        ckj = frac1*cmag[narg1] + frac2*cmag[narg2];
        sum += ckj*costab[(int)phasek];

    /* compute and save phase for next sample time */
    /* use dfr one frame behind */
        if(j==0)					    /* jwb 12/06/99 */
          phase[k] = phasek + 				    /* jwb 12/06/99 */
            sifac*(k*fa + frac2*dfr[narg1]);		    /* jwb 12/06/99 */
        else /* j>=1 */					    /* jwb 12/06/99 */
          phase[k] = phasek +
            sifac*(k*fa+frac1*dfr[narg0]+frac2*dfr[narg1]); /* jwb 12/06/99 */
      } /* end for */
    } /* end else if */					    /* jwb 12/06/99 */
/* play out last 2 frames based on npts -1 data */
    else /* t>=tl */					    /* jwb 03/01/00 */
    {							    /* jwb 03/01/00 */
      narg0 = (npts-1)*nhar1;				    /* jwb 03/02/00 */
      for(k=1; k <= nhar; k++)
      {
        narg0++; 					    /* jwb 03/01/00 */
        phasek = amod(phase[k], TABSIZ);	  	    /* jwb 03/01/00 */
        if(phasek < 0.) phasek += TABSIZ;		    /* jwb 03/01/00 */
        sum += cmag[narg0]*costab[(int)phasek];		    /* jwb 03/01/00 */
    /* compute and save phase for next sample time */
        phase[k] = phasek + sifac*k*fa;			    /* jwb 03/01/00 */
      }
    } /* end else */					    /* jwb 03/01/00 */

    if(samptype == SHORT)
    {
      if(sum > 32767.) 					    /* jwb 04/25/98 */
      {							    /* jwb 04/30/98 */
        P("\namplitude out of range (%.0f) at sample %d ",  /* jwb 04/30/98 */
            sum, samp);					    /* jwb 04/30/98 */
        P("clipped to 32767.\n");			    /* jwb 04/30/98 */
        isum = 32767;			    		    /* jwb 04/25/98 */
      }							    /* jwb 04/30/98 */
      else if(sum < -32767.) 				    /* jwb 04/25/98 */
      {							    /* jwb 04/30/98 */
        P("\namplitude out of range (%.0f) at sample %d ",  /* jwb 04/30/98 */
            sum, samp);					    /* jwb 04/30/98 */
        P("clipped to -32767.\n");			    /* jwb 04/30/98 */
        isum = - 32767;		    			    /* jwb 04/25/98 */
      }							    /* jwb 04/30/98 */
      else isum = sum; 					    /* jwb 04/30/98 */
      if(isum==-32768) 					    /* jwb 06/19/08 */
      {							    /* jwb 06/19/08 */
        P("at sample no. %d, %f detected\n", samp, sum);    /* jwb 06/19/08 */
        isum = 0;					    /* jwb 06/19/08 */
        P("set sample value to %d\n", isum);		    /* jwb 06/19/08 */
      }							    /* jwb 06/19/08 */
      if (byte_reverse) byteswap2((short*)(&isum));	  /* zheng 11/29/97 */
      write(uno,&isum,sizeof(short));                       /* jwb 11/02/91 */
    } /* end if(samptype */
    else
    {
      if (byte_reverse) byteswap4((int*)(&sum));     	  /* zheng 11/29/97 */
      write(uno,&sum,sizeof(float));                        /* jwb 11/02/91 */
    }
    tm = amod(t, .1);                                       /* jwb 12/05/94 */
    if (tm < tmp) {P( "%.2f  ",t); fflush(stdout);}         /* jwb 12/05/94 */
    tmp = tm;					            /* jwb 12/05/94 */
  } /* end for(t -- sample generation loop */

  /* Adjust the file and data size values in the SND or WAVE header. */
  if (filetype == WAVE) 				   /* jjjm 02/10/98 */
    fixWavHdr(uno);					   /* jjjm 02/10/98 */
  else if (filetype == SND)                                 /* jwb 04/18/98 */
    fixSndHdr(uno);					   /* jjjm 02/10/98 */
  P("\n");						    /* jwb 03/01/00 */
} /* end addsyn() */

int getfiltype(name) char *name;
{ int i, len;
  len = strlen(name);
  for(i=0;i<MAXFTYPES;i++)
  {
    if(!strcmp(&name[len-3],tail[i])) return(i);
  }
  return(MAXFTYPES);
}

getout()						    /* jwb 04/19/98 */
{							    /* jwb 04/19/98 */
  int i;						    /* jwb 04/19/98 */
  P("Sorry, only ");					    /* jwb 04/19/98 */
  for(i=0;i<MAXFTYPES-1;i++) P("%s, ", tail[i]);	    /* jwb 04/19/98 */
  P("and %s extensions are supported\n", tail[MAXFTYPES-1]);/* jwb 04/19/98 */
  exit(-1);						    /* jwb 04/19/98 */
}							    /* jwb 04/19/98 */

int plabel(double xpos, double ypos)
{
    return 0;
}
