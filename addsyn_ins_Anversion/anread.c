/*****************************************************************************
		anread.c
	Routine to read in analysis file data into cmag and dfr arrays.
	Modified for use in M4C additive synthesis instrument.

	Programmer:  James Beauchamp, 1987
	All Rights Reserved.
        mb = Mert Bay
   Changes:
	12/14/93 jwb	Modify to return file pointer and critical analysis
			data.  Allocate new cmag, dfr arrays only if data
			has not been read in before.  Eliminate references
			to br and time.
	11/01/95 jwb	Update to read 'compact' files and insert byte swapping
			for reverse byte order machines (e.g., Dec Alpha)
	08/15/02 mb	Modify to swap bytes for little endian machines.
******************************************************************************/
/*    includes and defines   */
#include <strings.h>
#include "macro.h"
#include "header.c"

#define P printf
#define NO_ANFILES 100					/* jwb 12/14/93 */
#define SH_MAX 32768.					/* jwb 11/01/95 */
typedef struct						/* jwb 12/14/93 */
{
  char *anfile;
  int fptr, nhar1, npts;
  float fa, *phas0, *cmag, *dfr;
} ANDAT;						/* jwb 12/14/93 */

/*    global variables:     */
int nhar1, npts;
float fa, *phas0, *cmag, *dfr;
ANDAT andat[NO_ANFILES];				/* jwb 12/14/93 */
int nofilesread = 0;					/* jwb 12/14/93 */
float tl, dt, smax;

int anread(char *anfile) 				/* jwb 12/14/93 */
{
  int i,nbpb,nbph, k,k1,k1inc, narg;
  int unit, nhar, type, curfileno;			/* jwb 12/14/93 */
  int j, nhar2, fftlen;                                 /* jwb 11/01/95 */
  unsigned short int *tempbufs;				/* jwb 11/01/95 */
  float *tempbuf, cm, df, dfscale;			/* jwb 11/01/95 */
  float ampscale;					/* jwb 11/01/95 */
  HEADER header;

/*	determine whether anfile data has been read yet  */
  curfileno = nofilesread;				/* jwb 12/14/93 */
  if(nofilesread > 0)					/* jwb 12/14/93 */
  {
    for(i=0;i<nofilesread;i++)				/* jwb 12/14/93 */
    {
      if(!strcmp(anfile,andat[i].anfile)) 		/* jwb 12/14/93 */
      {
	curfileno = i;
	break;
      }
    }
  }
  // if(curfileno < nofilesread)	/* file already read */ /* jwb 12/14/93 */
  // {
  //   P("File %s has already been read. Return same data.\n", anfile);
  //   nhar1 = andat[curfileno].nhar1;
  //   npts = andat[curfileno].npts;
  //   fa = andat[curfileno].fa;
  //   phas0 = andat[curfileno].phas0;
  //   cmag = andat[curfileno].cmag;
  //   dfr = andat[curfileno].dfr;
  //   return(andat[curfileno].fptr);
  // }
  if(curfileno >= NO_ANFILES)
  {
    P("Max. no. analysis files (%d) exceeded. Abort.\n", NO_ANFILES);
    return(-1);
  }
  if((unit = open(anfile,0)) == -1)  return(-1);  /* file not available */

  nofilesread++;  /* new file has been read, so increase cntr --jwb 12/14/93 */

  rdat(unit,&header);
  P("Read file %s \nDate Recorded: %s\n",anfile,header.date);
  P("Instrument is %s\n", header.instrument);		    /* jwb 11/01/95 */
  P("Analysis frequency = %.1f Hz\n", header.fa);	    /* jwb 11/01/95 */
  P("Time between frames = %.5f sec\n", header.dt);	    /* jwb 11/01/95 */
  P("Analysis format is type %s\n", header.type);	    /* jwb 11/01/95 */
  P("Number of analysis frames is %d\n", header.npts);
  P("Number of harmonics per frame = %d\n", header.nhar);
  P("Duration of sound is %5.3f\n",header.tl);
  fflush(stdout);
  dt = header.dt;
  fa = header.fa;
  tl = header.tl;
  npts = header.npts;
  nhar = header.nhar;
  nhar2 = 2*nhar;					    /* jwb 11/01/95 */
  fftlen = header.fftlen;                             	    /* jwb 11/01/95 */
  if(nhar <= 0){					    /* jwb 11/01/95 */
    P("No. of harmonics in file = %d. Cannot synthesize.\n", nhar);
    exit(1);						    /* jwb 11/01/95 */
  }
  nhar1 = nhar + 1;

  /*   allocate space for cmag and dfr arrays  */
  cmag = (float *) calloc(npts*nhar1,sizeof(float));
  dfr  = (float *) calloc(npts*nhar1,sizeof(float));

  if(!strcmp(header.type,"full"))       type = 'f';
  if(!strcmp(header.type,"compressed")) type = 'c';
  if(!strcmp(header.type,"simple"))     type = 's';
  if(!strcmp(header.type,"compact"))    type = 'p';	    /* jwb 11/01/95 */

  P("\nType : %s\n",header.type);
  if(type=='f')  /* "full" */
  {
    /* read in and save initial phase data    */
     phas0 = (float *) calloc(header.fftlen/2,sizeof(float));
     nbph = fftlen*sizeof(float)/2;			    /* jwb 11/01/95 */
     if(read(unit,phas0,nbph) != nbph)
     {P("Cannot read initial phases\n");  return(1);}
     if(byte_order())                                        /* mb  08/02/02*/
      for(i=0;i<fftlen/2;i++) byteswap4((int*)(phas0+i));    /* jwb 08/02/02*/

   /* allocate space to read one analysis frame  */
     nbpb = fftlen*sizeof(float);			    /* jwb 11/01/95 */
     tempbuf = (float *) malloc(nbpb);			    /* jwb 11/01/95 */
  }
    else if((type=='c')||(type=='s')||(type=='p'))	    /* jwb 11/01/95 */
	        /* "compressed" or "simple" or "compact" data type */
  {
    /* read in and save initial phase data    */
     phas0 = (float *) calloc(nhar1,sizeof(float));
     nbph = nhar*sizeof(float);
     if(read(unit,phas0+1,nbph) != nbph)
     { P("Cannot read initial phases\n");  return(1);}
    if(byte_order())                                        /* mb  08/12/02*/
      for(j=1;j<=nhar;j++) byteswap4((int*)(phas0+j));      /* mb  08/12/02*/

   /* allocate space to read one analysis frame  */
     nbpb = nhar2*sizeof(float);			    /* jwb 11/01/95 */
     tempbuf = (float *) malloc(nbpb);			    /* jwb 11/01/95 */
  }
  else {P("Data type not recognized.  Cannot process."); return(1);}

  if((type=='f')||(type=='c'))  /*  "full" or "compressed"  */
     ampscale = 1./(header.fftlen*0.5*.54);
  else if(type=='p') ampscale = header.smax;		    /* jwb 11/01/95 */
  else ampscale = 1.;					    /* jwb 11/01/95 */

  if(type=='p')  /* for short int files */		    /* jwb 11/01/95 */
  {
    tempbufs = (unsigned short *)tempbuf;		    /* jwb 11/01/95 */
    nbpb /= 2;						    /* jwb 11/01/95 */
    if(fa > 0) dfscale = fa/SH_MAX;			    /* jwb 11/01/95 */
    else if(fa == 0) dfscale = 1./8.;			    /* jwb 11/01/95 */
  }

    /* read in analysis data for npts frames and put in cmag & dfr arrays */
  for(i=0;i<npts;i++)
  {
  /* read one frame of data */
    if(read(unit,tempbuf,nbpb) != nbpb)
    { P("Cannot read frame %d\n",i);
      npts = i-2;
      break;
    }
   if(byte_order())                                          /*mb  08/02/02*/
    {
     if(type=='p') /* compact */
      for(j=0;j<nhar2;j++) byteswap2((short*)(tempbufs+j)); /* jwb 11/01/95 */
    else /* full, compressed, or simple */
      for(j=0;j<nhar2;j++) byteswap4((int*)(tempbuf+j));    /* jwb 11/01/95 */
    }
    if(type=='f'){k1=4; k1inc = 4;} /* for type "full"  */
    else  {k1=0; k1inc = 2;} /* for type "compressed" or "simple"  */
    narg = i*nhar1;
    for(k=1;k<=nhar;k++,k1 += k1inc)
    {
      narg++;
      if(type!='p')    /* NOT compact format */		    /* jwb 11/01/95 */
      {
        cm = ampscale*tempbuf[k1];
        df = tempbuf[k1+1];
      }
      else             /* compact format */		    /* jwb 11/01/95 */
      {
        cm = ampscale*tempbufs[k1];			    /* jwb 11/01/95 */
        df = dfscale*k*(tempbufs[k1+1]-SH_MAX);		    /* jwb 11/01/95 */
      }
      cmag[narg] = cm;					    /* jwb 11/01/95 */
      dfr[narg]  = df;					    /* jwb 11/01/95 */
    }
  } /* end for(i=0;i<npts ... loop  */
  free(tempbuf);					    /* jwb 11/01/95 */
/* store analysis data for this file for subsequent calls to anread */
  andat[curfileno].nhar1 = nhar1;
  andat[curfileno].npts = npts;
  andat[curfileno].fa = fa;
  andat[curfileno].phas0 = phas0;
  andat[curfileno].cmag = cmag;
  andat[curfileno].dfr = dfr;
  andat[curfileno].fptr = unit;
  andat[curfileno].anfile = anfile;

  return(unit);
} /* end anread() */
