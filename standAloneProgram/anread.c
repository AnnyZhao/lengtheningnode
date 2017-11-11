/******************************************************************************
		anread.c	Latest edit: 02/11/98

	Routine to read in analysis file data into cmag and dfr arrays.
	Version for addsyn routines

	Programmer:  James Beauchamp
	All Rights Reserved
	Changes:
	04/28/94	jwb	Modify to read 'compact' analysis files
	08/17/94	jwb	Port for DEC Alpha requires byte swapping
 	11/29/97 zheng	Change #ifdef alpha to if (byte_reverse)
        01/29/98 jjjm   'byte_order()' now returns a value (endianness)
	02/11/98 jjjm   Fixed error message re: harmonics in file.
                        Made 'byte_reverse' local to anread() so its setting
                        affects only the byte ordering of the input file.
******************************************************************************/
#include <stdio.h>
#include "macro.h"
#include "header.h"
#include "byteorder.h"
#define P printf
#define SH_MAX 32768.

/*    global variables:     */
extern HEADER header;
extern int nhar, nhar1, npts;
extern float *cmag, *dfr, *phase, tl, dt, fa;

int anread(char *filnam)
{
    int i,nbpb,nbph, k,k1,k1inc, narg;
    int unit, type;
    int j, nhar2, fftlen;				    /* jwb 08/17/94 */
    unsigned short int *tempbufs;
    float *tempbuf, cm, df, dfscale;
    float ampscale;

    /* If the host is little-endian, byte-reverse the data. */
    int byte_reverse = byte_order();	 /* zheng 11/29/97, rev. jjjm 98 */

    if((unit = open(filnam,0)) == -1)  return(1);  /* file not available */

    rdat(unit,&header);

    P("Read file %s \n",filnam);
    P("Instrument is %s\n", header.instrument);
    P("Analysis frequency = %.1f Hz\n", header.fa);
    P("Time between frames = %.5f sec\n", header.dt);
    P("Analysis format is type %s\n", header.type);
    P("Number of analysis frames is %d\n", header.npts);
    P("Number of harmonics per frame = %d\n", header.nhar);
    P("Duration of sound is %.3f\n", header.tl);
    dt = header.dt; 
    fa = header.fa;
    tl = header.tl; 
    npts = header.npts;
    nhar = header.nhar;
    nhar2 = 2*nhar;					    /* jwb 08/17/94 */
    fftlen = header.fftlen;				    /* jwb 08/17/94 */
    if(nhar <= 0){
      P("No. of harmonics in file = %d. Cannot synthesize.\n",nhar); /*02/98*/
      exit(1);
    }
    nhar1 = nhar + 1; 

    cmag = (float *) calloc(npts*nhar1,sizeof(float));
    dfr  = (float *) calloc(npts*nhar1,sizeof(float));

    if(!strcmp(header.type,"full")) 	  type = 'f';
    if(!strcmp(header.type,"compressed")) type = 'c';
    if(!strcmp(header.type,"simple")) 	  type = 's';
    if(!strcmp(header.type,"compact")) 	  type = 'p';

    if(type=='f')  /* "full" */
    {
      /* read in and save initial phase data    */
       phase = (float *) calloc(fftlen/2,sizeof(float));
       nbph  = fftlen*sizeof(float)/2;
       if(read(unit,phase,nbph) != nbph) 
       {P("Cannot read initial phases\n");  return(1);}
       if (byte_reverse)				   /* zheng 11/29/97 */
         for(i=0;i<fftlen/2;i++) byteswap4((int*)(phase+i));  /* jwb 8/17/94 */
     /* allocate space to read one analysis frame  */
       nbpb    = fftlen*sizeof(float);
       tempbuf = (float *) malloc(nbpb); 
    }
    else if((type=='c')||(type=='s')||(type=='p'))  
		/* "compressed" or "simple" or "compact" data type */
    {
      /* read in and save initial phase data    */
      phase = (float *) calloc(nhar1,sizeof(float));
      nbph  = nhar*sizeof(float);
      if(read(unit,phase+1,nbph) != nbph) 
      { P("Cannot read initial phases\n");  return(1); }
      if (byte_reverse)					   /* zheng 11/29/97 */
        for(j=1;j<=nhar;j++) byteswap4((int*)(phase+j));   /* jwb 08/17/94 */
     /* allocate space to read one analysis frame */
      nbpb = nhar2*sizeof(float);			   /* jwb 08/17/94 */
      tempbuf = (float *) malloc(nbpb);	
    }
    else {P("Data type not recognized.  Cannot process."); return(1);}

    if((type=='f')||(type=='c'))  /*  "full" or "compressed"  */
       ampscale = 1./(fftlen*0.5*.54);
    else if(type=='p') ampscale = header.smax;
    else	ampscale = 1;
    
    if(type=='p')  /* for short int files */ 
    {
      tempbufs = (unsigned short *)tempbuf;
      nbpb /= 2;			
      if(fa > 0) dfscale = fa/SH_MAX;
      else if(fa == 0) dfscale = 1./8.;
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
      if (byte_reverse)					   /* zheng 11/29/97 */
     {if(type=='p') /* compact */
        for(j=0;j<nhar2;j++) byteswap2((short*)(tempbufs+j));
							     /* jwb 08/19/94 */
      else /* full, compressed, or simple */
        for(j=0;j<nhar2;j++) byteswap4((int*)(tempbuf+j));   /* jwb 08/17/94 */
     }
      if(type=='f'){k1=4; k1inc = 4;} /* for type "full"  */
      else  {k1=0; k1inc = 2;} /* for type "compressed" or "simple"  */
      narg = i*nhar1;
      for(k=1;k<=nhar;k++,k1 += k1inc)
      { 
	  narg++;
	  if(type!='p')    /* NOT compact format */
	  {
	    cm = ampscale*tempbuf[k1];
	    df = tempbuf[k1+1];	
	  }
	  else		   /* compact format */
	  {
	    cm = ampscale*tempbufs[k1];
	    df = dfscale*k*(tempbufs[k1+1]-SH_MAX);
	  }
	  cmag[narg] = cm;
          dfr[narg]  = df;
        } 
      } /* end for(i=0;i<npts ... loop */
    free(tempbuf);
    return(0);
}
