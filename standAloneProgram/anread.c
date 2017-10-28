/******************************************************************************
		anread.c

        latest edit: 03/09/15

	Routine to read in analysis file data into cmag and dfr arrays.
	Also computes rms (into cmag), br, and dfr arrays.

	Programmer:  James Beauchamp, 1987
	University of Illinois at Urbana-Champaign
        Copyright 1987 All Rights Reserved 

Changes:
  09/23/87  jwb	Modify for 'simple' rather than 'compressed' file 
		storage.
  06/27/92  jwb	Insert fabs() function for BR calculation.
  10/03/92  jwb	anwrite(): Insert code to guard against overwriting 
		existing file.
  08/26/93  jwb	anread(): If nhar eq. 0, don't compute rms, dfr.wt.av.; 
		instead, read in stored values.
		anwrite(): If nhar eq. 0, write rms and dfr.wt.av to 
		file.
  12/30/93  jwb	anwrite(): Accept default of nhar harmonics to write to
		file.	
  01/06-7/94 jwb anread(), anwrite(): Install 'compact' form of data 
		using short ints for file storage.
  02/24/94  jwb anread(), anwrite(): Modify to scale amplitudes 
		according the max. value of cmag array members.
  08/19/94  jwb	Port to Dec Alpha requires byte swap operations
		in both anread() and anwrite().
  10/12/94  jwb	Declare globals as static.
  10/24/94  jwb Change fabs to abs (much faster)
  11/11/94  jwb	Fix bug introduced on 8/19/94.
  01/04/95  jwb	anwrite(): Save npts and tl in header.
  04/19/95  jwb	anwrite(): Fix bug wrt Alpha writing of phase array.
  05/19/95  jwb	anwrite(): Exclude the possibility of negative cmag.
  11/05/95  jwb	anwrite(): Restrict nharw to less than sr/(2*fa).
  11/17/95  jwb	anread(): No. of harmonics question goes to stderr
  05/06/96  jwb	anread(): Add nhreq arg to function. if nhreq == 1 
		the program requests nhar from the user. Otherwise,
		it just takes nhar = header.nhar.
  09/24/96 jwb	anread(): Additional header listing print out.
  10/30/97 zheng Replace #ifdef __alpha with if (byte_reverse)
  11/05/97 zheng Use gets and sscanf for input
  01/31/98 jjjm  Revised byte_order() to return a value.
  02/21/99 jwb	Install code to enable read and storage of data in
		auxiliary arrays.
  03/28/01 jwb  anread() and anwrite(): Fix bugs in reading and writing
		when nhar = 0.
  03/12/02 jwb  anwrite(): Reinstall 'simple' format output option.
  02/19/03 jwb  cosmetic
  11/16/07 jwb  Change symbol time to timet to avoid system conflict.
  11/05/13 jwb	Remove "How many harmonics to read"
  03/09/15 jwb	anwrite(): Fix nhard and nharw for fa=0 case.
******************************************************************************/
#include "monan.h"
#include "byteorder.h"		                           /* jjjm 01/31/98 */
#define USH_MAX 65535.					    /* jwb 02/24/94 */
#define SH_MAX 32768.					    /* jwb 02/24/94 */

int anread(char*, int);					    /* jwb 03/09/15 */
int anwrite();						    /* jwb 03/09/15 */
static char buffer[80];					    /* jwb 10/12/94 */

/* If anreadcode is 0, # harmonics read is taken from file header */
/* If anreadcode is 1, user is prompted for # harmonics to read */
/* If anreadcode is 2, user is prompted for # harmonics for data 
   to be read into auxiliary arrays */
int anread(char *filnamt, int anreadcode)		    /* jwb 02/21/99 */
{				                      /* rev. jjjm 01/31/98 */
  int i,nbpb,nbph, k,k1,k1inc, narg;			    /* jwb 02/21/99 */
  int unit, type;					    /* jwb 04/23/94 */
  int j, nhar2, fftlen, nitems;                             /* jwb 05/06/96 */
  unsigned short int *tempbufs;			            /* jwb 01/06/94 */
  float *tempbuf, cm, df, dfscale;			    /* jwb 01/06/94 */
  double sum1, sum2, sum3, sum4, sum5;
  double ampscale;					    /* jwb 02/21/99 */
  int byte_reverse;				 /* zheng 10/28/97,rev.jjjm */
/* temporary data variables */
  HEADER *headert;                                     	    /* jwb 02/21/99 */
  int nhart, nhar1t, nptst;                    	    	    /* jwb 02/21/99 */
  float *cmagt, *dfrt, *phaset, *brt, tlt, dtt, fat;   	    /* jwb 02/21/99 */

  /* If the host is little-endian, byte-reverse the data, since SNDAN
   * analysis files are big-endian by specification.                        */
  byte_reverse = byte_order();	               /* zheng 01/13/98, rev. jjjm */

/* file not available */
  if((unit = open(filnamt,0)) == -1)  return(1);  	    /* jwb 02/21/99 */

/* headera and header are allocated in monan.h. headert is a pointer */
  if(anreadcode != 2) headert = &header;		    /* jwb 02/21/99 */
  else headert = &headera;				    /* jwb 02/21/99 */
  rdat(unit,headert);					    /* jwb 02/21/99 */

  P("\nRead file %s \n", filnamt);			    /* jwb 02/21/99 */
  P("Date Recorded: %s\n",headert->date);		    /* jwb 02/21/99 */
  P("Instrument: %s\n", headert->instrument);		    /* jwb 02/21/99 */
  P("Performer: %s\n", headert->performer);		    /* jwb 02/21/99 */
  P("Comments: \n%s\n",headert->comments);		    /* jwb 02/21/99 */
  P("Original sample rate = %6.0f\n", headert->sr);	    /* jwb 02/21/99 */
  P("Pitch: %s; Dynamic: %s\n", 
        headert->pitch, headert->dyn);			    /* jwb 02/21/99 */
  P("Base analysis frequency = %.2f Hz\n", 
	fat = headert->fa);	    			    /* jwb 02/21/99 */
  P("Time between frames = %f sec\n", 
        dtt = headert->dt);		    		    /* jwb 02/21/99 */
  P("Duration of sound = %.4f sec\n",
	tlt = headert->tl);		    		    /* jwb 02/21/99 */
  P("Number of analysis frames = %d\n",
	nptst = headert->npts);			    	    /* jwb 02/21/99 */
  P("Number of harmonics per frame = %d\n", 
	headert->nhar);					    /* jwb 02/21/99 */
  P("Data format: Type %s\n",headert->type);		    /* jwb 02/21/99 */

/* determine nhart (temporary nhar value): */
pt:
#ifdef zilch
  if(anreadcode != 0)					    /* jwb 02/21/99 */
  {
    FP(stderr, 					    	    /* jwb 11/17/95 */
      "How many harmonics to read (default: %d)? ",	    /* jwb 05/06/96 */
		headert->nhar);	    			    /* jwb 02/21/99 */ 
    gets(buffer);					    /* zheng 11/05/97 */
    nitems = sscanf(buffer, "%d", &nhart);		    /* jwb 02/21/99 */ 
    if(nitems == -1) nhart = headert->nhar; 		    /* jwb 02/21/99 */
    else if(nitems == 0) 
    { FP(stderr,"Input not recognized.\n"); goto pt;}       /* jwb 05/06/96 */
    if(nhart > headert->nhar) 				    /* jwb 02/21/99 */
    {
      FP(stderr,"Only %d harmonics available.\n", 
		headert->nhar);				    /* jwb 02/21/99 */
      goto pt;
    }
    else if(nhart < 0)					    /* jwb 02/21/99 */
    {
      FP(stderr,
          "No. of harmonics must be non-negative.\n");      /* jwb 05/06/96 */
      goto pt;						    /* jwb 05/06/96 */
    }
    P("Use nhar = %d\n", nhart);			    /* jwb 02/21/99 */
  }
  else nhart = headert->nhar;				    /* jwb 02/21/99 */

/* note: nhart may be less than headert.nhar */
#endif

  nhart = headert->nhar;				    /* jwb 11/05/13 */
  nhar1t = nhart + 1; 
  nhar2 = max(2, 2*headert->nhar);			    /* jwb 03/28/01 */
  fftlen = headert->fftlen;                            	    /* jwb 02/21/99 */

  /*  allocate space for cmagt, brt, dfrt, and timet arrays: */
  cmagt = (float *) calloc(nptst*nhar1t,sizeof(float));     /* jwb 02/21/99 */
  dfrt = (float *) calloc(nptst*nhar1t,sizeof(float));      /* jwb 02/21/99 */
  brt  = (float *) calloc(nptst+1,sizeof(float));           /* jwb 02/21/99 */
  if(anreadcode != 2)					    /* jwb 02/21/99 */
    timet = (float *) calloc(nptst+1,sizeof(float));	    /* jwb 11/16/07 */

  if(!strcmp(headert->type,"full")) 	  type = 'f';	    /* jwb 02/21/99 */
  if(!strcmp(headert->type,"compressed")) type = 'c';	    /* jwb 02/21/99 */
  if(!strcmp(headert->type,"simple")) 	  type = 's';	    /* jwb 02/21/99 */
  if(!strcmp(headert->type,"compact")) 	  type = 'p';	    /* jwb 02/21/99 */

  if(type=='f')  /* "full" */
  {
    /* read in and save initial phase data    */
     phaset = (float *) calloc(fftlen/2,sizeof(float));     /* jwb 02/21/99 */
     nbph = fftlen*sizeof(float)/2;
     if(read(unit,phaset,nbph) != nbph) 
     {
       P("Cannot read initial phases\n");  
       return(1);
     }
     if (byte_reverse) 					  /* zheng 10/30/97 */
	for(i=0;i<fftlen/2;i++) 
              byteswap4((int*)(phaset+i));
  }

    /* "compressed" or "simple" or "compact" data type */
  else if((type=='c')||(type=='s')||(type=='p'))  	    /* jwb 01/06/94 */
  {
    /* read in and save initial phase data    */
    if(headert->nhar > 0)  /* normal original data */	    /* jwb 08/26/93 */
    {
      phaset = 						    /* jwb 02/21/99 */
          (float *)calloc(headert->nhar+1,sizeof(float));   
      nbph = headert->nhar*sizeof(float);		    /* jwb 02/21/99 */
      if(read(unit,phaset+1,nbph) != nbph) 
      { 
        P("Cannot read initial phases\n");  
        return(1); 
      }
      if(byte_reverse) 					   /* zheng 10/30/97 */
	for(i=1;i<=headert->nhar;i++) 
              byteswap4((int*)(phaset+i));

    }
   /*  headert->nhar eq. 0 for special format for storing only 
       rms & comp. freq. In this case, don't read initial phase data */
  }
  else 
  {
    P("Data type not recognized.  Cannot process."); 
    return(1);
  }

   /* allocate space to read one analysis frame 
                                  (based on original nhar or fftlen) */
  if(type=='f')					    	    /* jwb 03/28/01 */
    nbpb = fftlen*sizeof(float);
  else if((type=='c')||(type=='s'))			    /* jwb 03/28/01 */
    nbpb = nhar2*sizeof(float);				    /* jwb 03/28/01 */
  else /* type 'p' */					    /* jwb 03/28/01 */
    nbpb = nhar2*sizeof(unsigned short);		    /* jwb 03/28/01 */

  tempbuf = (float *) malloc(nbpb);			    /* jwb 03/28/01 */
  if(type=='p')						    /* jwb 03/28/01 */
    tempbufs = (unsigned short *) tempbuf;		    /* jwb 03/28/01 */

  /* scale factors */
  if((type=='f')||(type=='c'))  /*  "full" or "compressed"  */
    ampscale = 1./(fftlen*0.5*.54);
  else if(type=='s')            /* "simple */		    /* jwb 03/28/01 */
    ampscale = 1.; 
  else             /* type=='p', "compact" */
  {
    ampscale = headert->smax;				    /* jwb 02/24/94 */
    if(fat > 0) dfscale = fat/SH_MAX;			    /* jwb 02/21/99 */
    else if(fat == 0) dfscale = 1./8.;			    /* jwb 02/21/99 */
  }
  
  /* read in analysis data for nptst frames and put in cmagt & dfrt arrays */
  for(i=0;i<nptst;i++)  
  {
    if(read(unit,tempbuf,nbpb) != nbpb)
    { P("cannot read frame %d\n",i);
      nptst = i-2;
      break;
    }
    if (byte_reverse)					  /* zheng 11/29/97 */
    {
	if(type=='p') /* compact */
	  for(j=0;j<nhar2;j++) 
             byteswap2((short*)(tempbufs+j));  		    /* jwb 08/19/94 */
	else /* full, compressed, or simple */
	  for(j=0;j<nhar2;j++) 
             byteswap4((int*)(tempbuf+j));     		    /* jwb 08/17/94 */
    }
    if(nhart > 0)	/* normal data */		    /* jwb 03/28/01 */
    {
      if(type=='f') /* for type "full"  */
      {
         k1=4; k1inc = 4;
      } 
      else /* for type "compressed", "simple", or "compact"  */
      {
        k1=0; k1inc = 2;
      } 
      sum1 = sum2 = sum3 = sum4 = sum5 = 0.;
      narg = i*nhar1t;				    	    /* jwb 02/21/99 */
      for(k=1;k<=nhart;k++,k1 += k1inc)		    	    /* jwb 08/26/93 */
      { 
        narg++;
	if(type != 'p') /* not "compact" */		    /* jwb 01/06/94 */
	{
	  cm = ampscale*tempbuf[k1];			    /* jwb 01/06/94 */ 
	  df = tempbuf[k1+1];				    /* jwb 01/06/94 */
	}
	else /* type == 'p', "compact" */
	{
	  cm = ampscale*tempbufs[k1];			    /* jwb 01/06/94 */
	  df = dfscale*k*(tempbufs[k1+1]-SH_MAX);	    /* jwb 02/24/94 */
	}
	cmagt[narg] = cm;				    /* jwb 02/21/99 */
        dfrt[narg]  = df;				    /* jwb 02/21/99 */
/*  if(i==0) P("k=%d, cmag=%f, dfr=%f\n",k,cmagt[narg],dfrt[narg]);  */
 	sum1 += cm*cm;  sum4 += k*abs(cm);		    /* jwb 10/26/94 */ 
	sum5 += abs(cm); 				    /* jwb 10/26/94 */

   /* use NHARMIN harmonics to calculate composite frequency: */
	if(k<=NHARMIN) 
        {
          sum2 += cm;  sum3 += cm*df/k; 		    /* jwb 08/26/93 */
        }    
      } 
      narg = i*nhar1t;				    	    /* jwb 02/21/99 */
      if(sum2 > .0001) dfrt[narg] = sum3/sum2;		    /* jwb 02/21/99 */
      else dfrt[narg] = 0.;				    /* jwb 02/21/99 */

  /* store the rms points:  */
      cmagt[narg] = sqrt(sum1);  	  
  /* store the normalized spectral centroid points: */
      brt[i] = (sum4 + BRTHRESH)/(sum5 + BRTHRESH);
    } /* end if(nhart > 0) */ 
    else /* nhart == 0, special case for storing only rms & comp. dfr */
    {
      if(type != 'p')					    /* jwb 01/06/94 */
      {
        cmagt[i] = tempbuf[0];	 		    	    /* jwb 08/26/93 */
        dfrt[i]  = tempbuf[1];				    /* jwb 08/26/93 */
      }	
      else  /* type == 'p' */   			    /* jwb 01/06/94 */
      {
        cmagt[i] = ampscale*tempbufs[0]; 		    /* jwb 02/24/94 */
        dfrt[i]  = dfscale*(tempbufs[1] - SH_MAX);	    /* jwb 01/06/94 */
      }
    }
   /* store the time points:   */
    if(anreadcode != 2)	timet[i] = i*dtt;             	    /* jwb 11/16/07 */
  }  /* end for(i=0;i<nptst ... loop */
  headert->npts = nptst;				    /* jwb 08/19/94 */
  free(tempbuf);

  /* move data from temporary locations to permanent locations: */
  /* Remember: header and headera are already allocated via monan.h */
  if(anreadcode != 2)
  {
    nhar = nhart; nhar1 = nhar1t; npts = nptst;		    /* jwb 02/21/99 */
    tl = tlt; dt = dtt; fa = fat;			    /* jwb 02/21/99 */
    cmag = cmagt; dfr = dfrt; phase = phaset; br = brt;     /* jwb 02/21/99 */
  }
  else  /* anreadcode == 2 */				    /* jwb 02/21/99 */
  {
    nhara = nhart; nhar1a = nhar1t; nptsa = nptst;	    /* jwb 02/21/99 */
    tla = tlt; dta = dtt; faa = fat;			    /* jwb 02/21/99 */
    cmaga = cmagt; dfra = dfrt; phasea = phaset; bra = brt; /* jwb 02/21/99 */
    auxfile = 1;					    /* jwb 02/21/99 */
  }
  return(0);
} /* end anread() */

/******************************************************************************
 *
 *                   anwrite()
 *
 *         Routine to write mag fr analysis format file (.an file)
 *
 *                 J. Beauchamp    4/28/87
 *
 *****************************************************************************/
int anwrite()						    /* jwb 03/09/15 */
{
  int i,j,k, fd, narg, nbytes, nhard, nharw;		    /* jwb 11/05/95 */
  unsigned short int temp[2]; 				    /* jwb 01/07/94 */
  float dfscale, ampscale, cmax;			    /* jwb 02/24/94 */
  float ftemp[2];					    /* jwb 03/12/02 */
  char outfile[40], resp[10];				    /* jwb 10/03/92 */
  char dfmt;						    /* jwb 03/12/02 */
  int byte_reverse;			         /* zheng 10/28/97,rev.jjjm */

  /* If the host is little-endian we will byte-reverse the data. */
  byte_reverse = byte_order();	               /* zheng 01/13/98, rev. jjjm */

  P("Write (possibly modified) data as a mag dfr analysis file\n\n");
top:
  P("Give name for output file--");
  gets(outfile);					  /* zheng 11/05/97 */
  if (sscanf(outfile,"%s",outfile)==-1) goto top;
  if((fd = open(outfile,0)) != -1)			    /* jwb 10/03/92 */
      { P("This file already exists! Do you wish to overwrite it?(y/n) "); 
        gets(resp);					  /* zheng 11/06/97 */
        if(resp[0] == 'n')				    /* jwb 10/03/92 */
        {close(fd); goto top;}				    /* jwb 10/03/92 */
      }							    /* jwb 10/03/92 */
  fd = creat(outfile,0644);
  header.npts = npts; 
  header.tl = tl;	
pt0:
  if(fa==0) nhard = nhar;				    /* jwb 03/09/15 */
  else 							    /* jwb 03/09/15 */
    nhard = min((int)(header.sr/(2.*fa)),nhar);		    /* jwb 11/05/95 */
  P("How many harmonics to write? (%d max, default): ",
		nhard);					    /* jwb 11/05/95 */
  gets(buffer);						  /* zheng 11/06/97 */
  if(sscanf(buffer, "%d", &nharw) == -1) nharw=nhard;	    /* jwb 12/30/93 */
  if(nharw < 0)						    /* jwb 12/30/93 */
  {
    P("Non-negative number, please!\n");		    /* jwb 12/30/93 */
    goto pt0;						    /* jwb 12/30/93 */
  }
  if(nharw > nhard) 					    /* jwb 11/05/95 */
  { 
    P("Value must be <= %d.\n", nhard); 		    /* jwb 03/09/15 */
    goto pt0;						    /* jwb 08/26/93 */
  }		
  header.nhar = nharw;  
pt1:							    /* jwb 03/12/02 */
  P("compact or simple data format (c (default) or s)? ");  /* jwb 03/12/02 */
  gets(buffer);						    /* jwb 03/12/02 */
  if(sscanf(buffer, "%c", &dfmt) == -1) dfmt = 'c';	    /* jwb 03/12/02 */
  if((dfmt!='c')&&(dfmt!='s'))				    /* jwb 03/12/02 */
  {							    /* jwb 03/12/02 */
    P("Incorrect input! Try again!\n");			    /* jwb 03/12/02 */
    goto pt1;						    /* jwb 03/12/02 */
  }
/* set up parameters for compact file format output */
  if(dfmt=='c')						    /* jwb 03/12/02 */
  {							    /* jwb 03/12/02 */
    header.type = "compact";
    nbytes = 2*sizeof(short);				    /* jwb 01/07/94 */
    if(fa==0) dfscale = 8.;				    /* jwb 01/07/94 */
    else      dfscale = SH_MAX/fa;			    /* jwb 01/07/94 */
/*   find amplitude scale factor */
    cmax = 0.;						    /* jwb 02/24/94 */
    { 
      for(i=0;i<npts;i++)				    /* jwb 02/24/94 */
      {
        if(nharw > 0)					    /* jwb 02/24/94 */
        {
          for(k=1;k<=nharw;k++)				    /* jwb 02/24/94 */
            cmax = max(cmax,cmag[k + i*nhar1]);		    /* jwb 02/24/94 */
        }
        else	/* nhar == 0 */				    /* jwb 02/24/94 */
	  cmax = max(cmax,cmag[i*nhar1]); /* rms case */    /* jwb 02/24/94 */
      }
    }
    ampscale    = USH_MAX/cmax;				    /* jwb 02/24/94 */
    header.smax = 1./ampscale;				    /* jwb 02/24/94 */
    if(nharw > 0) 
      P("Write compact format with %d harmonics\n", nharw); /* jwb 12/30/93 */
    else /* if(nharw == 0) */				    /* jwb 03/12/02 */ 
    {							    /* jwb 03/12/02 */
      P("Write compact format with only rms amplitude ");   /* jwb 03/12/02 */
      P("and average fund. freq. data\n");		    /* jwb 03/12/02 */
    }							    /* jwb 03/12/02 */
  } /* end if(dfmt=='c') */				    /* jwb 03/12/02 */
  else /* dfmt == 's' */				    /* jwb 03/12/02 */
  {							    /* jwb 03/12/02 */
    header.type = "simple";				    /* jwb 03/12/02 */
    nbytes = 2*sizeof(float);				    /* jwb 03/12/02 */
    if(nharw > 0) 					    /* jwb 03/12/02 */
      P("Write simple format with %d harmonics\n", nharw);  /* jwb 03/12/02 */
    else /* if(nharw == 0) */ 				    /* jwb 03/12/02 */	
    {							    /* jwb 03/12/02 */
      P("Write simple format with only rms amplitude ");    /* jwb 03/12/02 */
      P("and average fund. freq. data\n");		    /* jwb 03/12/02 */
    }							    /* jwb 03/12/02 */
  }							    /* jwb 03/12/02 */

/*  write header  */
  wdat(fd,&header);

/* write initial phases, only if nharw is positive */
  if(nharw > 0) 
  {
    if (byte_reverse)
	for(k=1;k<=nharw;k++) byteswap4((int*)(phase+k));/* zheng 10/30/97 */

    write(fd,phase + 1,nharw*sizeof(float)); 	            /* jwb 8/26/93 */
  }
/* write mag, dfr values for nharw harmonics and npts frames */
  for(i=0;i<npts;i++)
  {
    if(nharw > 0)					    /* jwb 08/26/93 */
    {
      for(k=1;k<=nharw;k++)
      { 
        narg = k + i*nhar1;	/* nhar1 is a separator; do not modify */
        if(dfmt=='c')					    /* jwb 03/12/02 */
        {						    /* jwb 03/12/02 */
       	  temp[0] = ampscale*max(cmag[narg],0.);	    /* jwb 05/19/95 */
  	  temp[1] = dfscale*dfr[narg]/k + SH_MAX;	    /* jwb 01/07/94 */
  	  if (byte_reverse)				  /* zheng 10/30/97 */
  	  {						  /* zheng 10/30/97 */
            byteswap2((short*)temp); 			  /* zheng 10/30/97 */
  	    byteswap2((short*)temp+1);  		  /* zheng 10/30/97 */
  	  }
          write(fd,temp,nbytes);
        }						    /* jwb 03/12/02 */
        else /* dfmt=='s' */				    /* jwb 03/12/02 */
        {						    /* jwb 03/12/02 */
          ftemp[0] = cmag[narg];			    /* jwb 03/12/02 */
          ftemp[1] = dfr[narg];				    /* jwb 03/12/02 */
          if (byte_reverse)				    /* jwb 03/12/02 */
          {						    /* jwb 03/12/02 */
            byteswap4((int *)ftemp);			    /* jwb 03/12/02 */
            byteswap4((int *)(ftemp+1));		    /* jwb 03/12/02 */
          }						    /* jwb 03/12/02 */
          write(fd,ftemp,nbytes);			    /* jwb 03/12/02 */
        }						    /* jwb 03/12/02 */

    /*	  if(i==0) P("k=%d, cmag=%f, dfr=%f\n",k,cmag[narg],dfr[narg]);  */
      } /* end for(k... */
    } /* end if(nharw > 0) */

    else /* if(nharw==0): write only rms, dfr.wt.ave to file   jwb 08/26/93 */
    { 
      narg = i*nhar1;
      if(dfmt=='c')					    /* jwb 03/12/02 */
      {							    /* jwb 03/12/02 */
        temp[0] = ampscale*cmag[narg];			    /* jwb 02/24/94 */ 
        temp[1] = dfscale*dfr[narg] + SH_MAX;		    /* jwb 03/28/01 */
        if (byte_reverse)				  /* zheng 10/30/97 */
        {						  /* zheng 10/30/97 */
  	  byteswap2((short*)temp);			  /* zheng 10/30/97 */
  	  byteswap2((short*)(temp+1));  		  /* zheng 10/30/97 */
        }						  /* zheng 10/30/97 */
        write(fd,temp,nbytes);
      }
      else /* if(dfmt=='s') */				    /* jwb 03/12/02 */
      {							    /* jwb 03/12/02 */
        ftemp[0] = cmag[narg];				    /* jwb 03/12/02 */
        ftemp[1] = dfr[narg];				    /* jwb 03/12/02 */
        if (byte_reverse)				    /* jwb 03/12/02 */
        {						    /* jwb 03/12/02 */
            byteswap4((int *)ftemp);			    /* jwb 03/12/02 */
            byteswap4((int *)(ftemp+1));		    /* jwb 03/12/02 */
        }						    /* jwb 03/12/02 */
        write(fd,ftemp,nbytes);				    /* jwb 03/12/02 */
      }							    /* jwb 03/12/02 */
    } /* end else */					    /* jwb 03/12/02 */
  } /* end for(i ... */

  P("Written file is: %s\n",outfile);
  close(fd);
} /* end anwrite() */
