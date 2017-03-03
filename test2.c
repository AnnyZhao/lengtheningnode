/***************************************************************************
                                test.c

	Author: An Zhao  December, 2016

 ***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <fcntl.h>
#include "byteorder.h"
#include "sndhdr.h"
#include "wavhdr.h"
#include "macro.h"
#include "math.h"

#define READ 0
#define	READ	0
#define ERROR	(-1)
#define BUFSIZE	1024				/* Buffer size in floats */
#define BSHSIZE 2048				/*   "     "   in shorts */
#define BSIZE	4096        			/*   "     "   in bytes  */

#define FLOAT 	0
#define SHORT	1
#define SND	0
#define WAV     1
#define MAXTYPES  2
#define P  printf

int filetype;
float *audio, *time,dt;
char *filename;
int byte_reverse = 0;

float* samplesfloat;
float* times;
float* rmsar;
float* windowidx;
float* lengthenedsamples;

char *tail[MAXTYPES] =  {"snd", "wav"};
char graphlabel[200];

int openSoundFile(char*, int*, int* );
int getfiltype(char* );
int plotSamples(int, int, float*, float* );
int plotSpectrum(int, int, float* );
float calcRMS(float*, int, int, int );
int plotseg(float*, float*, int, char*, char*);		    /* jwb 02/04/17 */

int main(int argc, char** argv)
{
  int sampN, sr, fd, i, nwindow;
  float windowdur = 0.03;				    /* jwb 02/04/17 */

  if(argc < 2)
  {
    fprintf(stderr, "Usage: %s <soundfile>", argv[0]);
    P("try again!\n\n");				    /* jwb 02/04/17 */
    exit(1);
  }

//open file and retrieve sr and sampN
    int ret;
    filename = argv[1];
    P("input filename is %s \n", filename);
    fd = openSoundFile(filename, &sr, &sampN);
    if(fd == -1)
    {
      fprintf(stderr, "cannot open file %s\n", filename);
      exit(1);
    }
    P("file sample rate = %d, No. samples = %d\n", sr, sampN); /*jwb 2/4/17 */
    float origDur = (float)sampN/sr;			    /* jwb 02/04/17 */
    P("orignal sound duration is %.3f sec\n", origDur);

    P("For rms envelope calculation:\n");		    /* jwb 02/04/17 */
    int window = (float)sr*windowdur;			    /* jwb 02/04/17 */
    P("window length is %.3f sec, %d samples, ",	    /* jwb 02/04/17 */
	windowdur, window);				    /* jwb 02/04/17 */
    nwindow = sampN/window + 1;
    P("No. windows is %d\n", nwindow);
    rmsar = (float*)calloc(nwindow, sizeof(float));
    windowidx = (float*)calloc(nwindow, sizeof(float));

//  read all short int samples from sound file
    short int* samples = (short int*)calloc(sampN, sizeof(short int));
    ret = read(fd, samples, sampN*sizeof(short int));
//  P("ret %d \n",ret);
    if(ret == -1)
    {
      fprintf(stderr, "cannot read samples from %s\n", filename);
      exit(1);
    }
//  transfer short int samples to float array
    samplesfloat = (float*)calloc(sampN, sizeof(float));
    times = (float*)calloc(sampN, sizeof(float));
    for(i = 0; i < sampN; i++)
    {
      samplesfloat[i] = samples[i];
      times[i] = i/(float)sr;
      // if (i < 10)
      // {
      //     P("samples %d \n", samples[i]);
      //     P("samples %f \n", samplesfloat[i]);
      //     P("times %f \n", times[i]);
      // }
    }
//  plotSamples(sr, sampN, times, samplesfloat);

//  calculation of rms array
    i = 0;
    for(i = 0; i < nwindow; i++)
    {
      rmsar[i] = calcRMS(samplesfloat, i*window, window, sampN);
      windowidx[i] = i*window/(float)sr;
      // if (i < 72)
      // {
      //     P("rmsar[%d] is %f, windowidx[%d] is %f \n",
      //     i, rmsar[i], i, windowidx[i]);
      //
      // }
    }
//  plotSamples(sr, nwindow, windowidx, rmsar);
//  convert rms array to dB units and compute average dB value
    float* dbarr = (float*)calloc(nwindow, sizeof(float));
    i = 0;
    float sumdb = 0;
    for (i = 0; i < nwindow; i++)
    {
      dbarr[i] = 20.*log10f(rmsar[i]);
      // P("db[%d] is %.1f\n", i, dbarr[i]);
      sumdb += dbarr[i];
    }
    float avgdb = sumdb/nwindow;
    P("db average = %.1f\n", avgdb);
//  plotseg(windowidx,dbarr,nwindow, "Time(sec)","dB");

//  compute attack and decay frames
//  first pass calculation:
    float firstdiff = 0, seconddiff = 0;
    int attackf, decayf;
    for (i = 0; i < nwindow; i++)
    {
      firstdiff = seconddiff;
      seconddiff = dbarr[i] - avgdb;
      if(seconddiff > 0 && firstdiff < 0)
      {
        attackf = i;
        // P("attack frame is %d\n", attackf);
      }
      if(seconddiff < 0 && firstdiff > 0)
      {
        decayf = i;
        // P("decay frame is %d\n", decayf);
      }
    }
//  second pass calculation:
//  calculate average dB between first time attack and decay
    sumdb = 0;
    for(i = attackf; i < decayf; i++)
      sumdb += dbarr[i];
    avgdb = sumdb/(decayf - attackf);
    P("adjusted db average = %.1f\n", avgdb);

//  compute attack-end and decay-start frames
    firstdiff = 0; seconddiff = 0;
    for(i = 0; i < nwindow; i++)
    {
      firstdiff = seconddiff;
      seconddiff = dbarr[i] - avgdb;
      if(seconddiff > 0 && firstdiff < 0)
      {
        attackf = i;
        // P ("attack frame is %d\n", attackf);
      }
      if(seconddiff < 0 && firstdiff > 0)
      {
        decayf = i;
        // P("decay frame is %d\n", decayf);
      }
    }
//  calculate attack-end and decay-start times
    float attackt, decayt, mduration;
    attackt = (attackf-1+(avgdb-dbarr[attackf-1])/(dbarr[attackf]
      -dbarr[attackf-1]))*windowdur;
    decayt = (decayf-1+(dbarr[decayf-1]-avgdb)/(dbarr[decayf-1]
      -dbarr[decayf]))*windowdur;
//  time between attack-end and decay-start
    mduration = decayt - attackt;
    P("Based on the adjusted dB average:\n");
    P("attack-end time = %.4f, decay-start time = %.4f, timediff = %.4f\n",
       attackt, decayt, mduration);

//  elongate the note
    float totalt, extendt, ratio;			    /* jwb 02/03/17 */
    while(1)						    /* jwb 02/03/17 */
    {							    /* jwb 02/03/17 */
      P("Give new time duration: ");			    /* jwb 02/03/17 */
      scanf("%f", &totalt);
      extendt = totalt - origDur;
      ratio = extendt/mduration;
      if(ratio <= 2.0) break;				    /* jwb 02/03/17 */
      P("(time extension) > 2.*(time difference)\n");	    /* jwb 02/03/17 */
      P("This condition is not currently allowed.\n");	    /* jwb 02/03/17 */
      P("Try again.\n");				    /* jwb 02/03/17 */
    }							    /* jwb 02/03/17 */
    int numtotalt = totalt*sr;
    P("total samples in elongated file = %d\n", numtotalt);
    float reverset = decayt - 0.5*extendt;		    /* jwb 02/03/17 */
    P("Time added to extend duration is %.4f\n",extendt);   /* jwb 02/03/17 */
    P("Time is added by looping between decay-start\n");    /* jwb 02/03/17 */
    P("and decay-start minus 0.5*time-extension\n");	    /* jwb 02/03/17 */
    P("Based the original file, turn-around times will ");  /* jwb 02/03/17 */
    P("be at t = %.4f and t = %.4f\n", reverset, decayt);   /* jwb 02/03/17 */
    P("Based on the elongated file these points will ");    /* jwb 02/03/17 */
    P("occur at t = %.4f and t = %.4f\n\n", 		    /* jwb 02/03/17 */
        decayt, decayt + 0.5*extendt);			    /* jwb 02/03/17 */
    lengthenedsamples = (float*)calloc(numtotalt,sizeof(float));
    float* newtimes = (float*)calloc(numtotalt,sizeof(float));
    int attacksample = attackt*sr;
    int decaysample = decayt*sr;
    int extendsample = extendt*sr;
    for(i = 0; i < decaysample; i++)
      lengthenedsamples[i] = samplesfloat[i];
    if(ratio <= 2.)
    {
      int j;
      for(j = decaysample; j > decaysample - extendsample/2; j--)
        lengthenedsamples[i++] = samplesfloat[j];
      for(j = decaysample - extendsample/2; j < decaysample; j++)
        lengthenedsamples[i++] = samplesfloat[j];
      for(j = decaysample; j < sampN; j++)
        lengthenedsamples[i++] = samplesfloat[j];
      for(i = 0; i < numtotalt; i++)
        newtimes[i] = (float)i/sr;
    }
    short int* lengthenedsamplesint =
      (short int*)calloc(numtotalt,sizeof(short int));
//  P("total samples of elongated file = %d\n", numtotalt);
    for(i = 0; i < numtotalt; i++)
       lengthenedsamplesint[i] = lengthenedsamples[i];
//  plotSamples(sr, numtotalt, newtimes, lengthenedsamples);
    char* outfile = "elongated.wav";			    /* jwb 02/03/17 */
    int fdnew = creat(outfile, 0644);			    /* jwb 02/03/17 */
//  P ("new fd is %d\n", fdnew);
    P("Write elongated file %s\n\n", outfile);		    /* jwb 02/03/17 */
    writeWavHdr(fdnew, sr, 1, numtotalt, 16);
    write(fdnew, lengthenedsamplesint, 2*numtotalt);
}  /* end main() */

float calcRMS(float* samplesfloat, int begintime, int duration, int sampN)
{
    float rms, curr;
    int i;
    float sum = 0;
    for(i = 0; i < duration-1; i++)
    {
      if(sampN > begintime + i)
      {
        curr = sq(samplesfloat[begintime + i]);
        sum += curr;
      }
      else
      {
        duration = i;
        break;
      }
    }
    rms = sqrt(sum/duration);
    return rms;
}

int getfiltype(char*name)
{
    int i, len;
    len = strlen(name);
    for(i=0;i<MAXTYPES;i++)
    {
        if(!strcmp(&name[len-3],tail[i]))  return(i);
    }
    return(MAXTYPES);
}

int openSoundFile(char* filename, int* samplerate, int* sampN)
{
    int fd, sr, nsamps, nchans, samptype,sampsize;
//    P("in function openSoundFile :\n");
//    P("line2 %s \n",filename);
    fd = open(filename, READ);
    if (fd == (-1))
    {
        fprintf(stderr, "cannot open file %s\n", filename);
        exit(1);
    }
    filetype = getfiltype(filename);
//  P("filetype = %d\n", filetype);


    if(filetype == SND)
	{
	  readSndHdr(fd, &sr, &nchans, &nsamps, &samptype);
      *samplerate = sr;
      *sampN = nsamps;
//    P ("sr = %d, nchans = %d, sampN = %d, samptype = %d\n",sr,nchans,*sampN,samptype);
	  if (samptype == SND_FORMAT_LINEAR_16)
	    filetype = SHORT;
	  else
	    {
	      fprintf( stderr, "Unsupported sample type in '.snd' file.\n");
	      exit(1);
	    }
        byte_reverse = (byte_order() == little_endian);
    }
    else if (filetype == WAV)
    {
        readWavHdr(fd, &sr, &nchans, &nsamps, &sampsize);
        /* WAVE files in PCM format only support integers, not floats. */
        *samplerate = sr;
        *sampN = nsamps;
//      P ("sr = %d, nchans = %d, sampN = %d, sampsize = %d\n",sr,nchans,*sampN,sampsize);
        if (sampsize == 16)
            filetype = SHORT;
        else
        {
            fprintf(stderr,"WAVE file samples in %d-bit format.\n", sampsize);
            exit(1);
        }
        /* WAVE data is little-endian, so if the host is big-endian,
        byte-reverse the data. */
        byte_reverse = (byte_order() != little_endian);
    }
    return fd;
}

int plotSamples(int sr, int sampN, float* times, float* samplesfloat)
{
    sprintf(graphlabel, "signal amplitude for file ");
    int i1, i2;
    float t1, t2;
    P("plot some samples of the input sound file\n");	    /* jwb 02/02/17 */
    P("Give time range (t1 t2): ");			    /* jwb 02/02/17 */
    scanf("%f%f%*c", &t1, &t2);				    /* jwb 02/02/17 */
    i1 = t1 * sr;
    P("i1 =  %d\n", i1);
    i2 = t2 * sr;
    P("i2 = %d\n",i2);
    int length = i2 - i1;
    P("length is %d\n", length);
    plotseg(times+i1,samplesfloat+i1,length, "TIME(SEC)","AMPLITUDE");
    return 0;
}

int plabel(double xpos, double ypos)
{
    g_move_abs(xpos,ypos);
    g_text(graphlabel);
    return 0;
}
