/***************************************************************************
        addsynextend.c (aka extend.c)

	Author: An Zhao  March, 2017

	Latest edit:  11/25/17

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
#include <math.h>
#include "header.h"

#define PI 3.1415926534898
#define ERROR	(-1)
#define TABSIZ	5000
#define FLOAT 	0
#define SHORT	1
#define SND	0
#define WAVE     1
#define MAXFTYPES  2
#define P  printf

const char *tail[MAXFTYPES] =  {"snd", "wav"};

/* function declarations */
char *gettime();
int getfiltype();
void getout();
int anread(char*);
void cutSilence(float** cmag, float** dfr);
void calcRMS(float* cmag);
void extendsyn(float** cmag, float** dfr, int nhar1, float length, float extension, float frameDuration);
void reddursyn(float** Cmag, float** Dfr, int nhar1, float length, float origDur, float frameDuration);
void addsyn(char* filnam, int byte_reverse);
float blend(float x);
void findattackdecay(int *attackf, int *decayf, float frameDuration);
void timescalemin(float** Cmag, float** Dfr, int nhar1, float length, float mintime, float origDur, float frameDuration);
int plabel(double, double);
// void plotseg(float*, float*, int, char*, char*);
int plotseg(float*, float*, int, char*, char*);

/* global variables: */
HEADER header;
int nhar, nhar1, npts;
float *cmag, *dfr, *phase, dt, tl, smax, fs, fa, scalefac;
char *anfil, *synfil;					    /* jwb 11/25/17 */
char graphlabel[400];					    /* jwb 11/25/17 */
int axskip,doplot,research;
int barplot;
int color_on, color, linestyle;

/* main function */
int main(int argc, char** argv )
{
  int byte_reverse = 0;
  if(argc < 3)
  {
    fprintf(stderr, "Usage: %s <analysis file> <output file>\n", argv[0]);
    P("try again!\n\n");				    /* jwb 02/04/17 */
    exit(1);
  }
  if( *argv[1] == '-' ) 	/* check for option flag */
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
  //cutSilence(&cmag, &dfr);
  P("number of harmonics = %d, number of frames = %d\n", 
     nhar, npts);

  float frameDuration = dt;
  float origDur = tl;
  P("original sound duration is %.3f sec:\n", origDur);
  // increment for traversing the cmagi and dfri data

  calcRMS(cmag);
  int attackf, decayf;
  findattackdecay(&attackf, &decayf, frameDuration);

  float attackt = attackf * dt;
  float decayt = decayf * dt;
  P("dt is %f\n", dt);
  P("attackt is %f and decayt is %f\n", attackt, decayt);

  float totalt, extendt;
  float ratio;
  float minimumT = attackt + (origDur - decayt);
  P("minimum duration is %.3f\n", minimumT);
  while(1)
  {
    P("\nGive new time duration: ");
    scanf("%f", &totalt);
    if (totalt > minimumT) break;
    P("time is too short to reduce\n");
    P("Try again.\n");
  }
  //extend
  if(totalt > origDur)
  {
    // decayf = decayf - (decayf - attackf) * 0.05;
    // attackf = attackf + (decayf - attackf) * 0.05;
    float mduration = (decayf - attackf) * dt;
    P("mduration is %f\n", mduration);
    fflush(stdout);
    extendt = totalt - origDur;
    P("extendt is %f\n", extendt);
    fflush(stdout);
    ratio = extendt / mduration;
    P("ratio is %f\n", ratio);
    fflush(stdout);
    extendsyn(&cmag, &dfr, nhar1, totalt, extendt, frameDuration);
  }
  //shorten
  else if(totalt < origDur)
  {
    if(totalt < minimumT) // timescale compressed attack+decay
      timescalemin(&cmag,&dfr,nhar1,totalt,minimumT,origDur,frameDuration);
    else
      reddursyn(&cmag,&dfr,nhar1,totalt,origDur,frameDuration);
  }
  /*  resynthesize tone at sample rate fs */
  P("\nBegin synthesis\n");
  addsyn(synfil, byte_reverse);
  P("\nSynthesis of file %s complete\n\n", synfil);

} /*  end main() */

//functions
float square(float x) 
{
  return x*x;
}

void calcRMS(float* cmag)
{
  int i, k;
  float sum = 0.0;
  for(i = 0 ; i < npts; i++)
  {
    for(k = 1; k < nhar1; k++)
    {
      sum += square(cmag[k + i * nhar1]);
    }
    sum = sqrt(sum);
    cmag[i*nhar1] = sum;
  }
}

// Butterworth filter from http://baumdevblog.blogspot.com/2010/11/butterworth-lowpass-filter-coefficients.html

void getLPCoefficientsButterworth2Pole(const int samplerate, 
        const double cutoff, double* const ax, double* const by)
{
  double sqrt2 = 1.4142135623730950488;

  // Find cutoff frequency in [0..PI]
  double QcRaw  = (2 * PI * cutoff) / samplerate; 
  double QcWarp = tan(QcRaw); // Warp cutoff frequency

  double gain = 1 / (1+sqrt2/QcWarp + 2/(QcWarp*QcWarp));
  by[2] = (1 - sqrt2/QcWarp + 2/(QcWarp*QcWarp)) * gain;
  by[1] = (2 - 2 * 2/(QcWarp*QcWarp)) * gain;
  by[0] = 1;
  ax[0] = 1 * gain;
  ax[1] = 2 * gain;
  ax[2] = 1 * gain;
}

void ButterworthFilter(float* samples, float* samplespassed, int count, int sampleRate, double cutoff)
{
  double xv[3] = {samples[0]};
  double yv[3] = {samples[0]};
  double ax[3] = {0};
  double by[3] = {0};

  // cutoff = 20
  getLPCoefficientsButterworth2Pole(sampleRate, cutoff, ax, by);

  for (int i = 0; i < count; i++)
  {
    xv[2] = xv[1]; xv[1] = xv[0];
    xv[0] = samples[i];
    yv[2] = yv[1]; yv[1] = yv[0];

    yv[0] = (ax[0] * xv[0] + ax[1] * xv[1] + ax[2] * xv[2]
              - by[1] * yv[0] - by[2] * yv[1]);
    samplespassed[i] = yv[0];
  }
}  /* end ButterworthFilter() */
void findattackdecay(int *attackf, int *decayf, float frameDuration)
{
  //calculate the attack and decay
  //k = 0: rms amplitude
  //original no. of frames
  int i = 0;
  int debugattack = 10000;
  float frameRate = 1.0 / frameDuration;
  float* dbarr = (float*)calloc(npts+1, sizeof(float));
  float* dbarrlowpassed = (float*)calloc(npts+1, sizeof(float));
  P("npts in findattackdecay = %d\n", npts);
  float sumdb = 0;
  int npts_nonzero = 0;
  for(i = 0; i < npts; i++)
  {
    // if (i < 10)
    // {
    //     P("the first 0 harmonic in frame %d, %f\n",i, cmag[i * nhar1]);
    // }
    dbarr[i] = 20.*log10f(cmag[i * nhar1]);
    if (dbarr[i] <= 0)
    {
      dbarr[i] = 0.0001;
    }
    sumdb += dbarr[i];
    //if (dbarr[i] > 17.0) npts_nonzero += 1;
  }
  float avgdb = sumdb / npts;
  P("db average = %.1f\n",avgdb);

  float cutoff = 20;
  ButterworthFilter(dbarr, dbarrlowpassed, npts, frameRate, cutoff);

  //first pass calculation
  float firstdiff = 0, seconddiff = 0;

  P("npts=%d\n",npts);
  for (i = 0; i < npts; i++)
  {
    firstdiff = seconddiff;
    seconddiff = dbarr[i] - avgdb;
    if(seconddiff > 0 && firstdiff < 0)
    {
      if(debugattack > i)
      {
        *attackf = i;
        debugattack = i;
        P("attack frame is %d\n", *attackf);
      }
    }
    if(seconddiff < 0 && firstdiff > 0)
    {
      *decayf = i;
      P("decay frame is %d\n", *decayf);
    }
  }
} /* end findattackdecay() */

float blend(float x)
{
  return x;
}

float interpolation(float start, float end, float percentage)
{
  return start + (end - start) * percentage;
}

void timescalemin(float** Cmag, float** Dfr, int nhar1, float length, 
      float mintime, float origDur, float frameDuration)
{
  if (length > mintime) return; // not suitable for timescale
  //modify to minimumT
  reddursyn(Cmag,Dfr,nhar1,mintime*1.01,origDur,frameDuration);

  float* cmagold = *Cmag;
  float* dfrold = *Dfr;
  int nptsnew = length/dt;
  *Cmag = (float*) calloc(nptsnew * nhar1, sizeof(float));
  *Dfr = (float*) calloc(nptsnew * nhar1, sizeof(float));
  P("the size of cmag and dfr is %d\n", nptsnew*nhar1);// create memory space
  float coeff;
  coeff = mintime/length; // calculate number of samples to jump

  //interpolate between each sample
  int i, k;
  for(i = 0; i < nptsnew; i++)
  {
    float ptOld = i * coeff;
    int pt1 = floor(ptOld);
    int pt2 = ceil(ptOld);
    float percentage = ptOld - pt1;
    if (pt2 >= npts) pt2 = npts - 1;
    for(k = 1; k < nhar1; k++)
    {
      (*Cmag)[k + i * nhar1] = interpolation(cmagold[k + pt1 * nhar1],
          cmagold[k + pt2 * nhar1], percentage);
      (*Dfr)[k + i * nhar1] = interpolation(dfrold[k + pt1 * nhar1],
          dfrold[k + pt2 * nhar1], percentage);
    }
  }
  tl = length;
  npts = nptsnew;
}

void reddursyn(float** Cmag, float** Dfr, int nhar1, float length, 
       float origDur, float frameDuration)
{
  float x,w;
  int attackf, decayf;
  float frameRate = 1.0 / frameDuration;
  findattackdecay(&attackf, &decayf, frameDuration);
  decayf = decayf - (decayf - attackf) * 0.05;
  attackf = attackf + (decayf - attackf) * 0.05;
  int nptsnew = length / dt;
  P("nptsnew: %d\n", nptsnew);
  int overlapF = nptsnew - attackf - (npts - decayf);
  P("overlapF: %d\n", overlapF);
  int shiftF = npts - nptsnew;
  P("shiftF: %d\n",shiftF);
  x = 1./overlapF; //function factor

  float* cmagold = *Cmag;
  float* dfrold = *Dfr;

  *Cmag = (float*)calloc(nptsnew * nhar1, sizeof(float));
  *Dfr = (float*)calloc(nptsnew * nhar1, sizeof(float));
  P ("the new size of cmag and dfr is %d\n", nptsnew * nhar1);

  int i, j, k, narg1, narg2;
  // beginning to attackf
  for(i = 0; i < attackf; i++)
  {
    for (k = 1; k < nhar1; k++)
    {
        (*Cmag)[k + i * nhar1] = cmagold[k + i * nhar1];
        (*Dfr)[k + i * nhar1] = dfrold[k + i * nhar1];
    }
  }
  //start to blend
  for(j = attackf; j < attackf + overlapF + 1; j++)
  {
    w = blend(x*(j-attackf));
    int num = j - attackf;
    if (num < 10)
    {
        P("for frame %d the coeffienct is %f\n", num, w);
        P("crossfade frame %d with frame %d\n", j, j+shiftF);
    }
    int narg1, narg2;
    for (k = 1; k < nhar1; k++)
    {
        narg1 = k + j * nhar1;
        narg2 = k + (j + shiftF) * nhar1;
        (*Cmag)[k + i * nhar1] = (1-w) * cmagold[narg1] + w * cmagold[narg2];
        (*Dfr)[k + i * nhar1] = (1-w) * dfrold[narg1] + w * dfrold[narg2];
    }
    i++;
  }
  //decay to end
  for(j = decayf; j < npts; j++)
  {
    for(k = 1; k < nhar1; k++)
    {
      (*Cmag)[k + i * nhar1] = cmagold[k + j * nhar1];
      (*Dfr)[k + i * nhar1] = dfrold[k + j * nhar1];
    }
    i++;
    if (i >= nptsnew) break;
  }
  //update parameters for addsyn
  npts = nptsnew;
  tl = length;

} /* end reddursyn */

void extendsyn(float** cmag, float** dfr, int nhar1, float length, 
  float extension, float frameDuration)
{
  int i = 0;
  int k = 0;
  int attackf, decayf;
  int frameRate = 1.0 / frameDuration;
  findattackdecay(&attackf, &decayf, frameDuration);
  // decayf = decayf - (decayf - attackf) * 0.1;
  // attackf = attackf + (decayf - attackf) * 0.1;

  float mduration = (decayf - attackf) * dt;
  float ratio = extension / mduration;

  int nptsnew = length / dt;
  P("nptsnew: %d\n", nptsnew);
  int extendframes = nptsnew - npts;

  float* cmagold = *cmag;
  float* dfrold = *dfr;
  // calculate and interpolate amplitude variation

  // allocate data arrays
  float** amporiginal = (float**) calloc(nhar1, sizeof(float*));
  float** amplowpassed = (float**) calloc(nhar1, sizeof(float*));
  //new data
  float** amplowpassednew = (float**) calloc(nhar1, sizeof(float*));
  float** microvariations = (float**) calloc(nhar1, sizeof(float*));
  for(k = 1; k < nhar1; k++)
  {
  amporiginal[k] = (float*) calloc(npts, sizeof(float));
  amplowpassed[k] = (float*) calloc(npts, sizeof(float));
  microvariations[k] = (float*) calloc(npts, sizeof(float));
  amplowpassednew[k] = (float*) calloc(nptsnew, sizeof(float));
  }
  for(k = 1; k < nhar1; k++)				    /* jwb 11/27/17 */
  {							    /* jwb 11/27/17 */
  for(i = 0; i < npts; i++)				    /* jwb 11/27/17 */
    amporiginal[k][i] = (*cmag)[k + i * nhar1];	    /* jwb 11/27/17 */
  }							    /* jwb 11/27/17 */
  int kplot;						    /* jwb 11/17/17 */
  P("\nGive harmonic number for amplitude-vs-time plots: ");
  scanf("%d%*c", &kplot);				    /* jwb 11/17/17 */
  float *Time = (float*)calloc(npts,sizeof(float));	    /* jwb 11/17/17 */
  float *Haramp = (float*)calloc(npts,sizeof(float));     /* jwb 11/17/17 */
  for(i=0; i<npts; i++)				    /* jwb 11/17/17 */
  {							    /* jwb 11/17/17 */
  Haramp[i] = amporiginal[kplot][i];		    /* jwb 11/17/17 */
  Time[i] = i*dt;					    /* jwb 11/17/17 */
  }							    /* jwb 11/17/17 */
  P("\nHit Return to plot the original harmonic %d amplitude curve: ",
  kplot);						    /* jwb 11/28/17 */
  getchar();						    /* jwb 11/17/17 */
  sprintf(graphlabel,"file: %s, original amplitude",	    /* jwb 11/25/17 */
  anfil);						    /* jwb 11/25/17 */
  char svlabel[80];					    /* jwb 11/17/17 */
  sprintf(svlabel,"HARMONIC %d AMPL", kplot);		    /* jwb 11/29/17 */
  plotseg(Time,Haramp,npts,"TIME (SEC)",svlabel);	    /* jwb 11/17/17 */

  P("\nApply lowpass filter to smooth original harmonic amplitudes\n");
  P("except for attack and decay epochs\n");
  float cutoff;					    /* jwb 11/10/17 */
  P("\nGive lowpass filter cutoff frequency: ");	    /* jwb 11/10/17 */
  scanf("%f%*c", &cutoff);				    /* jwb 11/10/17 */
  for (k = 1; k < nhar1; k++)
  {
  ButterworthFilter(amporiginal[k]+attackf,amplowpassed[k]+attackf, 
  npts-attackf-(npts-decayf),frameRate,cutoff);
  for(i = 0; i < attackf; i++)
    amplowpassed[k][i] = amporiginal[k][i];
  for(i = decayf; i < npts; i++)
    amplowpassed[k][i] = amporiginal[k][i];
  }
  //  compute frames shift for low-passed amplitude-vs-time
  P("\nShift harmonic amplitudes to compensate for delay\n");
  int shiftamount=1./(sqrt(2.)*PI*cutoff*frameDuration)+0.5; /*jwb 11/17/17*/
  P("No. frames to compensate for LPF delay: %d,\n", shiftamount);
  P("corresponding to t_del = %.3f\n", shiftamount*dt);    /* jwb 11/28/17 */

  int shift;
  for(k = 1; k < nhar1; k++)
  {
  for (i = shiftamount + attackf; i < npts; i++)
  {
    shift = i - shiftamount;
    amplowpassed[k][shift] = amplowpassed[k][i];
  }
  for(shift = npts - shiftamount; shift < npts; shift++)
    amplowpassed[k][shift] = amplowpassed[k][npts - shiftamount - 1];
  for(i = 0; i < attackf; i++)
    amplowpassed[k][i] = amporiginal[k][i];
  for(i = decayf; i < npts; i++)
    amplowpassed[k][i] = amporiginal[k][i];
  }
// printCSV("amplowpassed1.csv", "Low-pass (1st partial)", amplowpassed[1], npts, frameDuration);
// printCSV("amplowpassed2.csv", "Low-pass (2st partial)", amplowpassed[2], npts, frameDuration);

  float* Harampagain = (float*)calloc(npts,sizeof(float));
  float* Timeagain = (float*)calloc(npts,sizeof(float));
  for(i=0; i<npts; i++)
  {
    Harampagain[i] = amplowpassed[kplot][i];
    Timeagain[i] = i*dt;
  }
  P("\nHit Return to plot the smoothed harmonic %d amplitude: ",
  kplot);						    /* jwb 11/28/17 */
  getchar();						    /* jwb 11/17/17 */
  char *exlab1, *exlab2; 				    /* jwb 11/25/17 */
  exlab1 = "smoothed amplitude";			    /* jwb 11/28/17 */
  exlab2 = "lowpass cutoff=";			    	    /* jwb 11/25/17 */
  sprintf(graphlabel,"file: %s, %s, %s%.1f Hz",	    /* jwb 11/28/17 */
  anfil,exlab1,exlab2,cutoff);			    /* jwb 11/28/17 */
  plotseg(Timeagain,Harampagain,npts,"TIME (SEC)",svlabel);

  P("\nSubtract the smoothed from the original harmonic amplitudes\n"); 
  P("to get the microvariations\n");			    /* jwb 11/28/17 */
  for(k = 1; k < nhar1; k++)
  {
    for(i = 0 ; i < npts; i++)
    {
      microvariations[k][i] = amporiginal[k][i] - amplowpassed[k][i];
      if(amporiginal[k][i] == 0) microvariations[k][i] = 0;
      if(amplowpassed[k][i] == 0) microvariations[k][i] = 1;
    }
  }
  float* zigzags = (float*)calloc(npts,sizeof(float));
  float* Timezig = (float*)calloc(npts,sizeof(float));
  for(i=0; i<npts; i++)
  {
    zigzags[i] = microvariations[kplot][i];
    Timezig[i] = i*dt;
  }
  P("\nHit Return to plot the harmonic %d microvariations: ",kplot);
  getchar();
  exlab1 = "microvariations";				    /* jwb 11/28/17 */
  exlab2 = "sans attack and decay";		 	    /* jwb 11/25/17 */
  sprintf(graphlabel,"file: %s, %s %s",		    /* jwb 11/28/17 */
  anfil,exlab1,exlab2);				    /* jwb 11/28/17 */
  plotseg(Timezig,zigzags,npts,"TIME (SEC)",svlabel);

  // crossfade into unshifted decay
  for(k = 1; k < nhar1; k++)
  {
    int lowerReach = min(decayf, 100);
    int length = lowerReach;
    for(i = decayf - lowerReach; i < decayf; i++) 
    {
      float percentShifted = 1 - ((float)(i - decayf + lowerReach))/length;
      float percentUnshifted = 1 - percentShifted;
      amplowpassed[k][i] = amplowpassed[k][i]* 
        percentShifted + amporiginal[k][i]*percentUnshifted;
    }
    for(i = decayf; i < npts; i++)
      amplowpassed[k][i] = amporiginal[k][i];
  }
  //time scale amplowpassed
  P("\nElongate the smoothed harmonic amplitudes\n\n");
  for(k = 1; k < nhar1; k++)
  {
    for(i = 0; i < nptsnew; i++)
    {
      int intoSustain = i - attackf;
      if(intoSustain <= 0)
        amplowpassednew[k][i] = amporiginal[k][i];
      else if((nptsnew - i) <= (npts - decayf))
        amplowpassednew[k][i] = amporiginal[k][(npts - (nptsnew - i))];
      else
      {
        float fullratio = ratio + 1;
        int a = floor((((float)i) - attackf)/fullratio) + attackf;
        int b = ceil((((float) i) - attackf)/fullratio) + attackf;
        float percentage = ((((float)i) - attackf)/fullratio) + attackf - a;
        if(a >= npts) a = npts - 1;
        if (b >= npts) b = npts - 1;
        amplowpassednew[k][i] 
        = amplowpassed[k][b]*percentage + amplowpassed[k][a]*(1 - percentage);
      }
    }
  }
  float* Harampnew = (float*)calloc(nptsnew,sizeof(float));
  float* Timenew = (float*)calloc(nptsnew,sizeof(float));
  for(i=0; i<nptsnew; i++)
  {
    Harampnew[i] = amplowpassednew[kplot][i];
    Timenew[i] = i*dt;
  }
  P("Hit Return to plot the new smoothed harmonic %d amplitude: ", kplot);
  getchar();
  exlab1 = "smoothed & elongated amplitude";		    /* jwb 11/28/17 */
  sprintf(graphlabel,"file: %s, %s", anfil,exlab1);	    /* jwb 11/28/17 */
  plotseg(Timenew,Harampnew,nptsnew,"TIME (SEC)",svlabel);

  // printCSV("amporiginal1.csv", "Original amplitude (1st partial)", amporiginal[1], npts, frameDuration);
  // printCSV("amporiginal2.csv", "Original amplitude (2st partial)", amporiginal[2], npts, frameDuration);
  //
  // printCSV("amplowpassednew1.csv", "Low-pass and stretch (1st partial)", amplowpassednew[1], nptsnew, frameDuration);
  // printCSV("amplowpassednew2.csv", "Low-pass and stretch (2st partial)", amplowpassednew[2], nptsnew, frameDuration);
  //
  // printCSV("amplowpassed1cross.csv", "Low-pass and crossfade (1st partial)", amplowpassed[1], npts, frameDuration);
  // printCSV("amplowpassed2cross.csv", "Low-pass and crossfade (2st partial)", amplowpassed[2], npts, frameDuration);

  // loop microvariations
  P("\nExtend the microvariations using the zig-zag method\n");
  *cmag = (float*)calloc(nptsnew * nhar1, sizeof(float));
  *dfr = (float*)calloc(nptsnew * nhar1, sizeof(float));

  int samplePointer = 0;
  int fullloop = ratio/2;
  float amplitudeScale = 0;
  P("\nnumber of full loops: %d\n", fullloop);
  // first beginning to decayf
  int j;
  for(samplePointer = 0; samplePointer < decayf; samplePointer++)
  {
    for(k = 1; k < nhar1; k++)
    {
      if(samplePointer <= attackf)
        (*cmag)[k + samplePointer * nhar1] = cmagold[k + samplePointer * nhar1];
      else
      {
        // amplitudeScale = 0;
        // if(amplowpassed[k][samplePointer] == 0)
        //     amplitudeScale = 1;
        // else
        // {
        //   amplitudeScale = amplowpassednew[k][samplePointer]/ 
        //     amplowpassed[k][samplePointer];
        // }
        //(*cmag)[k + samplePointer*nhar1] = cmagold[k + samplePointer*nhar1]* 
        //(amplowpassednew[k][samplePointer] / amplowpassed[k][samplePointer]);

        (*cmag)[k + samplePointer*nhar1] = amplowpassednew[k][samplePointer] 
            + amporiginal[k][samplePointer] - amplowpassed[k][samplePointer];
        (*dfr)[k + samplePointer * nhar1] = dfrold[k + samplePointer * nhar1];
      }
    }
  }
  //loop between decay and attack
  int counter = fullloop;
  while(counter > 0)
  {
    //decay to attack
  for(j = decayf; j > attackf; j--)
  {
    for(k = 1; k < nhar1; k++)
    {
      // if(amplowpassed[k][j] == 0) amplitudeScale = 1;
      // else amplitudeScale = 
      //   amplowpassednew[k][samplePointer]/amplowpassed[k][j];
      //(*cmag)[k + samplePointer*nhar1] = 
      //   cmagold[k + j*nhar1]*(amplowpassednew[k][samplePointer]/ 
      //     amplowpassed[k][j]);
      (*cmag)[k + samplePointer*nhar1] =
         amplowpassednew[k][samplePointer] + amporiginal[k][j] 
           - amplowpassed[k][j];
      (*dfr)[k + samplePointer*nhar1] = dfrold[k + j*nhar1];
    }
    samplePointer++;
  }
  //attack to decay
  for(j = attackf; j < decayf; j++)
  {
    for(k = 1; k < nhar1; k++)
    {
      // amplitudeScale = 0;
      // if (amplowpassed[k][j] == 0) amplitudeScale = 1;
      // else
      //   amplitudeScale = 
      //     amplowpassednew[k][samplePointer]/amplowpassed[k][j];
      //(*cmag)[k + samplePointer * nhar1] = cmagold[k + j * nhar1]
      //  *(amplowpassednew[k][samplePointer] / amplowpassed[k][j]);

      (*cmag)[k + samplePointer * nhar1] = 
          amplowpassednew[k][samplePointer] + amporiginal[k][j] 
            - amplowpassed[k][j];
      (*dfr)[k + samplePointer*nhar1] = dfrold[k + j * nhar1];
    }
    samplePointer++;
  }
  counter--;
  }
  //last round, decay to reverse to decay to end
  P("mduration is %f\n", mduration);
  float finalLoopT = extension - fullloop * mduration * 2;
  P ("final loop time is %f\n", finalLoopT);
  int reversef = roundf(decayf - 0.5 * finalLoopT / dt);
  P ("the reversef is %d\n", reversef);
  for(j = decayf; j > reversef; j--)
  {
  for(k = 1; k < nhar1; k++)
  {
  // amplitudeScale = 0;
  // if (amplowpassed[k][j] == 0) amplitudeScale = 1;
  // else
  //   amplitudeScale = 
  //     amplowpassednew[k][samplePointer] / amplowpassed[k][j];
    (*cmag)[k + samplePointer * nhar1] = 
       amplowpassednew[k][samplePointer] + amporiginal[k][j] 
         - amplowpassed[k][j];
    (*dfr)[k + samplePointer * nhar1] = dfrold[k + j * nhar1];
  }
  samplePointer++;
  }
  samplePointer -= 1;
  for(j = reversef; j < decayf; j++)
  {
    for(k = 1; k < nhar1; k++)
    {
      // amplitudeScale = 0;
      // if(amplowpassed[k][j] == 0) amplitudeScale = 1;
      // else
      //   amplitudeScale = 
      //     amplowpassednew[k][samplePointer] / amplowpassed[k][j];

      (*cmag)[k + samplePointer * nhar1] = 
         amplowpassednew[k][samplePointer] + amporiginal[k][j] 
           - amplowpassed[k][j];
      (*dfr)[k + samplePointer * nhar1] = dfrold[k + j * nhar1];
    }
    samplePointer++;
  }
  for(j = decayf; j < npts; j++)
  {
    for(k = 1; k < nhar1; k++)
    {
      // amplitudeScale = 0;
      // if (amplowpassed[k][j] == 0) amplitudeScale = 1;
      // else
      //  amplitudeScale = amplowpassednew[k][samplePointer] 
      //  /amplowpassed[k][j];
  
      (*cmag)[k + samplePointer * nhar1] = amporiginal[k][j];
      (*dfr)[k + samplePointer * nhar1] = dfrold[k + j * nhar1];
    }
    samplePointer++;
    if (samplePointer >= nptsnew) break;
  }
  P("the extension is done\n");
  P("nptsnew = %d, npts = %d\n", nptsnew, npts);

  P("\nAfter adding the extended zig-zagged microvariations\n");/*jwb 11/29/17*/
  P("to the extended smooth amplitude curve,\n");	    /* jwb 11/29/17 */
  P("Hit Return to plot the final harmonic %d amplitude curve: ",kplot);
  getchar();						    /* jwb 11/17/17 */
  float* HarampResult = (float*)calloc(nptsnew,sizeof(float));
  float* TimeResult = (float*)calloc(nptsnew,sizeof(float));
  for(i=0; i<nptsnew; i++)
  {
    HarampResult[i] = (*cmag)[kplot + i * nhar1];	    /* jwb 11/17/17 */
    TimeResult[i] = i*dt;
  }
  exlab1 = "final amplitude";				    /* jwb 11/28/17 */
  sprintf(graphlabel,"file: %s, %s", anfil,exlab1);	    /* jwb 11/28/17 */
  plotseg(TimeResult,HarampResult,nptsnew,"TIME (SEC)",svlabel);

  npts = nptsnew;
  tl = length;

  for (k = 1; k < nhar1; k++)
  {
    free(amporiginal[k]);
    free(amplowpassed[k]);
    free(amplowpassednew[k]);
  }

  free(amporiginal);
  free(amplowpassed);
  free(amplowpassednew);
} // end extendsyn()

void addsyn(char* filnam, int byte_reverse) /* jwb 12/06/99 */
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

//  P("\nOriginal time in input file:  ");		    /* jwb 12/05/94 */
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
//  if (tm < tmp) {P( "%.2f  ",t); fflush(stdout);}         /* jwb 12/05/94 */
    tmp = tm;					            /* jwb 12/05/94 */
  } /* end for(t -- sample generation loop */

  /* Adjust the file and data size values in the SND or WAVE header. */
  if (filetype == WAVE) 				   /* jjjm 02/10/98 */
    fixWavHdr(uno);					   /* jjjm 02/10/98 */
  else if (filetype == SND)                                 /* jwb 04/18/98 */
    fixSndHdr(uno);					   /* jjjm 02/10/98 */
  P("\n");						    /* jwb 03/01/00 */
} /* end addsyn() */

int getfiltype(char *name)
{ int i, len;
  len = strlen(name);
  for(i=0;i<MAXFTYPES;i++)
    if(!strcmp(&name[len-3],tail[i])) return(i);
  return(MAXFTYPES);
}
void getout()						    /* jwb 04/19/98 */
{							    /* jwb 04/19/98 */
  int i;						    /* jwb 04/19/98 */
  P("Sorry, only ");					    /* jwb 04/19/98 */
  for(i=0;i<MAXFTYPES-1;i++) P("%s, ", tail[i]);	    /* jwb 04/19/98 */
  P("and %s extensions are supported\n", tail[MAXFTYPES-1]);/* jwb 04/19/98 */
  exit(-1);						    /* jwb 04/19/98 */
}							    /* jwb 04/19/98 */
int plabel(double xpos,double ypos)
{
    g_move_abs(xpos,ypos);
    g_text(graphlabel);
    return(0);
}

