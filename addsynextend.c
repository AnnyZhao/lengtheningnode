/***************************************************************************
                                addsynextend.c

	Author: An Zhao  March, 2017

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
#include "header.h"


#define ERROR	(-1)
#define TABSIZ	5000
#define FLOAT 	0
#define SHORT	1
#define SND	0
#define WAVE     1
#define MAXFTYPES  2
#define P  printf

char *tail[MAXFTYPES] =  {"snd", "wav"};

/* function declarations */
char *gettime(); int getfiltype();
void cutSilence(float** cmag, float** dfr);
void calcRMS(float* cmag);
void extendsyn(float length, float extension, int ratio);
void reddursyn(float length, float origDur);
void addsyn(char* filnam, int byte_reverse);
float blend(float x);
void findattackdecay(int *attackf, int *decayf);
void getLPCoefficientsButterworth2Pole(const int samplerate, const double cutoff, double* const ax, double* const by);
void ButterworthFilter(float* samples, float* samplespassed, int count, int sampleRate);
/* global variables: */
HEADER header;
int nhar, nhar1, npts;
float *cmag, *dfr, *phase, dt, tl, smax, fs, fa, scalefac;

int main ( int argc, char** argv )
{

    char *anfil, *synfil;
    int byte_reverse = 0;
    if(argc < 3)
    {
      fprintf(stderr, "Usage: %s <soundfile> <output file>\n", argv[0]);
      P("try again!\n\n");				    /* jwb 02/04/17 */
      exit(1);
    }

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
    //cutSilence(&cmag, &dfr);
    P("the number of harmonics is %d, the number of frames is %d\n", nhar, npts);
    int aa;
    for (aa = 0; aa < 10; aa++)
    {
        P("0th harmonic value for frame %d is %f\n", aa, cmag[aa * nhar1]);
    }
    int attackf, decayf;
    calcRMS(cmag);
    for (aa = 0; aa < 10; aa++)
    {
        P("0th harmonic value for frame %d is %f\n", aa, cmag[aa * nhar1]);
    }

    float origDur = tl;
    P("original sound duration is %.3f sec:\n", origDur);

    findattackdecay(&attackf, &decayf);
    // decayf = decayf - (decayf - attackf) * 0.3;
    // attackf = attackf + (decayf - attackf) * 0.3;
    float attackt = attackf * dt;
    float decayt = decayf * dt;
    P ("attackt is %f and decayt is %f\n", attackt, decayt);

    float totalt, extendt;
    int ratio;
    float minimumT = attackt + (origDur - decayt);
    P ("minimum time is %f\n", minimumT);
    while (1)
    {
        P("Give new time duration: ");			    /* jwb 02/03/17 */
        scanf("%f", &totalt);
        if (totalt > minimumT) break;
        P("time is too short to reduce\n");	    /* jwb 02/03/17 */
        P("Try again.\n");
    }

    if (totalt > origDur)
    {
        decayf = decayf - (decayf - attackf) * 0.01;
        attackf = attackf + (decayf - attackf) * 0.01;
        float mduration = (decayf - attackf) * dt;
        P("mduration is %f\n", mduration);
        extendt = totalt - origDur;
        P("extendt is %f\n", extendt);
        ratio = extendt / mduration;
        P("ratio is %d\n",ratio);
        extendsyn(totalt, extendt, ratio);
        P("the extension is done");
    }

    if (totalt < origDur)
    {
        int f1, f2;
        float t1, t2;
        reddursyn(totalt, origDur);
        P("the reduction is done\n");
    }


    /*  resynthesize tone at sample rate fs */
    P("Begin synthesis\n");
    addsyn(synfil, byte_reverse);				    /* jwb 12/06/99 */
    P("\nSynthesis of file %s complete \n", 		    /* jwb 04/02/98 */
    			synfil);		    /* jwb 04/02/98 */
					    /* jwb 02/03/17 */
} /* end main() */


void calcRMS(float* cmag)
{
    int i, k;
    float sum = 0.0;
    for (i = 0 ; i < npts; i++)
    {
        for (k = 1; k < nhar1; k++)
        {
            sum += sq(cmag[k + i * nhar1]);
        }
        sum = sqrt(sum);
        cmag[i* nhar1] = sum;
    }
}

void cutSilence(float** cmag, float** dfr)
{
    calcRMS(*cmag);
    int i;
    int beginF = 0, endF = 0;
    float* dbarr = (float*)calloc(npts, sizeof(float));
    //find out the cut point
    for (i = 0; i < npts; i++)
    {
        dbarr[i] = 20.*log10f((*cmag)[i * nhar1]);
        if (i > 2000)
        {
            P("DB value for frame %d is %f\n", i , dbarr[i]);
        }
        if ((dbarr[i-1] < 30.0) && (dbarr[i] > 30.0) && beginF == 0)
        {
            beginF = i;
            P("begin frame is %d\n", i);
        }

        if ((dbarr[i-1] > 29.0) && (dbarr[i] < 29.0) && endF == 0)
        {
            endF = i - 1;
            P("end frame is %d\n", i-1);
        }
    }
    //update information
    npts = endF - beginF + 1;
    P("npts is now %d\n",npts);
    tl = npts * dt;
    float * cmagold = *cmag;
    float * dfrold = *dfr;
    *cmag = (float*)calloc(npts * nhar1, sizeof(float));
    *dfr = (float*)calloc(npts * nhar1, sizeof(float));
    int j = 0;
    for (i = beginF; i < endF + 1; i++)
    {
        int k;
        for (k = 1; k < nhar1; k++)
        {
            (*cmag)[k + j * nhar1] = cmagold[k + i * nhar1];
            (*dfr)[k + j * nhar1] = dfrold[k + i * nhar1];
        }
        j++;
    }
    //
}

void findattackdecay(int *attackf, int *decayf)
{
    //calculate the attack and decay
    //k = 0: rms amplitude
    //original no. of frames
    int i = 0;
    float* dbarr = (float*)calloc(npts+1, sizeof(float));
    P("npts in findattackdecay is %d\n",npts);
    float sumdb = 0;
    int npts_nonzero = 0;
    for (i = 0; i < npts; i++)
    {
        // if (i < 10)
        // {
        //     P("the first 0 harmonic in frame %d, %f\n",i, cmag[i * nhar1]);
        // }

        dbarr[i] = 20.*log10f(cmag[i * nhar1]);
        sumdb += dbarr[i];
        //if (dbarr[i] > 17.0) npts_nonzero += 1;
    }
    float avgdb = sumdb/npts;
    P("db average = %.1f\n",avgdb);

    //first pass calculation
    float firstdiff = 0, seconddiff = 0;

    printf("npts=%d\n",npts);
    for (i = 0; i < npts; i++)
    {
        //printf("i=%d\n",i);
      firstdiff = seconddiff;
      seconddiff = dbarr[i] - avgdb;
      if(seconddiff > 0 && firstdiff < 0)
      {
         *attackf = i;
         P("attack frame is %d\n", *attackf);
      }
      if(seconddiff < 0 && firstdiff > 0)
      {
         *decayf = i;
         P("decay frame is %d\n", *decayf);
      }
     //printf("i=%d\n",i);
    }
}

float blend(float x)
{
    return x;
}

void reddursyn(float length, float origDur)
{
    float x,w;
    int attackf, decayf;
    findattackdecay(&attackf, &decayf);
    // decayf = decayf - (decayf - attackf) * 0.3;
    // attackf = attackf + (decayf - attackf) * 0.3;
    int nptsnew = length/dt;
    P("nptsnew: %d\n", nptsnew);
    int overlapF = nptsnew - attackf - (npts - decayf);
    P("overlapF: %d\n", overlapF);
    int shiftF = npts - nptsnew;
    P("shiftF: %d\n",shiftF);
    x = 1./overlapF; //function factor

    float* cmagold = cmag;
    float* dfrold = dfr;

    cmag = (float*)calloc(nptsnew * nhar1, sizeof(float));
    dfr = (float*)calloc(nptsnew * nhar1, sizeof(float));
    P ("the size of cmag and dfr is %d\n", nptsnew * nhar1);

    int i, j, k, narg1, narg2;
    // beginning to attackf
    for (i = 0; i < attackf; i++)
    {
        for (k = 1; k < nhar1; k++)
        {
            cmag[k + i * nhar1] = cmagold[k + i * nhar1];
            dfr[k + i * nhar1] = dfrold[k + i * nhar1];
        }
    }
    //start to blend
    for (j = attackf; j < attackf + overlapF + 1; j++)
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
            cmag[k + i * nhar1] = (1-w) * cmagold[narg1] + w * cmagold[narg2];
            dfr[k + i * nhar1] = (1-w) * dfrold[narg1] + w * dfrold[narg2];
        }
        i++;
    }
    //decay to end
    for (j = decayf; j < npts; j++)
    {
        for (k = 1; k < nhar1; k++)
        {
            cmag[k + i * nhar1] = cmagold[k + j * nhar1];
            dfr[k + i * nhar1] = dfrold[k + j * nhar1];
        }
        i++;
    }
    //update parameters for addsyn
    npts = nptsnew;
    tl = length;

    //debug printing
    for (k = 1; k <= 5; k ++)
    {
        for (i = attackf; i < attackf + 20; i++)
        {
            narg1 = k + i * nhar1;
            P("cmag[%d] is %f, dfr[%d] is %f\n", narg1, cmag[narg1], narg1, dfr[narg1]);
        }
    }
}

// Butterworth filter from http://baumdevblog.blogspot.com/2010/11/butterworth-lowpass-filter-coefficients.html
// Cutoff is 40 Hz

void getLPCoefficientsButterworth2Pole(const int samplerate, const double cutoff, double* const ax, double* const by)
{
    double sqrt2 = 1.4142135623730950488;

    double QcRaw  = (2 * 3.141516 * cutoff) / samplerate; // Find cutoff frequency in [0..PI]
    double QcWarp = tan(QcRaw); // Warp cutoff frequency

    double gain = 1 / (1+sqrt2/QcWarp + 2/(QcWarp*QcWarp));
    by[2] = (1 - sqrt2/QcWarp + 2/(QcWarp*QcWarp)) * gain;
    by[1] = (2 - 2 * 2/(QcWarp*QcWarp)) * gain;
    by[0] = 1;
    ax[0] = 1 * gain;
    ax[1] = 2 * gain;
    ax[2] = 1 * gain;
}

void ButterworthFilter(float* samples, float* samplespassed, int count, int sampleRate)
{
    double xv[3];
    double yv[3];
    double ax[3];
    double by[3];

    getLPCoefficientsButterworth2Pole(sampleRate, 2, ax, by);

    for (int i=0;i<count;i++)
    {
        xv[2] = xv[1]; xv[1] = xv[0];
        xv[0] = samples[i];
        yv[2] = yv[1]; yv[1] = yv[0];

        yv[0] =   (ax[0] * xv[0] + ax[1] * xv[1] + ax[2] * xv[2]
                     - by[1] * yv[0]
                     - by[2] * yv[1]);
        samplespassed[i] = yv[0];
    }
}

// void extendsyn(float** cmag, float** dfr, int nhar1, float length, float extension, float frameDuration)
// {
//     int i = 0;
//     int k = 0;
//     int attackf, decayf;
//     findattackdecay(&attackf, &decayf);
//     decayf = decayf - (decayf - attackf) * 0.05;
//     attackf = attackf + (decayf - attackf) * 0.05;
//
//     float mduration = (decayf - attackf) * dt;
//     float ratio = extension / mduration;
//
//     int nptsnew = length / dt;
//     P("nptsnew: %d\n", nptsnew);
//     int extendframes = nptsnew - npts;
//
//     int frameRate = 1.0 / frameDuration;
//
//     float* cmagold = *cmag;
//     float* dfrold = *dfr;
//     // calculate and interpolate amplitude variation
//
// //original data
//     float** amporiginal = (float**) calloc(nhar1, sizeof(float*));
//     float** amplowpassed = (float**) calloc(nhar1, sizeof(float*));
//     float** microvaria = (float**) calloc(nhar1, sizeof(float*));
// //new data
//     float** amplowpassednew = (float**) calloc(nhar1, sizeof(float*));
//     float** microvarialooped = (float**) calloc(nhar1, sizeof(float*));
//
//     for (k = 1; k < nhar1; k++)
//     {
//     //original data
//         amporiginal[k] = (float*) calloc(npts, sizeof(float));
//         amplowpassed[k] = (float*) calloc(npts, sizeof(float));
//         microvaria[k] = (float*) calloc(npts, sizeof(float));
//     //new data
//         amplowpassednew[k] = (float*) calloc(nptsnew, sizeof(float));
//         microvarialooped[k] = (float*) calloc(nptsnew, sizeof(float));
//     }
//
// // doing lowpass
//     for (k = 1; k < nhar1; k++)
//     {
//         for (i = 0; i < npts; i++)
//         {
//             amporiginal[k][i] = (*cmag)[k + i * nhar1];
//         }
//         // apply low-pass filtering to dboriginal
//         // to permit some variation in db level deviating from original
//         ButterworthFilter(amporiginal[k], amplowpassed[k], npts, frameRate);
//     }
//
//
// //shift amplowpassed
// int shift;
// int shiftamount = frameRate / 2;
//     for (k = 1; k < nhar1; k++)
//     {
//
//         for (i = shiftamount; i < npts; i++)
//         {
//             shift = i - shiftamount;
//             amplowpassed[k][shift] = amplowpassed[k][i];
//         }
//         for (shift = npts - shiftamount; shift < npts; shift++)
//         {
//             amplowpassed[k][shift] = amplowpassed[k][npts - shiftamount - 1];
//         }
//     }
// //time scale amplowpassed
//     for (k = 1; k < nhar1; k++)
//     {
//         for (i = 0; i < nptsnew; i++)
//         {
//             int intoSustain = i - attackf;
//             if (intoSustain <= 0)
//             {
//                 amplowpassednew[k][i] = amporiginal[k][i];
//             }
//             else if (i - extendframes >= decayf)
//             {
//                 amplowpassednew[k][i] = amporiginal[k][i - extendframes];
//             }
//             else
//             {
//                 float fullratio = ratio + 1;
//                 int a = floor((((float) i) - attackf) / fullratio) + attackf;
//                 int b = ceil((((float) i) - attackf) / fullratio) + attackf;
//                 float percentage = ((((float) i) - attackf) / fullratio) + attackf - a;
//                 if (a >= npts)
//                 {
//                     a = npts - 1;
//                 }
//                 if (b >= npts)
//                 {
//                     b = npts - 1;
//                 }
//                 amplowpassednew[k][i] = amporiginal[k][b] * percentage + amporiginal[k][a] * (1 - percentage);
//             }
//         }
//     }
//
//     //loop microvaria
//     *cmag = (float*)calloc(nptsnew * nhar1, sizeof(float));
//     *dfr = (float*)calloc(nptsnew * nhar1, sizeof(float));
//
//     int samplePointer = 0;
//     int fullloop = ratio / 2;
//     P("number of full loops: %d\n", fullloop);
//     // first beginning to decayf
//     int j;
//     for (samplePointer = 0; samplePointer < decayf; samplePointer++)
//     {
//         for (k = 1; k < nhar1; k++)
//         {
//             if (samplePointer <= attackf)
//             {
//                 //microvarialooped[k][samplePointer] = cmagold[k + samplePointer * nhar1];
//                 (*cmag)[k + samplePointer * nhar1] = cmagold[k + samplePointer * nhar1];
//             }
//             else
//             {
//                 // float result = pow(10.0, ((20.0 * logamplitude + dbscalefactor[k][samplePointer] - dbnewscalefactor[k][samplePointer]) / 20.0));
//                 (*cmag)[k + samplePointer * nhar1] = cmagold[k + samplePointer * nhar1] / amplowpassed[k][samplePointer] * amplowpassednew[k][samplePointer];
//             }
//             (*dfr)[k + samplePointer * nhar1] = dfrold[k + samplePointer * nhar1];
//         }
//     }
//     //loop between decay and attack
//     int counter = fullloop;
//     while (counter > 0)
//     {
//         //decaty to attack
//         for (j = decayf; j > attackf; j--)
//         {
//             for (k = 1; k < nhar1; k++)
//             {
//                 //float logamplitude = log10(cmagold[k + j * nhar1]);
//                 (*cmag)[k + samplePointer * nhar1] = cmagold[k + j * nhar1] / amplowpassed[k][j] * amplowpassednew[k][samplePointer];
//                 (*dfr)[k + samplePointer * nhar1] = dfrold[k + j * nhar1];
//             }
//             samplePointer++;
//         }
//         //attack to decay
//         for (j = attackf; j < decayf; j++)
//         {
//             for (k = 1; k < nhar1; k++)
//             {
//                 (*cmag)[k + samplePointer * nhar1] = cmagold[k + j * nhar1] / amplowpassed[k][j] * amplowpassednew[k][samplePointer];
//                 (*dfr)[k + samplePointer * nhar1] = dfrold[k + j * nhar1];
//             }
//             samplePointer++;
//         }
//         counter--;
//     }
//
//     //last round, decay to reverse to end
//     P("mduration is %f\n", mduration);
//     float finalLoopT = extension - fullloop * mduration * 2;
//     P ("final loop time is %f\n", finalLoopT);
//     int reversef = decayf - 0.5 * finalLoopT / dt;
//     P ("the reversef is %d\n", reversef);
//     for (j = decayf; j > reversef; j--)
//     {
//         for (k = 1; k < nhar1; k++)
//         {
//             (*cmag)[k + samplePointer * nhar1] = cmagold[k + j * nhar1] / amplowpassed[k][j] * amplowpassednew[k][samplePointer];
//             (*dfr)[k + samplePointer * nhar1] = dfrold[k + j * nhar1];
//         }
//         samplePointer++;
//     }
//     for (j = reversef; j < npts; j++)
//     {
//         for (k = 1; k < nhar1; k++)
//         {
//             (*cmag)[k + samplePointer * nhar1] = cmagold[k + j * nhar1] / amplowpassed[k][j] * amplowpassednew[k][samplePointer];
//             (*dfr)[k + samplePointer * nhar1] = dfrold[k + j * nhar1];
//         }
//         samplePointer++;
//         if (samplePointer >= nptsnew)
//             break;
//     }
//
//     npts = nptsnew;
//     tl = length;
//
//     for (k = 1; k < nhar1; k++)
//     {
//         free(amporiginal[k]);
//         free(amplowpassed[k]);
//         free(amplowpassednew[k]);
//     }
//
//     free(amporiginal);
//     free(amplowpassed);
//     free(amplowpassednew);
// }

void extendsyn(float length, float extension, int ratio)
{
    int i = 0;
    int attackf, decayf;
    findattackdecay(&attackf, &decayf);
    decayf = decayf - (decayf - attackf) * 0.43;
    attackf = attackf + (decayf - attackf) * 0.43;
    //
    float* cmagold = cmag;
    float* dfrold = dfr;
    int nptsnew = length/dt;
    P("nptsnew: %d\n", nptsnew);
    cmag = (float*)calloc(nptsnew * nhar1, sizeof(float));
    dfr = (float*)calloc(nptsnew * nhar1, sizeof(float));
    int fullloop;
    fullloop = ratio / 2;
    P("number of full loops: %d\n",fullloop);
    // first beginning to decayf
    int j,k;
    for (i = 0; i < decayf; i++)
    {
        for (k = 1; k < nhar1; k++)
        {
            cmag[k + i * nhar1] = cmagold[k + i * nhar1];
            dfr[k + i * nhar1] = dfrold[k + i * nhar1];
        }
    }
    //loop between decay and attack
    int counter = fullloop;
    while (counter > 0)
    {
        //decaty to attack
        for (j = decayf; j > attackf; j--)
        {
            for (k = 1; k < nhar1; k++)
            {
                cmag[k + i * nhar1] = cmagold[k + j * nhar1];
                dfr[k + i * nhar1] = dfrold[k + j * nhar1];
            }
            i++;
        }
        //attack to decay
        for (j = attackf; j < decayf; j++)
        {
            for (k = 1; k < nhar1; k++)
            {
                cmag[k + i * nhar1] = cmagold[k + j * nhar1];
                dfr[k + i * nhar1] = dfrold[k + j * nhar1];
            }
            i++;
        }
        counter--;
    }

    //last round, decay to reverse to end
    float mduration = (decayf - attackf) * dt;
    P("mduration is %f\n", mduration);
    float finalLoopT = extension - fullloop * mduration * 2;
    P ("final loop time is %f\n", finalLoopT);
    int reversef = decayf - 0.5 * finalLoopT / dt;
    P ("the reversef is %d\n", reversef);
    for (j = decayf; j > reversef; j--)
    {
        for (k = 1; k < nhar1; k++)
        {
            cmag[k + i * nhar1] = cmagold[k + j * nhar1];
            dfr[k + i * nhar1] = dfrold[k + j * nhar1];
        }
        i++;
    }
    for (j = reversef; j < npts; j++)
    {
        for (k = 1; k < nhar1; k++)
        {
            cmag[k + i * nhar1] = cmagold[k + j * nhar1];
            dfr[k + i * nhar1] = dfrold[k + j * nhar1];
        }
        i++;
    }

    npts = nptsnew;
    tl = length;
}

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
  // if(filetype == AU)                                        /* jwb 04/15/98 */
  // {                                                         /* jwb 04/15/98 */
  //   P("\nThis program does not support the .au extension.\n");
  //   P("Extension .snd will be used instead\n");             /* jwb 04/15/98 */
  //   filetype == SND;                                        /* jwb 04/18/98 */
  // }                                                         /* jwb 04/15/98 */
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
