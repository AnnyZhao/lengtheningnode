#include <math.h>
#include <stdio.h>

#define ERROR	(-1)
#define TABSIZ	5000
#define FLOAT 	0
#define SHORT	1
#define SND	0
#define WAVE     1
#define MAXFTYPES  2
#define P  printf

extern float dt, tl;
extern float* cmag;

float sq(float x) { return x * x; }

void reddursyn(float** Cmag, float** Dfr, int nhar1, float length, float origDur);

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

void findattackdecay(int *attackf, int *decayf)
{
    //calculate the attack and decay
    //k = 0: rms amplitude
    //original no. of frames
    int i = 0;
    int debugattack = 10000;
    float* dbarr = (float*)calloc(npts+1, sizeof(float));
    P("npts in findattackdecay is %d\n", npts);
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
         if (debugattack > i)
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
     //printf("i=%d\n",i);
    }
}

float blend(float x)
{
    return x;
}

float interpolation(float start, float end, float percentage)
{
    return start + (end - start) * percentage;
}

void timescalemin(float** Cmag, float** Dfr, int nhar1, float length, float mintime, float origDur)
{
    if (length > mintime)
        return; // not suitable for timescale
    reddursyn(Cmag, Dfr, nhar1, mintime * 1.01, origDur);     //modify to minimumT

    float* cmagold = *Cmag;
    float* dfrold = *Dfr;

    int nptsnew = length / dt;
    *Cmag = (float*) calloc(nptsnew * nhar1, sizeof(float));
    *Dfr = (float*) calloc(nptsnew * nhar1, sizeof(float));
    P ("the size of cmag and dfr is %d\n", nptsnew * nhar1); // create memory space

    float coeff;
    coeff = mintime / length; // calculate number of samples to jump

    //interpolate between each sample
    int i, k;
    for (i = 0; i < nptsnew; i++)
    {
        float ptOld = i * coeff;
        int pt1 = floor(ptOld);
        int pt2 = ceil(ptOld);
        float percentage = ptOld - pt1;
        if (pt2 >= npts) pt2 = npts - 1;
        for (k = 1; k < nhar1; k++)
        {
            (*Cmag)[k + i * nhar1] = interpolation(cmagold[k + pt1 * nhar1],
                                                   cmagold[k + pt2 * nhar1],
                                                   percentage);
            (*Dfr)[k + i * nhar1] = interpolation(dfrold[k + pt1 * nhar1],
                                                  dfrold[k + pt2 * nhar1],
                                                  percentage);
        }
    }

    tl = length;
    npts = nptsnew;
}

void reddursyn(float** Cmag, float** Dfr, int nhar1, float length, float origDur)
{
    float x,w;
    int attackf, decayf;
    findattackdecay(&attackf, &decayf);
    // decayf = decayf - (decayf - attackf) * 0.3;
    // attackf = attackf + (decayf - attackf) * 0.3;
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
    P ("the size of cmag and dfr is %d\n", nptsnew * nhar1);

    int i, j, k, narg1, narg2;
    // beginning to attackf
    for (i = 0; i < attackf; i++)
    {
        for (k = 1; k < nhar1; k++)
        {
            (*Cmag)[k + i * nhar1] = cmagold[k + i * nhar1];
            (*Dfr)[k + i * nhar1] = dfrold[k + i * nhar1];
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
            (*Cmag)[k + i * nhar1] = (1-w) * cmagold[narg1] + w * cmagold[narg2];
            (*Dfr)[k + i * nhar1] = (1-w) * dfrold[narg1] + w * dfrold[narg2];
        }
        i++;
    }
    //decay to end
    for (j = decayf; j < npts; j++)
    {
        for (k = 1; k < nhar1; k++)
        {
            (*Cmag)[k + i * nhar1] = cmagold[k + j * nhar1];
            (*Dfr)[k + i * nhar1] = dfrold[k + j * nhar1];
        }
        i++;
    }
    //update parameters for addsyn
    npts = nptsnew;
    tl = length;

    //debug printing
}


void extendsyn(float** cmag, float** dfr, int nhar1, float length, float extension, int ratio)
{
    int i = 0;
    int attackf, decayf;
    findattackdecay(&attackf, &decayf);
    decayf = decayf - (decayf - attackf) * 0.3;
    attackf = attackf + (decayf - attackf) * 0.3;
    //
    float* cmagold = *cmag;
    float* dfrold = *dfr;
    int nptsnew = length/dt;
    P("nptsnew: %d\n", nptsnew);
    *cmag = (float*)calloc(nptsnew * nhar1, sizeof(float));
    *dfr = (float*)calloc(nptsnew * nhar1, sizeof(float));
    int fullloop;
    fullloop = ratio / 2;
    P("number of full loops: %d\n",fullloop);
    // first beginning to decayf
    int j,k;
    for (i = 0; i < decayf; i++)
    {
        for (k = 1; k < nhar1; k++)
        {
            (*cmag)[k + i * nhar1] = cmagold[k + i * nhar1];
            (*dfr)[k + i * nhar1] = dfrold[k + i * nhar1];
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
                (*cmag)[k + i * nhar1] = cmagold[k + j * nhar1];
                (*dfr)[k + i * nhar1] = dfrold[k + j * nhar1];
            }
            i++;
        }
        //attack to decay
        for (j = attackf; j < decayf; j++)
        {
            for (k = 1; k < nhar1; k++)
            {
                (*cmag)[k + i * nhar1] = cmagold[k + j * nhar1];
                (*dfr)[k + i * nhar1] = dfrold[k + j * nhar1];
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
            (*cmag)[k + i * nhar1] = cmagold[k + j * nhar1];
            (*dfr)[k + i * nhar1] = dfrold[k + j * nhar1];
        }
        i++;
    }
    for (j = reversef; j < npts; j++)
    {
        for (k = 1; k < nhar1; k++)
        {
            (*cmag)[k + i * nhar1] = cmagold[k + j * nhar1];
            (*dfr)[k + i * nhar1] = dfrold[k + j * nhar1];
        }
        i++;
    }

    npts = nptsnew;
    tl = length;
}
