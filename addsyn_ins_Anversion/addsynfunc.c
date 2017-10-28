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
        if (dbarr[i] <= 0)
        {
            dbarr[i] = 0.0001;
        }
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
        if (i >= nptsnew) break;
    }
    //update parameters for addsyn
    npts = nptsnew;
    tl = length;

    //debug printing
}

// Butterworth filter from http://baumdevblog.blogspot.com/2010/11/butterworth-lowpass-filter-coefficients.html
// Cutoff is 40 Hz

void getLPCoefficientsButterworth2Pole(const int samplerate, const double cutoff, double* const ax, double* const by)
{
    double sqrt2 = 1.4142135623730950488;

    double QcRaw  = (2 * PI * cutoff) / samplerate; // Find cutoff frequency in [0..PI]
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

void extendsyn(float** cmag, float** dfr, int nhar1, float length, float extension, float frameDuration)
{
    int i = 0;
    int k = 0;
    int attackf, decayf;
    findattackdecay(&attackf, &decayf);
    decayf = decayf - (decayf - attackf) * 0.05;
    attackf = attackf + (decayf - attackf) * 0.05;

    float mduration = (decayf - attackf) * dt;
    float ratio = extension / mduration;

    int nptsnew = length / dt;
    P("nptsnew: %d\n", nptsnew);
    int extendframes = nptsnew - npts;

    int frameRate = 1.0 / frameDuration;

    float* cmagold = *cmag;
    float* dfrold = *dfr;
    // calculate and interpolate amplitude variation

//original data
    float** amporiginal = (float**) calloc(nhar1, sizeof(float*));
    float** amplowpassed = (float**) calloc(nhar1, sizeof(float*));
    float** microvaria = (float**) calloc(nhar1, sizeof(float*));
//new data
    float** amplowpassednew = (float**) calloc(nhar1, sizeof(float*));
    float** microvarialooped = (float**) calloc(nhar1, sizeof(float*));

    for (k = 1; k < nhar1; k++)
    {
    //original data
        amporiginal[k] = (float*) calloc(npts, sizeof(float));
        amplowpassed[k] = (float*) calloc(npts, sizeof(float));
        microvaria[k] = (float*) calloc(npts, sizeof(float));
    //new data
        amplowpassednew[k] = (float*) calloc(nptsnew, sizeof(float));
        microvarialooped[k] = (float*) calloc(nptsnew, sizeof(float));
    }

// doing lowpass
    for (k = 1; k < nhar1; k++)
    {
        for (i = 0; i < npts; i++)
        {
            amporiginal[k][i] = (*cmag)[k + i * nhar1];
        }
        // apply low-pass filtering to dboriginal
        // to permit some variation in db level deviating from original
        ButterworthFilter(amporiginal[k], amplowpassed[k], npts, frameRate);
    }


//shift amplowpassed
int shift;
int shiftamount = 1./sqrt(2.) * 4. * atan(1.) * 2 * frameDuration;
    for (k = 1; k < nhar1; k++)
    {

        for (i = shiftamount; i < npts; i++)
        {
            shift = i - shiftamount;
            amplowpassed[k][shift] = amplowpassed[k][i];
        }
        for (shift = npts - shiftamount; shift < npts; shift++)
        {
            amplowpassed[k][shift] = amplowpassed[k][npts - shiftamount - 1];
        }
    }
//time scale amplowpassed
    for (k = 1; k < nhar1; k++)
    {
        for (i = 0; i < nptsnew; i++)
        {
            int intoSustain = i - attackf;
            if (intoSustain <= 0)
            {
                amplowpassednew[k][i] = amporiginal[k][i];
            }
            else if (i - extendframes >= decayf)
            {
                amplowpassednew[k][i] = amporiginal[k][i - extendframes];
            }
            else
            {
                float fullratio = ratio + 1;
                int a = floor((((float) i) - attackf) / fullratio) + attackf;
                int b = ceil((((float) i) - attackf) / fullratio) + attackf;
                float percentage = ((((float) i) - attackf) / fullratio) + attackf - a;
                if (a >= npts)
                {
                    a = npts - 1;
                }
                if (b >= npts)
                {
                    b = npts - 1;
                }
                amplowpassednew[k][i] = amporiginal[k][b] * percentage + amporiginal[k][a] * (1 - percentage);
            }
        }
    }

    //loop microvaria
    *cmag = (float*)calloc(nptsnew * nhar1, sizeof(float));
    *dfr = (float*)calloc(nptsnew * nhar1, sizeof(float));

    int samplePointer = 0;
    int fullloop = ratio / 2;
    P("number of full loops: %d\n", fullloop);
    // first beginning to decayf
    int j;
    for (samplePointer = 0; samplePointer < decayf; samplePointer++)
    {
        for (k = 1; k < nhar1; k++)
        {
            if (samplePointer <= attackf)
            {
                //microvarialooped[k][samplePointer] = cmagold[k + samplePointer * nhar1];
                (*cmag)[k + samplePointer * nhar1] = cmagold[k + samplePointer * nhar1];
            }
            else
            {
                // float result = pow(10.0, ((20.0 * logamplitude + dbscalefactor[k][samplePointer] - dbnewscalefactor[k][samplePointer]) / 20.0));
                (*cmag)[k + samplePointer * nhar1] = cmagold[k + samplePointer * nhar1] / amplowpassed[k][samplePointer] * amplowpassednew[k][samplePointer];
            }
            (*dfr)[k + samplePointer * nhar1] = dfrold[k + samplePointer * nhar1];
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
                //float logamplitude = log10(cmagold[k + j * nhar1]);
                (*cmag)[k + samplePointer * nhar1] = cmagold[k + j * nhar1] / amplowpassed[k][j] * amplowpassednew[k][samplePointer];
                (*dfr)[k + samplePointer * nhar1] = dfrold[k + j * nhar1];
            }
            samplePointer++;
        }
        //attack to decay
        for (j = attackf; j < decayf; j++)
        {
            for (k = 1; k < nhar1; k++)
            {
                (*cmag)[k + samplePointer * nhar1] = cmagold[k + j * nhar1] / amplowpassed[k][j] * amplowpassednew[k][samplePointer];
                (*dfr)[k + samplePointer * nhar1] = dfrold[k + j * nhar1];
            }
            samplePointer++;
        }
        counter--;
    }

    //last round, decay to reverse to end
    P("mduration is %f\n", mduration);
    float finalLoopT = extension - fullloop * mduration * 2;
    P ("final loop time is %f\n", finalLoopT);
    int reversef = decayf - 0.5 * finalLoopT / dt;
    P ("the reversef is %d\n", reversef);
    for (j = decayf; j > reversef; j--)
    {
        for (k = 1; k < nhar1; k++)
        {
            (*cmag)[k + samplePointer * nhar1] = cmagold[k + j * nhar1] / amplowpassed[k][j] * amplowpassednew[k][samplePointer];
            (*dfr)[k + samplePointer * nhar1] = dfrold[k + j * nhar1];
        }
        samplePointer++;
    }
    for (j = reversef; j < npts; j++)
    {
        for (k = 1; k < nhar1; k++)
        {
            (*cmag)[k + samplePointer * nhar1] = cmagold[k + j * nhar1] / amplowpassed[k][j] * amplowpassednew[k][samplePointer];
            (*dfr)[k + samplePointer * nhar1] = dfrold[k + j * nhar1];
        }
        samplePointer++;
        if (samplePointer >= nptsnew)
            break;
    }

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
}
