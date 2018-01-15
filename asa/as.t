.comment        as.t    Additive Synthesis Instrument
.prefix
ADDSYN

.scorecard
    filename    (0 0 string)
    pch        (4 12 octpitch 8)
    ascale        (.1 10 log 1)
    lfrac        (0 1 lfrac 1)

.instrument
    float sifaci, ascalei, fai, *cmagi, *dfri, *phasei, dxi, xi, lfraci;
    float dti, tli, nptsi;
    int nhar1i, fptri, itrap, inoteno;

.globals
extern int nhar1, npts;
extern float fa, *phas0, *cmag, *dfr;
#include <strings.h>
#include "anread.c"
#include "addsynfunc.c"

.note
int k;
float freq, fac, OFF;
char *anFile = GetScorecardString(filename);

fptri = anread(anFile);
if(fptri == -1)
{
    P("Can not open file %s. Abort job.\n",anFile);
    exit(-1);
}
P("\nAdd syn of file %s at start_time = %.3f, duration = %.3f\n",
    anFile,STIME,DUR);
P("pitch = %.3f, ampscale= %.3f\n", pch, ascale); fflush(stdout);


freq = cpspitch(pch);
ascalei = ascale;
sifaci = (FFUN_LEN/SR)*(freq/fa);    /* provides frequency scaling */
phasei = (float *) malloc(nhar1 * sizeof(float));
fac = FFUN_LEN/(8. * atan(1)); OFF = .25 * FFUN_LEN;
for (k = 1; k < nhar1; k++) {
    phasei[k] = OFF + fac * phas0[k];
}
/* provides time-scaling: */

//start to Modify
float frameDuration = dt;
float origDur = tl;
P("original sound duration is %.3f sec:\n", origDur);
// increment for traversing the cmagi and dfri data
dxi = (npts - 1) / ((float) origDur * SR);
fflush(stdout);
calcRMS(cmag);
int attackf, decayf;
findattackdecay(&attackf, &decayf);

float attackt = attackf * dt;
float decayt = decayf * dt;
P("dt is %f\n", dt);
P("attackt is %f and decayt is %f\n", attackt, decayt);

float totalt, extendt;
totalt = DUR;
float ratio;
float minimumT = attackt + (origDur - decayt);
P("minimum duration is %.3f\n", minimumT);
fflush(stdout);

cmagi = cmag;
dfri = dfr;

//extend
if (DUR > origDur)
{
    decayf = decayf - (decayf - attackf) * 0.05;
    attackf = attackf + (decayf - attackf) * 0.05;
    float mduration = (decayf - attackf) * dt;
    P("mduration is %f\n", mduration);
    fflush(stdout);
    extendt = totalt - origDur;
    P("extendt is %f\n", extendt);
    fflush(stdout);
    ratio = extendt / mduration;
    P("ratio is %f\n", ratio);
    fflush(stdout);
    extendsyn(&cmagi, &dfri, nhar1, DUR, extendt, frameDuration);
    P("the extension is done\n");
}
//shorten
else if (DUR < origDur)
{
    if (DUR < minimumT) // timescale compressed attack+decay
    {
        timescalemin(&cmagi, &dfri, nhar1, DUR, minimumT, origDur);
    }
    else
    {
        reddursyn(&cmagi, &dfri, nhar1, DUR, origDur);
    }
}
// do nothing if DUR = origDur

dxi = (npts - 1) / ((float) DUR * SR);
printf("dxi: %f\n", dxi);

/* save parameters for sample comp */
fai = fa;
nhar1i = nhar1;
nptsi = npts;
dti = dt;
tli = tl;
xi = 0.;
lfraci = lfrac;
itrap = 0;
inoteno++;
P("Complete addsyn setup\n");
fflush(stdout);

.endnote
free(phasei);

.sample
int i, k, narg1, narg2, iphasek;
float sum, w1,w2,x1,x2, phasek, ampk;
// if(itrap==0)
// {
//     P("in .sample, noteno = %d\n", inoteno);
// }
// if (i>2170)
// {
//     P("i is %d\n",i);
// }
// if(i >= npts)
// {
//   P("i is greater then npts, it is going to crash\n");
// }
//   compute interpolation weights from current xi
i = xi; w2 = xi - i; w1 = 1.0 - w2;
if ((i + 1) < nptsi)
{
    narg1 = i * nhar1i;
    narg2 = (i + 1) * nhar1i;
    sum = 0;                /* add up to form sample value sum */
    for (k = 1; k < nhar1i; k++)
    {
        /* compute Sine table lookup phase for each harmonic */
        phasek = phasei[k];
        phasek += sifaci * (k * fai + w1 * dfri[narg1 + k] + w2 * dfri[narg2 + k]);
        phasek = amod(phasek, FFUN_LEN);
        if (phasek < 0.) phasek += FFUN_LEN;
        // save phase for next sample
        phasei[k] = phasek;
        iphasek = phasek; x2 = phasek - iphasek; x1 = 1.0 - x2;
        ampk = w1*cmagi[narg1 + k] + w2*cmagi[narg2 + k];
        sum += ampk*(x1*Sine[iphasek] + x2*Sine[iphasek +1]);
    }
}
else sum = 0.;
NoQuad(ascalei * sum, lfraci);        /* output sample */
xi += dxi;                /* update an array lookup index */
