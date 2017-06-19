.comment		as.t    Additive Synthesis Instrument
.prefix
ADDSYN

.scorecard
	filename	(0 0 string)
	pch		(4 12 octpitch 8)
	ascale		(.1 10 log 1)
	lfrac		(0 1 lfrac 1)

.instrument
  float sifaci, ascalei, fai, *cmagi, *dfri, *phasei, dxi, xi, lfraci;
  int nhar1i, fptri, itrap, inoteno;


.globals

#include "anread.c"
#include "macro.h"

float sq(float x) { return x * x; }

#include "addsynfunc.c"

/******************************************************************************
	12/04/96 jwb	fix click at end of some notes
 *****************************************************************************/
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

//

freq = cpspitch(pch); ascalei = ascale;
sifaci = (FFUN_LEN/SR)*(freq/fa);	/* provides frequency scaling */
phasei = (float *) malloc(nhar1*sizeof(float));
fac = FFUN_LEN/(8.*atan(1)); OFF = .25*FFUN_LEN;
for(k=1;k<nhar1;k++) phasei[k] = OFF + fac*phas0[k];
/* provides time-scaling: */
dxi = (npts-1)/((float)SR);				    /* jwb 12/04/96 */
/* save parameters for sample comp */
fai = fa;
cmagi = cmag;
dfri = dfr;
nhar1i = nhar1;
xi = 0.;
lfraci = lfrac;

//start to Modify
float origDur = tl;
P("original sound duration is %.3f sec:\n", origDur);
calcRMS(cmagi);
int attackf, decayf;
findattackdecay(&attackf, &decayf);

float attackt = attackf * dt;
float decayt = decayf * dt;
P ("dt is %f\n", dt);
P ("attackt is %f and decayt is %f\n", attackt, decayt);

float totalt, extendt;
//totalt = DUR;
int ratio;
float minimumT = attackt + (origDur - decayt);
P ("minimum time is %f\n", minimumT);

while (DUR < minimumT)
{
   DUR = minimumT;
}
//extend
if (totalt > origDur)
{
	decayf = decayf - (decayf - attackf) * 0.3;
	attackf = attackf + (decayf - attackf) * 0.3;
	float mduration = (decayf - attackf) * dt;
	P("mduration is %f\n", mduration);
	extendt = totalt - origDur;
	P("extendt is %f\n", extendt);
	ratio = extendt / mduration;
	P("ratio is %d\n",ratio);
	extendsyn(cmagi, dfri, nhar1i, DUR, extendt, ratio);
	P("the extension is done");
}
//shorten
if (totalt < origDur)
{
	reddursyn(cmagi, dfri, nhar1i, DUR, origDur);
}
itrap = 0;
inoteno++;
P("Complete AddSyn setup\n"); fflush(stdout);

.endnote

free(phasei);

.sample
int i, k, narg1, narg2, iphasek;
float sum, w1,w2,x1,x2, phasek, ampk;
if(itrap==0)
{
	P("in .sample, noteno = %d\n", inoteno);
}
i = xi; w2 = xi - i; w1 = 1.0 - w2;	/* compute interpolation weights */
if (i>2170)
{
	P("i is %d\n",i);
}
if (i >= npts)
{
	P("i is greater then npts, it is going to crash\n");
}
narg1 = i*nhar1i; narg2 = (i+1)*nhar1i;
sum = 0;				/* add up to form sample value sum */
for(k=1;k<nhar1i;k++)
{	/* compute Sine table lookup phase for each harmonic */
  phasek = phasei[k];
  phasek += sifaci*(k*fai + w1* dfri[narg1 + k] + w2*dfri[narg2 + k]);
  phasek = amod(phasek, FFUN_LEN);
  if(phasek < 0.) phasek += FFUN_LEN;
  phasei[k] = phasek;
  iphasek = phasek; x2 = phasek - iphasek; x1 = 1.0 - x2;
  ampk = w1*cmagi[narg1 + k] + w2*cmagi[narg2 + k];
  sum +=  ampk*(x1*Sine[iphasek] + x2*Sine[iphasek +1]);
  itrap++;
}

NoQuad(ascalei*sum, lfraci);		/* output sample */
xi += dxi;				/* update an array lookup index */
