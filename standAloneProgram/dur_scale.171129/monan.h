/******************************************************************************    
 *				monan.h
 *
 * 			includes and defines for the monan routines
 *
 *	04/06/92 jwb	change P definition from fprintf to printf
 *	08/26/93 jwb	insert NHARMIN define
 *	11/01/94 jwb	change BRTHRESH to float value
 *	04/10/96 tjm    add convert_3D_to_2D()
 *	06/03/96 jwb	color_on, color global variables
 *	02/21/99 jwb	Introduce auxiliary analysis variables.
 *			smax, nchans, & ampscale eliminated from extern list.
 *	03/04/03 jwb	Eliminate size definition for filnam and elabel.
 *			These are defined in mainmenu.c
 *	11/16/07 jwb	Change array time to timet to avoid system conflict.
 *	06/15/08 jwb	Include stdlib.h and strings.h
 *****************************************************************************/
#include "stdlib.h"
#include <stdio.h>
#include <math.h>
#include <strings.h>
#include "g_raph.h"
#include "macro.h"
#include "header.h"
#define P printf
#define FP fprintf					    /* jwb 02/20/99 */
#define BRTHRESH 100.
#define NHARMIN 5					    /* jwb 08/26/93 */

/* function type declarations  */
void plotat(),plotaa(),plotbt(),plotft(),plotfi();
void plabel(double,double);				    /* jwb 04/27/94 */

void convert_3D_to_2D(					    /* tjm 04/10/96 */
 	float x,
	float y,
	float z,
	float *x2d,
	float *y2d,
	float yaw,
	float pitch,
	float roll,
	float size_factor,
        int perspective);

/*    global variables:     */
/* main analysis data */
extern HEADER header;
extern int nhar, nhar1, npts;				     /* jwb 02/21/99 */
extern float *cmag, *dfr, *phase, *br, *timet, tl, dt, fa;   /* jwb 11/16/07 */
/* auxiliary analysis data: */
extern HEADER headera;					     /* jwb 02/21/99 */
extern int nhara, nhar1a, nchansa, nptsa, auxfile;	     /* jwb 02/21/99 */
extern float *cmaga, *dfra, *phasea, *bra, tla, dta, faa;    /* jwb 02/21/99 */

extern char filnam[];					     /* jwb 03/04/03 */
extern char elabel[];					     /* jwb 03/04/03 */
extern int extra_label;					    /* jwb 12/20/95 */
extern int color_on, color;                                 /* jwb 06/03/96 */
