/******************************************************************************
  		plotsegp.c

        latest edit: 10/15/14
  
  	contains functions plotseg()  plotseg1(), plotseg2(), plotseg3(), 
	plotbar(), plotbar1(), approve(), plabel(), getgoodnums()
 	set_subwindow(), reset_window(), plotseg_log, plotseg2_log,
	set_blackbackground, colorset()
  
  	Programmer: J. Beauchamp  1987  All Rights Reserved.
Changes:
11/09/88 jwb	Install ref to axskip, allowing skip of axes plots.
07/03/89 jwb	Add plotseg2().
07/25/89 jwb	Fix string out in plotbar().
10/11/89 jwb	Change g_inquire_ndc to g_inquire_viewport,
  		 	fix time axis bug, and install plotseg3().
12/17/89 jwb	Remove useless statement in plotbar() (before for).
11/04/90 jwb	Fix no. digits horiz. axis for large xinc.
11/04/90 jwb	Extend range of if(!axskip) and make variables static
  		 	for plotseg() and plotseg1().
11/06/90 jwb	Insert research mode code (approve axis values)
  		including approve() routine.
02/28/91 jwb	Replace HUGE with HUGEI for IEEE huge value given in 
  		macro.h.
02/19/92 jwb	Enlarge bottom Y border from .25 to .32.
11/14/92 jwb	Fix plabel() calls .
12/04/93 jwb	Fix nhdigits in plotseg1().
12/30/93 jwb	plotbar: remove npts from arg list; it's now computed 
  	        from n1 and n2. Improve xinc, nhtick calculation.
02/02/94 jwb	Approve: install user restrictions.
01/04/04 jwb	plotseg2,3: HUGE replaced by HUGEI.
04/27/94 jwb	plabel: move from plotamp.c
10/12/94 jwb	Declare global variables to be static
11/20/95 jwb	Install set_subwindow(), reset_window(), getgoodnums().
		Use set_subwindow() to limit graph to proper graph area
		in plotseg1(), plotseg2(), plotseg3().
11/26/95 jwb	Extra digits for axes with small numbers.
12/20/95 jwb	plotbar(): Eliminate 'string' argument; 
		plabel(): Handle nhars and elabel extra_label output.
05/07/96 jwb	plabel(): Fix extra label output.
		plotbar(): Remove extra elabel output.
  		plotseg2(): Allow user to correct vertical scale under 
  	    	research mode.
05/07/96 tjm	plotseg2(): Introduce color (T. J. Madden). 
		Install colorset() and setblackbackground().
06/03/96 jwb	plotseg(),plotseg1(),plotseg3(): Install color.
06/07/96 tjm	Install plotseg_log(), plotseg2_log().
08/01/96 jwb	Consolidate tjm changes since 6/96
08/07/96 jwb	Fix bug in plotseg2_log().
11/05/96 jwb	Change variable pause to doplot.
11/16/97 zheng	Use gets and sscanf for input.
10/08/98 jwb	Most plot routines: Add right side and top axes (without nos.).
10/09/98 jwb	Research mode: Provide line thickness adjustment.
02/09/03 jwb	plotbar1(): Like plotbar() except that abscissa values
		are floats rather than ints.
		Throughout: Replace npts by nopts to avoid conflict with
		monan.h global variable. 
02/19/03 jwb	Default line thickness constant MINLINETHICK
		plotseg() & plotseg2(): Introduce barplot global variable 
		for bar plot option.
03/04/03 jwb	plabel(): increase size of texout array and have elabel
                put out directly.
06/07/04 jwb	plotseg() and plotseg2(): Fix first array item plot bug.
06/09/04 jwb	plotseg_log() and plotseg2_log(): 
		  Add right side and top axes (without nos.).
06/13/04 jwb	plotbar1(): Fix plotting bugs.
07/05/04 jwb	plotseg(): Use set_subwindow to limit graph to proper graph
		area.
07/12/04 jwb	plotseg(), etc.: Make all default line thicknesses equal to
		MINLINETHICK.
07/13/04 jwb	approve(): Slightly reduce entered plot increment to
		compensate for inaccurate data translation by sscanf.
07/22/04 jwb	Eliminate plabel() so that main program file can have it.
		Replace 'include monan.h' with two other includes.
		Insert externs for compatibility with main program.
		All routines get MINLINETHICK default line thickness.
		All line plot routines set linestyle via global linestyle 
		variable.
04/16/14 jwb	plotseg(): Fix bug that prevented equal vertical values 
		from being plotted. However, plotseg() doesn't seem to work.
		plotbar1() needs research = 1 to override a problem with the
		horizontal increment if 1.0 is not to be used.
10/15/14 jwb	plotseg() didn't work because doplot = 1 was not set from
		the main program. Setting doplot = 1 cured the problem.
		Another problem the equal vertical value fix didn't work
		because the data entry was using %f instead of %lf for a
		double. After making this change it works ok for equal values.
10/15/14 jwb	plotbar1(): There were two problems with plotbar1. First,
		the horizontal increment was assumed to 1.0, and that is 
		the only value that is needed. The code was changed so that
		the increment was based on xmin and xmax, using getgoodnums().
		Second it couldn't handle the case where ymin = ymax. Code
		to ask the user the preferred vertical range was installed.
10/15/14 jwb	Introduce MAXLINETHICK and make all default line thicknesses
		equal to that.
11/25/17 jwb	Declare all functions to be void.
 *****************************************************************************/
#include <math.h>
#include <stdio.h>
#include "g_raph.h"					    /* jwb 07/22/04 */
#include "macro.h"					    /* jwb 07/22/04 */
#define NHTICKS 11
#define NVTICKS 11
#define MINLINETHICK 5					    /* jwb 02/19/03 */
#define MAXLINETHICK 10					    /* jwb 10/15/14 */
#define A .15    /* left X border  */
#define B .05    /* right X border */
#define C .32    /* bottom Y border */
#define D .05    /* top Y border   */
#define P printf
extern int axskip,doplot,research;			    /* jwb 11/05/96 */
extern int barplot;					    /* jwb 02/19/03 */
extern int color_on, color, linestyle;			    /* jwb 07/22/04 */
double xpmin, xpmax, ypmin, ypmax;			    /* jwb 11/20/95 */
static double hinc,vinc;				    /* jwb 12/12/94 */
char resp[80];						    /* jwb 10/09/98 */
void set_subwindow(double, double, double, double,
		   double, double, double, double);	    /* jwb 11/20/95 */
void reset_window(double, double, double, double);	    /* jwb 11/20/95 */
void setblackbackground(double,double,double,double);       /* tjm 05/07/96 */
void colorset(int);				            /* tjm 05/07/96 */
void getgoodnums(double*, double*, double*, int*, int*);    /* jwb 11/22/95 */
void approve(double*,double*,double*,int*,char*);	    /* jwb 11/25/17 */ 

/******************************************************************************
			plotseg()
  This function plots nopts of y vs. x data with axis labels and automatic
  scaling.
 *****************************************************************************/
void plotseg(x,y,nopts,xlabel,ylabel)			    /* jwb 11/25/17 */
float x[],y[]; int nopts; char *xlabel,*ylabel;
{
/*  vertical and horizontal scales determined by program  */
  double xwmin, xwmax, ywmin, ywmax;
  double xmin,xmax,xinc,ymin,ymax,yinc,omag;
  int nhdigits,nvdigits,nvticks,nhticks;
  int i,j, linethick;					    /* jwb 10/09/98 */

  if(!axskip)
  {
    /*	 find the min and max values of x and y  for the plot ranges  */
    xmin = HUGEI; xmax = -HUGEI;  ymin = HUGEI; ymax = -HUGEI;
    for(i=0;i<nopts;i++)
    {
	xmin = min(xmin,(double)x[i]);
	xmax = max(xmax,(double)x[i]);
	ymin = min(ymin,(double)y[i]);
	ymax = max(ymax,(double)y[i]);
    }
    if(ymin == ymax)
    {
      P("All vertical values are %f.\n", ymin);
      while(1)
      { 
        P("What vertical range (ymin ymax) do you want? ");
        gets(resp);					  /* zheng 11/06/97 */
        if (sscanf(resp,"%lf%lf",&ymin,&ymax)==2) break;   /* jwb 10/15/14 */
        P("Got ymin=%lf, ymax=%lf\n\n", ymin, ymax);	   /* jwb 10/15/14 */
      }
    }

    /* compute horizontal and vertical axis parameters */

    xinc = (xmax - xmin)/(NHTICKS -1);
    nhdigits = max(.845 - log10(xinc),0.); /* set for .7 threshold  */
    if(nhdigits) xinc = 
       pow(10.,(double)-nhdigits)*floor(xinc*pow(10.,(double)nhdigits)+.4999);
    else
    { omag = floor(log10(xinc));
      xinc = .5*floor(2.*xinc*pow(10.,-omag)+.4999);
      if((xinc - floor(xinc) > .25)&&(xinc > 2.5)) xinc = ceil(xinc);
      xinc *= pow(10.,omag);
      if(xinc - floor(xinc) > .25) nhdigits = 1;  /* handles X.5 case */
    }
    xmin = floor(xmin/xinc)*xinc;
    xmax = ceil(xmax/xinc)*xinc;
    nhticks = 1 + (xmax - xmin)/xinc + .5;
    if((nhdigits == 0)&&((xinc - floor(xinc)) > .01)) nhdigits = 1;
    if(research) 
    {
      approve(&xmin,&xmax,&xinc,&nhticks,"horizontal");
      if(xinc < 1) nhdigits++; 				    /* jwb 11/26/95 */
    }
    yinc = (ymax - ymin)/(NVTICKS -1);
    nvdigits = max(.99 - log10(yinc),0.);
    if(nvdigits) yinc =
       pow(10.,(double)-nvdigits)*floor(yinc*pow(10.,(double)nvdigits)+.4999);
    else
    { omag = floor(log10(yinc));
      yinc = .5*floor(2.*yinc*pow(10.,-omag)+.4999);
      if((yinc - floor(yinc) > .25)&&(yinc > 2.5)) yinc = ceil(yinc);
      yinc *= pow(10.,omag);
      if(yinc - floor(yinc) > .25) nvdigits = 1;  /* handles Y.5 case */
    }
/*   P("before floor/ceil: ymin=%f, ymax=%f, yinc=%f\n",ymin,ymax,yinc); */
    ymin = floor(ymin/yinc)*yinc;
    ymax = ceil(ymax/yinc)*yinc;
    nvticks = 1 + (ymax - ymin)/yinc +.5;
    if(research) 
    {
      approve(&ymin,&ymax,&yinc,&nvticks,"vertical");
    }
/*    P("after floor/ceil: ymin=%f, ymax=%f, nvticks=%d\n",ymin,ymax,nvticks);
    getchar();  */

    /*	set window within current viewport  */

    xwmin = (1.+A)*xmin - A*xmax;  xwmax = (1.+B)*xmax - B*xmin;
    ywmin = (1.+C)*ymin - C*ymax;  ywmax = (1.+D)*ymax - D*ymin;

    if(research)					    /* jwb 10/09/98 */
    {							    /* jwb 10/09/98 */
      P("Give line thickness (in pixels, default %d): ",    /* jwb 07/12/04 */
           10);						    /* jwb 04/16/14 */
      gets(resp);					    /* jwb 10/09/98 */
    }							    /* jwb 10/09/98 */
    g_init();
    g_set_line_thickness(10);			    	    /* jwb 04/16/04 */
    if(research)					    /* jwb 10/09/98 */
    {							    /* jwb 10/09/98 */
      if(sscanf(resp,"%d",&linethick)==1) 		    /* jwb 10/09/98 */
        g_set_line_thickness(linethick);		    /* jwb 10/09/98 */
    }							    /* jwb 10/09/98 */
    g_set_window(xwmin,xwmax,ywmin,ywmax);
    if(color_on)					    /* jwb 06/03/96 */
    {
      setblackbackground(xwmin,xwmax,ywmin,ywmax);	    /* jwb 06/03/96 */
      g_set_gray(1.0);					    /* jwb 06/03/96 */
    }
    g_axis(xmin,ymin,xmax-xmin,xlabel,1,2,0,xmin,xinc,nhticks,nhdigits);
    g_axis(xmin,ymin,ymax-ymin,ylabel,1,2,90,ymin,yinc,nvticks,nvdigits);
    g_axis(xmin,ymax,xmax-xmin,"",1,2,0,xmin,0.,nhticks,0); /* jwb 10/08/98 */
    g_axis(xmax,ymin,ymax-ymin,"",1,2,90,ymin,0.,nvticks,0);/* jwb 10/08/98 */
    g_text_extent(" ",&hinc,&vinc); 			    /* jwb 11/14/92 */
    plabel(.6*xmin+.4*xwmin, ymin - 4.*vinc);   	    /* jwb 11/14/92 */	

   /*  reset viewport and window for graph */
    set_subwindow(xmin,xmax,ymin,ymax,xwmin,xwmax,ywmin,ywmax);/*jwb 07/05/04*/
  }

   /*  draw graph  */

    g_move_abs(x[0],y[0]);				    /* jwb 08/07/96 */
    if(color_on) colorset(color);			    /* jwb 06/03/96 */
    if(!barplot) /* line graph */			    /* jwb 02/19/03 */
    {							    /* jwb 07/22/04 */
      g_set_linestyle(linestyle);			    /* jwb 07/22/04 */
      for(i=1;i<nopts;i++) g_draw_abs(x[i],y[i]);	    /* jwb 08/07/96 */
    }							    /* jwb 07/22/04 */
    else         /* bar graph */
      for(i=0;i<nopts;i++) 				    /* jwb 06/07/04 */
      {							    /* jwb 02/19/03 */
        g_move_abs(x[i],ymin);				    /* jwb 02/19/03 */
        g_draw_rel(0.,y[i]-ymin);			    /* jwb 02/19/03 */
      }							    /* jwb 02/19/03 */
// fprintf(stderr,"\ndiagnostic: doplot= %d\n\n", doplot);  /* jwb 10/15/14 */

    if(doplot) g_term();				    /* jwb 11/05/96 */
/* doplot = 1: complete graph (wait for CR on TEK)
   doplot = 0: go directly to another plot on same view  */

}  /*  end plotseg()  */

/******************************************************************************
			plotseg_log()
  This function plots nopts of y vs. x data with axis labels and automatic
  scaling. plotseg_log() is same as plotseg() except that a horizontal log 
  axis is used. Installed by Tim Madden, -06/07/96 
 *****************************************************************************/
void plotseg_log(x,y,nopts,xlabel,ylabel) 		    /* jwb 11/25/17 */
float x[],y[]; int nopts; char *xlabel,*ylabel;
{
/*  vertical and horizontal scales determined by program  */
  double xwmin, xwmax, ywmin, ywmax;
  double xmin,xmax,xinc,ymin,ymax,yinc,omag;
  int nhdigits,nvdigits,nvticks,nhticks;
  int i,j, linethick;					    /* jwb 10/09/98 */
  double x_plot;					    /* tjm 06/07/96 */
  int approving;					    /* tjm 06/07/96 */
  static char buffer[80], temp[80];			    /* tjm 06/07/96 */

  if(!axskip)
  {
    /*	 find the min and max values of x and y  for the plot ranges  */
    xmin = HUGEI; xmax = -HUGEI;  ymin = HUGEI; ymax = -HUGEI;
    for(i=0;i<nopts;i++)
    {
	xmin = min(xmin,(double)x[i]);
	xmax = max(xmax,(double)x[i]);
	ymin = min(ymin,(double)y[i]);
	ymax = max(ymax,(double)y[i]);
    }
    if(ymin == ymax)
    {P("All vertical values are %f.\n", ymin);
     while(1)
     { P("What vertical range (ymin ymax) do you want? ");
       gets(resp);					  /* zheng 11/06/97 */
       if (sscanf(resp,"%f%f",&ymin,&ymax)==2) break;
     }
    }

    /* compute horizontal and vertical axis parameters */

    xinc = (xmax - xmin)/(NHTICKS -1);
    nhdigits = max(.845 - log10(xinc),0.); /* set for .7 threshold  */
    if(nhdigits) xinc = 
       pow(10.,(double)-nhdigits)*floor(xinc*pow(10.,(double)nhdigits)+.4999);
    else
    { omag = floor(log10(xinc));
      xinc = .5*floor(2.*xinc*pow(10.,-omag)+.4999);
      if((xinc - floor(xinc) > .25)&&(xinc > 2.5)) xinc = ceil(xinc);
      xinc *= pow(10.,omag);
      if(xinc - floor(xinc) > .25) nhdigits = 1;  /* handles X.5 case */
    }
    nhticks = 1 + (xmax - xmin)/xinc + .5;
    if((nhdigits == 0)&&((xinc - floor(xinc)) > .01)) nhdigits = 1;
    if(research)          /* code in this if statement by   -- tjm 06/07/96 */
    {		
      P("Horizontal Scale: min_val. = %f max_val. = %f \n",xmin,xmax);
      approving = 1;
      while(approving)
      {
        P("Approve? (y/n) ");
        gets(buffer);					  /* zheng 11/06/97 */
        if(sscanf(buffer, "%s", temp) == 1) 
     	{
          if(!strcmp(temp,"y"))      approving = 0;
          else if(!strcmp(temp,"n")) 
	  { while(1)
            { P("Enter New Values: ");
              gets(resp);				  /* zheng 11/06/97 */
	      if (sscanf(resp,"%lf%lf",&xmin,&xmax)==2) break;
	    }
            P("Horizontal Scale: min_val. = %f, max_val. = %f \n",xmin,xmax);
	    approving = 0;
	  }
	} 
      }	/* end while(aproving) */
    }/* end if(research) */

    yinc = (ymax - ymin)/(NVTICKS -1);
    nvdigits = max(.99 - log10(yinc),0.);
    if(nvdigits) yinc =
       pow(10.,(double)-nvdigits)*floor(yinc*pow(10.,(double)nvdigits)+.4999);
    else
    { omag = floor(log10(yinc));
      yinc = .5*floor(2.*yinc*pow(10.,-omag)+.4999);
      if((yinc - floor(yinc) > .25)&&(yinc > 2.5)) yinc = ceil(yinc);
      yinc *= pow(10.,omag);
      if(yinc - floor(yinc) > .25) nvdigits = 1;  /* handles Y.5 case */
    }
/*   P("before floor/ceil: ymin=%f, ymax=%f, yinc=%f\n",ymin,ymax,yinc); */
    ymin = floor(ymin/yinc)*yinc;
    ymax = ceil(ymax/yinc)*yinc;
    nvticks = 1 + (ymax - ymin)/yinc +.5;
    if(research) 
    {
      approve(&ymin,&ymax,&yinc,&nvticks,"vertical");
    }
/*    P("after floor/ceil: ymin=%f, ymax=%f, nvticks=%d\n",ymin,ymax,nvticks);
    getchar();  */

    /*	set window within current viewport  */

    xwmin = (1.+A)*xmin - A*xmax;  xwmax = (1.+B)*xmax - B*xmin;
    ywmin = (1.+C)*ymin - C*ymax;  ywmax = (1.+D)*ymax - D*ymin;

    if(research)					    /* jwb 10/09/98 */
    {							    /* jwb 10/09/98 */
      P("Give line thickness (in pixels, default %d): ",    /* jwb 07/12/04 */
           MAXLINETHICK);                                   /* jwb 10/15/14 */
      gets(resp);
    }							    /* jwb 10/09/98 */
    g_init();
    g_set_line_thickness(MAXLINETHICK);			    /* jwb 10/15/14 */
    if(research)					    /* jwb 10/09/98 */
    {							    /* jwb 10/09/98 */
      if(sscanf(resp,"%d",&linethick)==1) 		    /* jwb 10/09/98 */
        g_set_line_thickness(linethick);		    /* jwb 10/09/98 */
    }							    /* jwb 10/09/98 */
    g_set_window(xwmin,xwmax,ywmin,ywmax);
    if(color_on)					    /* jwb 06/03/96 */
    {
      setblackbackground(xwmin,xwmax,ywmin,ywmax);	    /* jwb 06/03/96 */
      g_set_gray(1.0);					    /* jwb 06/03/96 */
    }
   if(research)						    /* tjm 06/07/96 */
     g_axis_log(xmin,ymin,xmax-xmin,xlabel,
		1,2,0, xmin, xmax, 0, nhdigits);	    /* tjm 06/07/96 */
   else							    /* tjm 06/07/96 */
     g_axis_log(xmin,ymin,xmax-xmin,xlabel,                 /* tjm 06/07/96 */
		1,2,0, xmin, xmax, 1, nhdigits);	    /* tjm 06/07/96 */

    g_axis(xmin,ymin,ymax-ymin,ylabel, 1,2,90,ymin, yinc, nvticks, nvdigits);
    g_axis(xmax,ymin,ymax-ymin,"",1,2,90,ymin,0.,nvticks,0);/* jwb 06/09/04 */
    g_text_extent(" ",&hinc,&vinc); 	
    plabel(.6*xmin+.4*xwmin, ymin - 4.*vinc);   
  }

   /*  draw graph  */
  if(research)   /* code in this if-else block by   --        tjm 06/07/96 */
  {
    g_plot_log(xmin,xmax-xmin,xmin,xmax,0,x[0],&x_plot);
    g_move_abs(x_plot,y[0]);				    /* tjm 07/22/96 */

    if(color_on) colorset(color);			    /* jwb 06/03/96 */
    for(i=1;i<nopts;i++)  
    {
      g_plot_log(xmin,xmax-xmin,xmin,xmax,0,x[i],&x_plot);
      g_draw_abs(x_plot,y[i]);
    }
  }
  else
  {
    g_plot_log(xmin,xmax-xmin,xmin,xmax,1,x[0],&x_plot);
    g_move_abs(x_plot,y[0]);			     	    /* tjm 7/22/96*/

    if(color_on) colorset(color);			   
    for(i=1;i<nopts;i++)  
    {
      g_plot_log(xmin, (xmax-xmin), xmin, xmax,1,x[i],&x_plot);
      g_draw_abs(x_plot,y[i]);
    }
  }

  if(doplot) g_term();					    /* jwb 11/05/96 */
/* doplot = 1: complete graph (wait for CR on TEK)
   doplot = 0: go directly to another plot on same view  */

}  /*  end plotseg_log()  */

/******************************************************************************
                        plotseg1()
  This function plots nopts of y vs. x data with axis labels. Max and min values
  of x and y and numbers of horizontal and vertical axis ticks are specified by 
  calling program.
 *****************************************************************************/

void plotseg1(x,y,nopts,xlabel,ylabel,x1,x2,y1,y2,nhticks,nvticks) 
	float x[],y[]; int nopts; char *xlabel,*ylabel;    /* jwb 11/25/17 */
    	float x1,x2,y1,y2; int nhticks,nvticks;
/*   vertical and horizontal scales predetermined.  call this when you
     know what you want. */
{
  double xwmin, xwmax, ywmin, ywmax;
  double xmin,xmax,ymin,ymax,xinc,yinc,ym;
  double omag;					    	    /* jwb 12/04/93 */
  int nhdigits,nvdigits, linethick;			    /* jwb 10/09/98 */
  int i,j;

  if(!axskip)
  {
     xmin=x1; xmax=x2; ymin=y1; ymax=y2;
    /* compute horizontal and vertical axis parameters */

    xinc = (xmax - xmin)/(nhticks -1);
    nhdigits = max(.845 - log10(xinc),0.);
    omag = pow(10.,(double)nhdigits);			/* jwb 12/04/93 */
    if((nhdigits >= 0)&&((xinc*omag - floor(xinc*omag)) > .01)) nhdigits++;

    if(research) 
    {
      approve(&xmin,&xmax,&xinc,&nhticks,"horizontal");
      if(xinc < 1) nhdigits++;				    /* jwb 11/26/95 */
    }
    yinc = (ymax - ymin)/(nvticks -1);
    nvdigits = max(.845 - log10(yinc),0.);
    if(research) 
    {
      approve(&ymin,&ymax,&yinc,&nvticks,"vertical");
    }
    /*	set mapping of our coords to viewport using g_set_window  */
    xwmin = (1.+A)*xmin - A*xmax;  xwmax = (1.+B)*xmax - B*xmin;
    ywmin = (1.+C)*ymin - C*ymax;  ywmax = (1.+D)*ymax - D*ymin;
    if(research)					    /* jwb 10/09/98 */
    {							    /* jwb 10/09/98 */
      P("Give line thickness (in pixels, default %d): ",    /* jwb 07/12/04 */
           MAXLINETHICK);                                   /* jwb 10/15/14 */
      gets(resp);
    }							    /* jwb 10/09/98 */
    g_init();
    g_set_line_thickness(MAXLINETHICK);			    /* jwb 10/15/14 */
    if(research)					    /* jwb 10/09/98 */
    {							    /* jwb 10/09/98 */
      if(sscanf(resp,"%d",&linethick)==1) 		    /* jwb 10/09/98 */
        g_set_line_thickness(linethick);		    /* jwb 10/09/98 */
    }							    /* jwb 10/09/98 */
    g_set_window(xwmin,xwmax,ywmin,ywmax);
    if(color_on)					    /* jwb 06/03/96 */
    {
      setblackbackground(xwmin,xwmax,ywmin,ywmax);	    /* jwb 06/03/96 */
      g_set_gray(1.0);					    /* jwb 06/03/96 */
    }
    g_axis(xmin,ymin,xmax-xmin,xlabel,1,2, 0,xmin,xinc,nhticks,nhdigits);
    g_axis(xmin,ymin,ymax-ymin,ylabel,1,2,90,ymin,yinc,nvticks,nvdigits);
    g_axis(xmin,ymax,xmax-xmin,"",1,2, 0,xmin,0.,nhticks,0); /* jwb 10/08/98 */
    g_axis(xmax,ymin,ymax-ymin,"",1,2,90,ymin,0.,nvticks,0); /* jwb 10/08/98 */
    g_text_extent(" ",&hinc,&vinc); 			    /* jwb 11/14/92 */
    plabel(.6*xmin+.4*xwmin, ymin - 4.*vinc);   	    /* jwb 11/14/92 */	

   /*  reset viewport and window for graph */
    set_subwindow(xmin,xmax,ymin,ymax,xwmin,xwmax,ywmin,ywmax);/*jwb 11/20/95*/
  }

/*	draw graph  */
    if(color_on) colorset(color);			    /* jwb 06/03/96 */
    g_move_abs(x[0],y[0]);				    /* jwb 07/22/04 */
    g_set_linestyle(linestyle);				    /* jwb 07/22/04 */
    for(i=1;i<nopts;i++)  g_draw_abs(x[i],y[i]);	    /* jwb 07/22/04 */

    if(doplot) g_term();				    /* jwb 11/05/96 */
/* doplot = 1: complete graph (wait for CR on TEK)
   doplot = 0: go directly to another plot on same view  */

} /* end plotseg1() */

/******************************************************************************
                        plotseg2()
  This function plots nopts of y vs. x data with axis labels. Automatic
  scaling is done for the horizontal axis, but min and max values of y and the
  number of vertical axis ticks are given by the calling program.
 *****************************************************************************/
void plotseg2(x,y,nopts,xlabel,ylabel,y1,y2,nvticks)	    /* jwb 11/25/17 */ 
	float x[],y[]; int nopts; char *xlabel,*ylabel;
    	float y1,y2; int nvticks;
/*  same as plotseg() except vertical scale is predetermined   */
{
    double xwmin, xwmax, ywmin, ywmax;
    double xmin,xmax,ymin,ymax,xinc,yinc,omag;
    int nhdigits,nvdigits,nhticks;
    int i,j, linethick;					    /* jwb 10/09/98 */

  if(!axskip)
  {
     ymin=y1; ymax=y2;
    /*	 find the min and max values of x for the plot range  */
    xmin = HUGEI; xmax = -HUGEI;
    for(i=0;i<nopts;i++)
    {
	xmin = min(xmin,(double)x[i]);
	xmax = max(xmax,(double)x[i]);
    }
    if(ymin == ymax)
    { while(1)
      { P("give vertical range (ymin ymax): ");
        gets(resp);					/* zheng 11/06/97 */
        if (sscanf(resp,"%f%f",&ymin,&ymax)==2) break;
      }
    }

    /* compute horizontal axis parameters */
    xinc = (xmax - xmin)/(NHTICKS -1);
    nhdigits = max(.845 - log10(xinc),0.); /* set for .7 threshold  */
    if(nhdigits) xinc = 
       pow(10.,(double)-nhdigits)*floor(xinc*pow(10.,(double)nhdigits)+.4999);
    else
    { omag = floor(log10(xinc));
      xinc = .5*floor(2.*xinc*pow(10.,-omag)+.4999);
      if((xinc - floor(xinc) > .25)&&(xinc > 2.5)) xinc = ceil(xinc);
      xinc *= pow(10.,omag);
      if(xinc - floor(xinc) > .25) nhdigits = 1;  /* handles X.5 case */
    }
    xmin = floor(xmin/xinc)*xinc;
    xmax = ceil(xmax/xinc)*xinc;
    nhticks = 1 + (xmax - xmin)/xinc + .5;
    if((nhdigits == 0)&&((xinc - floor(xinc)) > .01)) nhdigits = 1;
    if(research) 
    {
      approve(&xmin,&xmax,&xinc,&nhticks,"horizontal");
      if(xinc < 1) nhdigits++;				     /* jwb 11/26/95 */
    }
    /* compute vertical axis parameters */
    yinc = (ymax - ymin)/(nvticks -1);
    nvdigits = max(.99 - log10(yinc),0.);
    if(research)					     /* jwb 05/07/96 */
    {
      approve(&ymin,&ymax,&yinc,&nvticks,"vertical");        /* jwb 05/07/96 */
    }
    /*	set mapping of our coords to viewport using g_set_window  */
    xwmin = (1.+A)*xmin - A*xmax;  xwmax = (1.+B)*xmax - B*xmin;
    ywmin = (1.+C)*ymin - C*ymax;  ywmax = (1.+D)*ymax - D*ymin;
    if(research)					    /* jwb 10/09/98 */
    {							    /* jwb 10/09/98 */
      P("Give line thickness (in pixels, default %d): ",    /* jwb 07/12/04 */
           MAXLINETHICK);                                   /* jwb 10/15/14 */
      gets(resp);
    }							    /* jwb 10/09/98 */
/*
   P("ymin=%f, yinc=%f, nvticks=%d, nvdigits=%d\n",
          ymin, yinc, nvticks, nvdigits);   
*/
    g_init();
    g_set_line_thickness(MAXLINETHICK);			    /* jwb 10/15/14 */
    if(research)					    /* jwb 10/09/98 */
    {							    /* jwb 10/09/98 */
      if(sscanf(resp,"%d",&linethick)==1) 		    /* jwb 10/09/98 */
        g_set_line_thickness(linethick);		    /* jwb 10/09/98 */
    }							    /* jwb 10/09/98 */
    g_set_window(xwmin,xwmax,ywmin,ywmax);
    if(color_on)					    /* tjm 05/07/96 */
    {
      setblackbackground(xwmin,xwmax,ywmin,ywmax);	    /* tjm 05/07/96 */
      g_set_gray(1.0);					    /* tjm 05/07/96 */
    }
    g_axis(xmin,ymin,xmax-xmin,xlabel,1,2, 0,xmin,xinc,nhticks,nhdigits);
    g_axis(xmin,ymin,ymax-ymin,ylabel,1,2,90,ymin,yinc,nvticks,nvdigits);
    g_axis(xmin,ymax,xmax-xmin,"",1,2, 0,xmin,0.,nhticks,0);/* jwb 10/08/98 */
    g_axis(xmax,ymin,ymax-ymin,"",1,2,90,ymin,0.,nvticks,0);/* jwb 10/08/98 */
    g_text_extent(" ",&hinc,&vinc); 			    /* jwb 11/14/92 */
    plabel(.6*xmin+.4*xwmin, ymin - 4.*vinc);   	    /* jwb 11/14/92 */	

   /*  reset viewport and window and draw graph  */
    set_subwindow(xmin,xmax,ymin,ymax,xwmin,xwmax,ywmin,ywmax);/*jwb 11/20/95*/

  } /* end !axskip */

    if(color_on) colorset(color);			    /* tjm 05/07/96 */
      
    g_move_abs(x[0],y[0]);				    /* jwb 07/22/04 */
    if(!barplot)					    /* jwb 02/19/03 */
    {							    /* jwb 07/22/04 */
      g_set_linestyle(linestyle);			    /* jwb 07/22/04 */
      for(i=1;i<nopts;i++) g_draw_abs(x[i],y[i]);
    }							    /* jwb 07/22/04 */
    else						    /* jwb 02/19/03 */
      for(i=0;i<nopts;i++) 				    /* jwb 06/07/04 */
      {							    /* jwb 02/19/03 */
        g_move_abs(x[i],ymin);				    /* jwb 02/19/03 */
        g_draw_rel(0.,y[i]-ymin);			    /* jwb 02/19/03 */
      }							    /* jwb 02/19/03 */

    if(doplot) g_term();				    /* jwb 11/05/96 */
/* doplot = 1: complete graph (wait for CR on TEK)
   doplot = 0: go directly to another plot on same view  */

}  /* end plotseg2()  */

/******************************************************************************
                        plotseg2_log()
  This function plots nopts of y vs. x data with axis labels and automatic
  scaling of x data. Scaling for y data is specified by calling program.
  plotseg2_log() is same as plotseg2() except that a horizontal log axis is 
  used. It is the same as plotseg_log() except that the vertical scale is 
  predetermined.
  Installed by Tim Madden, 06/07/96
 *****************************************************************************/
void plotseg2_log(x,y,nopts,xlabel,ylabel,y1,y2,nvticks)    /* jwb 11/25/17 */ 
	float x[],y[]; int nopts; char *xlabel,*ylabel;
    	float y1,y2; int nvticks;
{
    double xwmin, xwmax, ywmin, ywmax;
    double xmin,xmax,ymin,ymax,xinc,yinc,omag;
    int nhdigits,nvdigits,nhticks;
    int i,j;
  double x_plot;					    /* tjm 06/07/96 */
  int approving;					    /* tjm 06/07/96 */
  static char buffer[80], temp[80];			    /* tjm 06/07/96 */
  int linethick;					    /* jwb 10/09/98 */

  if(!axskip)
  {
     ymin=y1; ymax=y2;
    /*	 find the min and max values of x for the plot range  */
    xmin = HUGEI; xmax = -HUGEI;
    for(i=0;i<nopts;i++)
    {
	xmin = min(xmin,(double)x[i]);
	xmax = max(xmax,(double)x[i]);
    }
    if(ymin == ymax)
    {while(1)
     { P("give vertical range (ymin ymax): ");
       gets(resp);					/* zheng 11/06/97 */
       if (sscanf(resp,"%f%f",&ymin,&ymax)==2) break;
     }
    }

    /* compute horizontal axis parameters */
    xinc = (xmax - xmin)/(NHTICKS -1);
    nhdigits = max(.845 - log10(xinc),0.); /* set for .7 threshold  */
    if(nhdigits) xinc = 
       pow(10.,(double)-nhdigits)*floor(xinc*pow(10.,(double)nhdigits)+.4999);
    else
    { omag = floor(log10(xinc));
      xinc = .5*floor(2.*xinc*pow(10.,-omag)+.4999);
      if((xinc - floor(xinc) > .25)&&(xinc > 2.5)) xinc = ceil(xinc);
      xinc *= pow(10.,omag);
      if(xinc - floor(xinc) > .25) nhdigits = 1;  /* handles X.5 case */
    }
    nhticks = 1 + (xmax - xmin)/xinc + .5;
    if((nhdigits == 0)&&((xinc - floor(xinc)) > .01)) nhdigits = 1;
    if(research) /*  code within this if statement by   --     tjm 06/07/96 */
    {
      P("Horizontal Scale: min_val. = %f max_val. = %f \n",xmin,xmax);
      approving = 1;
      while(approving)
      {
        P("Approve? (y/n) ");
        gets(buffer);					/* zheng 11/06/97 */
        if(sscanf(buffer, "%s", temp) == 1) 
     	{
          if(!strcmp(temp,"y"))      approving = 0;
          else if(!strcmp(temp,"n")) 
          { while(1)
            { P("Enter New Values: ");
              gets(resp);				/* zheng 11/06/97 */
	      if (sscanf(resp,"%lf%lf",&xmin,&xmax)==2) break;
	    }
            P("Horizontal Scale: min_val. = %f max_val. = %f \n",xmin,xmax);
            approving = 0;
          }
	}
     }
    }/* end if(research) */

    /* compute vertical axis parameters */
    yinc = (ymax - ymin)/(nvticks -1);
    nvdigits = max(.99 - log10(yinc),0.);
    if(research)
    {
      approve(&ymin,&ymax,&yinc,&nvticks,"vertical");
    }
    /*	set mapping of our coords to viewport using g_set_window  */
    xwmin = (1.+A)*xmin - A*xmax;  xwmax = (1.+B)*xmax - B*xmin;
    ywmin = (1.+C)*ymin - C*ymax;  ywmax = (1.+D)*ymax - D*ymin;
    if(research)					    /* jwb 10/09/98 */
    {							    /* jwb 10/09/98 */
      P("Give line thickness (in pixels, default %d): ",    /* jwb 07/12/04 */
           MAXLINETHICK);                                   /* jwb 10/15/14 */
      gets(resp);
    }							    /* jwb 10/09/98 */
    g_init();
    g_set_line_thickness(MAXLINETHICK);			    /* jwb 10/15/14 */
    if(research)					    /* jwb 10/09/98 */
    {							    /* jwb 10/09/98 */
      if(sscanf(resp,"%d",&linethick)==1) 		    /* jwb 10/09/98 */
        g_set_line_thickness(linethick);		    /* jwb 10/09/98 */
    }							    /* jwb 10/09/98 */
    g_set_window(xwmin,xwmax,ywmin,ywmax);
    if(color_on)
    {
      setblackbackground(xwmin,xwmax,ywmin,ywmax);
      g_set_gray(1.0);	
    }
    if(research) /* code in this if-else block by    --        tjm 06/07/96 */
    {							    /* jwb 06/09/04 */
      g_axis_log(xmin,ymin,xmax-xmin,xlabel,1,2,0,xmin,xmax,0,nhdigits);
      g_axis_log(xmin,ymax,xmax-xmin,"",1,2,0,xmin,xmax,0,nhdigits); /*6/9/04*/
    }							    /* jwb 06/09/04 */
    else
    {							    /* jwb 06/09/04 */
      g_axis_log(xmin,ymin,xmax-xmin,xlabel,1,2,0,xmin,xmax,1,nhdigits);
      g_axis_log(xmin,ymax,xmax-xmin,"",1,2,0,xmin,xmax,1,nhdigits); /*6/9/04*/
    }							    /* jwb 06/09/04 */
    g_axis(xmin,ymin,ymax-ymin,ylabel, 1,2,90,ymin, yinc, nvticks, nvdigits);
    g_axis(xmax,ymin,ymax-ymin,"",1,2,90,ymin,0.,nvticks,0);/* jwb 06/09/04 */
    g_text_extent(" ",&hinc,&vinc);
    plabel(.6*xmin+.4*xwmin, ymin - 4.*vinc);

   /*  reset viewport and window and draw graph  */
    set_subwindow(xmin,xmax,ymin,ymax,xwmin,xwmax,ywmin,ywmax);

  } /* end !axskip */

      
  if(research) /* code in this if-else block by       --       tjm 06/07/96 */
  {
    g_plot_log(xmin,xmax-xmin,xmin,xmax,0,x[0],&x_plot);
    g_move_abs(x_plot,y[0]);

    if(color_on) colorset(color);
    for(i=1;i<nopts;i++)  
    {
      g_plot_log(xmin,xmax-xmin,xmin,xmax,0,x[i],&x_plot);
      g_draw_abs(x_plot,y[i]);
    }
  }
  else
  {
    g_plot_log(xmin,xmax-xmin,xmin,xmax,1,x[0],&x_plot);    /* jwb 08/07/96 */
    g_move_abs(x_plot,y[0]);			    	    /* jwb 08/07/96 */

    if(color_on) colorset(color);
    for(i=1;i<nopts;i++)  
    {
      g_plot_log(xmin,xmax-xmin,xmin,xmax,1,x[i],&x_plot);
      g_draw_abs(x_plot,y[i]);
    }
  }  /* if(research)-else block */

  if(doplot) g_term();					    /* jwb 11/05/96 */
/* doplot = 1: complete graph (wait for CR on TEK)
   doplot = 0: go directly to another plot on same view  */

}  /* end plotseg2_log()  */

/******************************************************************************
                        plotseg3()
  This function plots nopts of y vs. x data with axis labels. Automatic
  scaling is done for the horizontal axis, but min and max values of x and the
  number of horizonatal axis ticks are given by the calling program.
 *****************************************************************************/
void plotseg3(x,y,nopts,xlabel,ylabel,x1,x2,nhticks)	    /* jwb 11/25/17 */ 
	float x[],y[]; int nopts; char *xlabel,*ylabel;
    	float x1,x2; int nhticks;
{
    double xwmin, xwmax, ywmin, ywmax;
    static double xmin,xmax,ymin,ymax,xinc,yinc,omag;
    static int nhdigits,nvdigits,nvticks;
    int i,j, linethick;					    /* jwb 10/09/98 */

  if(!axskip)
  {
     xmin = x1; xmax = x2;
    /*	 find the min and max values of y for the plot range  */
    ymin = HUGEI; ymax = -HUGEI;
    for(i=0;i<nopts;i++)
    {
	ymin = min(ymin,(double)y[i]);
	ymax = max(ymax,(double)y[i]);
    }

    /* compute horizontal axis parameters */
    xinc = (xmax - xmin)/(nhticks -1);
    nhdigits = max(.845 - log10(xinc),0.);
    omag = pow(10.,(double)nhdigits);
    if((nhdigits >= 0)&&((xinc*omag - floor(xinc*omag)) > .01)) nhdigits++;
    if(research) 
    {
      approve(&xmin,&xmax,&xinc,&nhticks,"horizontal");
      if(xinc < 1) nhdigits++;				    /* jwb 11/26/95 */
    }
    /* compute vertical axis parameters */
    yinc = (ymax - ymin)/(NVTICKS -1);
    nvdigits = max(.845 - log10(yinc),0.); /* set for .7 threshold  */
    if(nvdigits) yinc = 
       pow(10.,(double)-nvdigits)*floor(yinc*pow(10.,(double)nvdigits)+.4999);
    else
    { omag = floor(log10(yinc));
      yinc = .5*floor(2.*yinc*pow(10.,-omag)+.4999);
      if((yinc - floor(yinc) > .25)&&(yinc > 2.5)) yinc = ceil(yinc);
      yinc *= pow(10.,omag);
      if(yinc - floor(yinc) > .25) nvdigits = 1;  /* handles X.5 case */
    }
    ymin = floor(ymin/yinc)*yinc;
    ymax = ceil(ymax/yinc)*yinc;
    nvticks = 1 + (ymax - ymin)/yinc + .5;
    if(research) 
    {
      approve(&ymin,&ymax,&yinc,&nvticks,"vertical");
    }
    /*	set mapping of our coords to viewport using g_set_window  */
    xwmin = (1.+A)*xmin - A*xmax;  xwmax = (1.+B)*xmax - B*xmin;
    ywmin = (1.+C)*ymin - C*ymax;  ywmax = (1.+D)*ymax - D*ymin;
    if(research)					    /* jwb 10/09/98 */
    {							    /* jwb 10/09/98 */
      P("Give line thickness (in pixels, default %d): ",    /* jwb 07/12/04 */
           MAXLINETHICK);                                   /* jwb 10/15/14 */
      gets(resp);
    }							    /* jwb 10/09/98 */
    g_init();
    g_set_line_thickness(MAXLINETHICK);			    /* jwb 10/15/14 */
    if(research)					    /* jwb 10/09/98 */
    {							    /* jwb 10/09/98 */
      if(sscanf(resp,"%d",&linethick)==1) 		    /* jwb 10/09/98 */
        g_set_line_thickness(linethick);		    /* jwb 10/09/98 */
    }							    /* jwb 10/09/98 */
    g_set_window(xwmin,xwmax,ywmin,ywmax);
    if(color_on)					    /* jwb 06/03/96 */
    {
      setblackbackground(xwmin,xwmax,ywmin,ywmax);	    /* jwb 06/03/96 */
      g_set_gray(1.0);					    /* jwb 06/03/96 */
    }
    g_axis(xmin,ymin,xmax-xmin,xlabel,1,2, 0,xmin,xinc,nhticks,nhdigits);
    g_axis(xmin,ymin,ymax-ymin,ylabel,1,2,90,ymin,yinc,nvticks,nvdigits);
    g_axis(xmin,ymax,xmax-xmin,"",1,2, 0,xmin,0.,nhticks,0); /* jwb 10/08/98 */
    g_axis(xmax,ymin,ymax-ymin,"",1,2,90,ymin,0.,nvticks,0); /* jwb 10/08/98 */
    g_text_extent(" ",&hinc,&vinc); 			     /* jwb 11/14/92 */
    plabel(.6*xmin+.4*xwmin, ymin - 4.*vinc);   	     /* jwb 11/14/92 */	

   /*  reset viewport and window and draw graph  */
    set_subwindow(xmin,xmax,ymin,ymax,xwmin,xwmax,ywmin,ywmax);/*jwb 11/20/95*/
  } /* end !axskip */

    if(color_on) colorset(color);			    /* jwb 06/03/96 */

    g_move_abs(x[0],y[0]);				    /* jwb 07/22/04 */
    g_set_linestyle(linestyle);				    /* jwb 07/22/04 */
    for(i=1;i<nopts;i++)  g_draw_abs(x[i],y[i]);	    /* jwb 07/22/04 */

    if(doplot) g_term();				    /* jwb 11/05/96 */
/* doplot = 1: wait for CR to jetison plot 
   doplot = 0: go directly to another plot on same view  */

}  /* end plotseg3() */

/*****************************************************************************
				plotbar()
   This function plots array y[] from index n1 to n2 as a bar plot.
   It includes automatic scaling and axis labeling.
 ****************************************************************************/
void plotbar(n1,n2,y,xlabel,ylabel) 			    /* jwb 11/25/17 */
float y[]; int n1, n2; char *xlabel,*ylabel; 		    /* jwb 12/20/95 */
{
    double xwmin,xwmax,ywmin,ywmax;
    double xmin,xmax,ymin,ymax;
    double xinc,yinc,dx,ym,omag;
    int i,j,nhdigits,nvdigits,nhticks,nvticks,nopts;  	    /* jwb 12/30/93 */
    int linethick;					    /* jwb 10/09/98 */

    /*	 find the min and max values of x and y  for the plot ranges  */
    xmin = n1; xmax = n2;   ymin = ymax = 0.;
    nopts = n2 -n1 + 1;					    /* jwb 12/30/93 */
    for(i=0;i<nopts;i++)
    {
	ymin = min(ymin,(double)y[i]);
	if((double)y[i] > ymax) ymax = y[i]; 
    }
    ym = ymax;

    /* compute horizontal and vertical axis parameters */
    xinc = 1.;
    while(1)
    {
      nhticks = 1 + ceil((nopts-1)/xinc);		/* jwb 12/30/93 */
      if(nhticks <= 13) break;				/* jwb 12/30/93 */
      xinc += 1.;
    }
    xmax = xmin + (nhticks - 1)*xinc;
    nhdigits = 0;
    if(research) 
    {
      approve(&xmin,&xmax,&xinc,&nhticks,"horizontal");
      if(xinc < 1) nhdigits++;				    /* jwb 11/26/95 */
    }
    yinc = (ymax - ymin)/(NVTICKS -1);
    nvdigits = max(.99 - log10(yinc),0.);
    if(nvdigits) yinc =
       pow(10.,(double)-nvdigits)*floor(yinc*pow(10.,(double)nvdigits)+.4999);
    else
    { omag = floor(log10(yinc));
      yinc = .5*floor(2.*yinc*pow(10.,-omag)+.4999);
      if((yinc - floor(yinc) > .25)&&(yinc > 2.5)) yinc = ceil(yinc);
      yinc *= pow(10.,omag);
      if(yinc - floor(yinc) > .25) nvdigits = 1;  /* handles Y.5 case */
    }
/*   P("before floor/ceil: ymin=%f, ymax=%f, yinc=%f\n",ymin,ymax,yinc); */
    ymin = floor(ymin/yinc)*yinc;
    ymax = ceil(ymax/yinc)*yinc;
    nvticks = 1 + (ymax - ymin)/yinc +.5;
    if(research) 
    {
      approve(&ymin,&ymax,&yinc,&nvticks,"vertical");
    }
    /*	set window within current viewport  */
    xwmin = (1.+A)*xmin - A*xmax;  xwmax = (1.+B)*xmax - B*xmin;
    ywmin = (1.+C)*ymin - C*ymax;  ywmax = (1.+D)*ymax - D*ymin;
/*  printf("nvticks=%d,y[0]=%f,yinc=%f\n",nvticks,y[0],yinc);
    printf("xmin=%f,xmax=%f,ymin=%f,ymax=%f\n",xmin,xmax,ymin,ymax);
   printf("xwmin=%f,xwmax=%f,ywmin=%f,ywmax=%f=\n",xwmin,xwmax,ywmin,ywmax);
    getchar();
*/
    if(research)					    /* jwb 10/09/98 */
    {							    /* jwb 10/09/98 */
      P("Give line thickness (in pixels, default %d): ",    /* jwb 07/12/04 */
           MAXLINETHICK);                                   /* jwb 10/15/14 */
      gets(resp);
    }							    /* jwb 10/09/98 */
    g_init();
    g_set_line_thickness(MAXLINETHICK);			    /* jwb 10/15/14 */
    if(research)					    /* jwb 10/09/98 */
    {							    /* jwb 10/09/98 */
      if(sscanf(resp,"%d",&linethick)==1) 		    /* jwb 10/09/98 */
        g_set_line_thickness(linethick);		    /* jwb 10/09/98 */
    }							    /* jwb 10/09/98 */
    g_set_window(xwmin,xwmax,ywmin,ywmax);
    g_axis(xmin,ymin,xmax-xmin,xlabel,1,2, 0,xmin,xinc,nhticks,nhdigits);
    g_axis(xmin,ymin,ymax-ymin,ylabel,1,2,90,ymin,yinc,nvticks,nvdigits);
    g_axis(xmin,ymax,xmax-xmin,"",1,2, 0,xmin,0.,nhticks,0); /* jwb 10/08/98 */
    g_axis(xmax,ymin,ymax-ymin,"",1,2,90,ymin,0.,nvticks,0); /* jwb 10/08/98 */
    g_text_extent(" ",&hinc,&vinc); 
    plabel(.6*xmin+.4*xwmin, ymin - 4.*vinc);   	     /* jwb 04/06/92 */

   /*  draw graph  */

    for(i=1;i<nopts;i++)
    {
       g_move_abs(xmin+i,0.);
       g_draw_rel(0.,y[i]-ymin);			    /* jwb 02/19/03 */
    }
    g_term();
} /* end plotbar() */

/*****************************************************************************
				plotbar1()
   This function plots array y[] vs. x[] as a bar plot.
   It includes automatic scaling and axis labeling.
 ****************************************************************************/
void plotbar1(x,y,nopts,xlabel,ylabel)			    /* jwb 11/25/17 */
float x[], y[]; 					    /* jwb 02/09/03 */
int nopts;						    /* jwb 02/09/03 */
char *xlabel,*ylabel; 		    	    		    /* jwb 12/20/95 */
{
    double xwmin,xwmax,ywmin,ywmax;
    double xmin,xmax,ymin,ymax;
    double xinc,yinc,omag;
    int i,j,nhdigits,nvdigits,nhticks,nvticks;  	    /* jwb 02/09/03 */
    int linethick;					    /* jwb 10/09/98 */

    /*	 find the min and max values of x and y  for the plot ranges  */
    xmin = xmax = x[0];					    /* jwb 10/15/14 */
    for(i=0;i<nopts;i++)				    /* jwb 10/15/14 */
    {							    /* jwb 10/15/14 */
	xmin = min(xmin,(double)x[i]);			    /* jwb 10/15/14 */
	xmax = max(xmax,(double)x[i]);			    /* jwb 10/15/14 */
    }							    /* jwb 10/15/14 */
// P("before goodnums: xmin = %.1lf, xmax = %.1lf\n", xmin, xmax);
    getgoodnums(&xmin, &xmax, &xinc, &nhticks, &nhdigits);  /* jwb 10/15/14 */
// P("after goodnums: xmin = %.1lf, xmax = %.1lf ", xmin, xmax);
// P("xinc = %.1lf, nhticks = %d, nhdigits = %d\n\n", xinc, nhticks, nhdigits);

    ymin = ymax = y[0];					    /* jwb 10/15/14 */
    for(i=0;i<nopts;i++)
    {
	ymin = min(ymin,(double)y[i]);
	ymax = max(ymax,(double)y[i]);			    /* jwb 10/15/14 */
    }
    if(ymin==ymax)					    /* jwb 10/15/14 */
    {							    /* jwb 10/15/14 */
      P("ymin = ymax = %lf\n", ymin);			    /* jwb 10/15/14 */
      P("Give vertical axis range (ymin ymax): ");	    /* jwb 10/15/14 */
      scanf("%lf%lf%*c", &ymin, &ymax);			    /* jwb 10/15/14 */
      P("Got new ymin = %lf, ymax = %lf\n\n", ymin, ymax);  /* jwb 10/15/14 */
    }							    /* jwb 10/15/14 */
// P("before goodnums: ymin = %.1lf, ymax = %.1lf\n", ymin, ymax);
    getgoodnums(&ymin, &ymax, &yinc, &nvticks, &nvdigits);  /* jwb 10/15/14 */
// P("after goodnums: ymin = %.1lf, ymax = %.1lf ", ymin, ymax);
// P("yinc = %.1lf, nvticks = %d, nvdigits = %d\n\n", yinc, nvticks, nvdigits);

    /* approve and possible alter  horizontal and vertical axis parameters */
    if(research) 
    {
      approve(&xmin,&xmax,&xinc,&nhticks,"horizontal");
      if(xinc < 1.) nhdigits++;				    /* jwb 04/16/14 */
    }
    if(research) 
    {
      approve(&ymin,&ymax,&yinc,&nvticks,"vertical");
      if(yinc < 1) nvdigits++;				  
    }
    /*	set window within current viewport  */
    xwmin = (1.+A)*xmin - A*xmax;  xwmax = (1.+B)*xmax - B*xmin;
    ywmin = (1.+C)*ymin - C*ymax;  ywmax = (1.+D)*ymax - D*ymin;
/*
    P("nvticks=%d,y[0]=%f,yinc=%f\n",nvticks,y[0],yinc);
    P("xwmin=%f,xwmax=%f,ywmin=%f,ywmax=%f=\n",xwmin,xwmax,ywmin,ywmax);
    P("xmin=%f,xmax=%f,ymin=%f,ymax=%f\n",xmin,xmax,ymin,ymax);
    for(i=0;i<nopts;i++)
      P("x[%d]=%.3f, y[%d]=%.3f\n", i, x[i], i, y[i]);
    getchar();
*/
    if(research)					   
    {							    
      P("Give line thickness (in pixels, default %d): ", 
				MAXLINETHICK); 		    /* jwb 10/15/14 */
      gets(resp);					   
    }							    
    g_init();
    g_set_line_thickness(MAXLINETHICK);			    /* jwb 10/15/14 */
    if(research)					    
    {					
      if(sscanf(resp,"%d",&linethick)==1) 	
        g_set_line_thickness(linethick);		
    }							   
    g_set_window(xwmin,xwmax,ywmin,ywmax);
    g_axis(xmin,ymin,xmax-xmin,xlabel,1,2, 0,xmin,xinc,nhticks,nhdigits);
    g_axis(xmin,ymin,ymax-ymin,ylabel,1,2,90,ymin,yinc,nvticks,nvdigits);
    g_axis(xmin,ymax,xmax-xmin,"",1,2, 0,xmin,0.,nhticks,0);
    g_axis(xmax,ymin,ymax-ymin,"",1,2,90,ymin,0.,nvticks,0);
    g_text_extent(" ",&hinc,&vinc); 
    plabel(.6*xmin+.4*xwmin, ymin - 4.*vinc);   	    

   /*  draw graph  */

    g_move_abs(x[0],y[0]);				    /* jwb 06/13/04 */
    for(i=0;i<nopts;i++)				    /* jwb 06/13/04 */
    {
       g_move_abs(x[i],ymin);				    /* jwb 06/13/04 */
       g_draw_rel(0.,max(y[i]-ymin,0.));		    /* jwb 06/13/04 */
    }
    g_term();
} /* end plotbar1() */

void approve(xmin,xmax,inc,ticks,type)			    /* jwb 11/25/17 */ 
  double *xmin,*xmax,*inc;int *ticks;char *type;
{
  char resp[10],choice;
  P("%s min. val. = %f, max. val. = %f, incr. = %f, no. ticks = %d\n", 
				type,*xmin,*xmax,*inc,*ticks);
  while(1)
  { P("approve? (y/n) ");
    gets(resp);						   /* zheng 11/06/97 */
    if (sscanf(resp,"%c",&choice)==1) break;
  }
  if(choice == 'n') 
  {
    P("Give replacement values\n");
    while(1)
    { P("min. val: ");
      gets(resp);					   /* zheng 11/06/97 */
      if (sscanf(resp,"%lf",xmin)==1) break;
    }
incr:							    /* jwb 01/02/94 */
    while(1)
    { P("increment: ");
      gets(resp);					   /* zheng 11/06/97 */
      if (sscanf(resp,"%lf",inc)==1) break;
    }
    if(*inc <= 0) {P("increment must be positive!\n"); goto incr;} /* 1/2/94 */
/* correct slightly inaccurate data translatation of inc by sscanf */
    *inc *= .9999;					    /* jwb 07/13/04 */
noticks:						    /* jwb 01/02/94 */
    while(1)
    { P("no. ticks: ");
      gets(resp);					   /* zheng 11/06/97 */
      if (sscanf(resp,"%d",ticks)==1) break;
    }
    if(*ticks < 2) {P("Too few ticks!\n"); goto noticks;}    /* jwb 01/02/94 */
    if(*ticks > 21) {P("Too many ticks!\n"); goto noticks;}  /* jwb 01/02/94 */
    *xmax = *xmin + (*ticks - 1)*(*inc);
  } /* end do not approve choice */
} /* end approve() */

/*****************************************************************************
			getgoodnums()

	function to compute "good value" of x large than given value, no. of
	ticks (nt), and increment (dx).
 ****************************************************************************/
void getgoodnums(double *xmin,double *xmax,double *xinc,int *ntk,int *ndig)
{
  double dx, omag;
  int nd;
  dx = (*xmax - *xmin)/(NHTICKS-1);
  nd = max(.845 - log10(dx),0.);
  if(nd) dx = pow(10.,(double)-nd)*floor(dx*pow(10.,(double)nd) + .4999);
  else
  {
    omag = floor(log10(dx));
    dx = 0.5*floor(2.*dx*pow(10.,-omag) + .4999);
    if((dx - floor(dx) > .25)&&(dx > 2.5)) dx = ceil(dx);
    dx *= pow(10.,omag);
    if(dx - floor(dx) > .25) nd = 1;
  }
  *xmin = floor(*xmin/dx)*dx;
  *xmax = ceil(*xmax/dx)*dx;
  *xinc = dx;
  *ntk = 1 + (*xmax - *xmin)/dx + .5;
  if((nd == 0)&&((dx - floor(dx)) > .01)) nd = 1;
  *ndig = nd;
}

/******************************************************************************
  Make a subwindow inside of the full viewport window.
  Assumes that the full viewport has been set and not changed.
  However this routine changes the viewport size, so if this routine is
  called again, the viewport probably will need to be restored using
  reset_window().
 *****************************************************************************/
void set_subwindow(double x1,double x2,double y1,double y2,
  		double xw1,double xw2,double yw1,double yw2)
{
  double xpmin1, xpmax1, ypmin1, ypmax1, xratio, yratio;

  g_inquire_viewport(&xpmin,&xpmax,&ypmin,&ypmax);
  xratio = (xpmax - xpmin)/(xw2 - xw1);
  xpmin1 = xpmin + xratio*(x1 - xw1);
  xpmax1 = xpmin + xratio*(x2 - xw1); 
  yratio = (ypmax - ypmin)/(yw2 - yw1);
  ypmin1 = ypmin + yratio*(y1 - yw1);
  ypmax1 = ypmin + yratio*(y2 - yw1);
  g_set_viewport(xpmin1,xpmax1,ypmin1,ypmax1);
  g_set_window(x1,x2,y1,y2);
}

/******************************************************************************
  Restore the original window before call to set_subwindow.
 *****************************************************************************/
void reset_window(double x1,double x2,double y1,double y2)
{
  g_set_viewport(xpmin,xpmax,ypmin,ypmax);
  g_set_window(x1,x2,y1,y2);
}

void setblackbackground(double x1,double x2,double y1,double y2)
{
  g_new_polygon();/* start a new polygon, to be painted later*/
  g_move_abs(.99*x1+.01*x2,.99*y1+.01*y2);
  g_addline(0.99*x2+.01*x1,.99*y1+.01*y2);
  g_addline(0.99*x2+.01*x1,0.99*y2+.01*y1);
  g_addline(.99*x1+.01*x2,0.99*y2+.01*y1);
  g_set_gray(0.0);
  g_set_linestyle(SOLID);
  g_stroke();
  g_set_gray(0.0);
  g_fill();
}

void colorset(int colour)				    /* tjm 05/07/96 */
{
  switch(colour)
  {
    case 0:
      g_set_color(1.0,0.0,0.0);
      break;
    case 1:
      g_set_color(0.0,1.0,0.0);
      break;
    case 2:
      g_set_color(0.0,0.0,1.0);
      break;
    default:						    /* jwb 06/03/96 */
      g_set_color(1.0,1.0,1.0);				    /* jwb 06/03/96 */
      break;						    /* jwb 06/03/96 */
  }
}
