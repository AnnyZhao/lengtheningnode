#define abs(x)  (((x) > 0) ? (x) : (-(x)))
#define min(x,y)  (((x) < (y)) ? (x) : (y))
#define max(x,y) (((x) > (y)) ? (x) : (y))
#define amod(x,y) ((x) - ((int)((x)/(y)))*(y))

/* #   define HUGEI 1.79769313486231700e+308 */	/* max double in IEEE format */
#define HUGEI 2.e+99 /* above value didn't work with 2.0 op sys */
