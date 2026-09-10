
#define max( a , b )  ( (a) > (b) ? (a) : (b) )
#define min( a , b )  ( (a) < (b) ? (a) : (b) )

#define ETA 2.2204460492503131e-16
#define SQRTETA 1.4901161193847656e-08
#define CCMAX  0.3
#define MAXCOR 3
#define MSBP 20
#define MXNCF 10
#define RATIO 5.0

extern double   sm1[13];
/* newly added static variables */

struct lsoda_common_t {
	double **yh, **wm, *ewt, *savf, *acor;
	int     *ipvt;
	void * memory;

	/* static variables for lsoda() */

	double   h, hu, rc, tn;
	double   tsw, pdnorm;

	/* no static variable for prja(), solsy() */
	/* static variables for stoda() */

	double   crate, el[14];
#ifdef CFODE_STATIC
	double (*elco)[14], (*tesco)[4];
#else
	double elco[13][14], tesco[13][4];
#endif
	double hold, rmax;
	double   pdest, pdlast;

	/* static variables for various vectors and the Jacobian. */

	int      ialth, ipup, nslp;
	int      icount, irflag;
	int      imxer;
	int      illin, nhnil, nslast,
					jcur, meth, mused, nq, nst,
					ncf, nfe, nje, nqu, miter;
	/* LINE fix: jstart was a LOCAL of lsoda(), rebuilt as 1 by block d on every
	   continuation call. In the Fortran it lives in a common block: dstoda
	   leaves it at 1 and the driver overrides it with -1 when a method switch
	   has to be completed on the next step, so that -1 SURVIVES the return.
	   Dropping it is invisible under itask=1 over sparse output times, where
	   the switch is normally completed inside the same call, and fatal in a
	   stepping mode (itask=2/5): every step returns, the switch is never
	   completed, and the elco tables stay on the old method. */
	int      jstart;
};
#define _C(x) (ctx->common->x)
