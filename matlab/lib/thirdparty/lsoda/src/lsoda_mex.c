/*
 * lsoda_mex.c - MEX gateway for the LSODA ODE solver
 *
 * Solves dy/dt = f(t,y) with automatic stiff/nonstiff switching.
 *
 * Usage from MATLAB:
 *   [T, Y] = lsoda_mex(odefun, tspan, y0, rtol, atol)
 *   [T, Y] = lsoda_mex(odefun, tspan, y0, rtol, atol, mxstep)
 *   [T, Y] = lsoda_mex(odefun, tspan, y0, rtol, atol, mxstep, mxordn, mxords)
 *
 * Inputs:
 *   odefun - function handle: dydt = odefun(t, y)
 *   tspan  - [t0 t1 t2 ... tf] output time points (at least 2 elements)
 *   y0     - initial condition vector (column)
 *   rtol   - relative tolerance (scalar or vector)
 *   atol   - absolute tolerance (scalar or vector)
 *   mxstep - max internal steps per output interval (default 500)
 *   mxordn - max order for nonstiff (Adams) method (1-12, default 0=use max)
 *   mxords - max order for stiff (BDF) method (1-5, default 0=use max)
 *
 * Outputs:
 *   T - column vector of output times
 *   Y - matrix of solution, one row per time point
 *
 * Based on liblsoda by Yu Feng (MIT License), from the original
 * LSODA by Linda R. Petzold and Alan C. Hindmarsh.
 */

#include "mex.h"
#include "lsoda.h"
#include <string.h>

/* User data passed through the LSODA context */
typedef struct {
    mxArray *odefun;     /* MATLAB function handle */
    int neq;             /* number of equations */
    int eval_failed;     /* flag for MATLAB evaluation errors */
} mex_userdata_t;

static void raise_lsoda_failure(struct lsoda_context_t *ctx, double t,
                                double *rtol_arr, double *atol_arr, double *y,
                                double *T_buf, double *Y_buf)
{
    int state = ctx->state;
    char *error_msg = ctx->error;

    ctx->error = NULL;
    lsoda_free(ctx);
    mxFree(rtol_arr);
    mxFree(atol_arr);
    mxFree(y);
    if (T_buf != NULL) {
        mxFree(T_buf);
    }
    if (Y_buf != NULL) {
        mxFree(Y_buf);
    }

    mexErrMsgIdAndTxt("lsoda:failed",
        "LSODA failed at t=%g (state=%d): %s",
        t, state, error_msg ? error_msg : "unknown error");
}

/* ODE right-hand side callback: calls MATLAB function handle */
static int mex_odefun(double t, double *y, double *ydot, void *data)
{
    mex_userdata_t *ud = (mex_userdata_t *)data;
    mxArray *prhs[3], *plhs[1];
    int i;

    if (ud->eval_failed) return -1;

    /* Create MATLAB arrays for t and y */
    prhs[0] = ud->odefun;
    prhs[1] = mxCreateDoubleScalar(t);
    prhs[2] = mxCreateDoubleMatrix(ud->neq, 1, mxREAL);

    /* Copy y into MATLAB array (y is 0-indexed from caller) */
    double *yptr = mxGetPr(prhs[2]);
    for (i = 0; i < ud->neq; i++)
        yptr[i] = y[i];

    /* Call the MATLAB function: dydt = odefun(t, y) */
    if (mexCallMATLAB(1, plhs, 3, prhs, "feval") != 0) {
        ud->eval_failed = 1;
        mxDestroyArray(prhs[1]);
        mxDestroyArray(prhs[2]);
        return -1;
    }

    /* Copy result back */
    double *dydtptr = mxGetPr(plhs[0]);
    for (i = 0; i < ud->neq; i++)
        ydot[i] = dydtptr[i];

    /* Clean up */
    mxDestroyArray(prhs[1]);
    mxDestroyArray(prhs[2]);
    mxDestroyArray(plhs[0]);

    return 0;
}

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    int neq, ntout, i, iout;
    double *tspan, *y0, *rtol_in, *atol_in;
    int mxstep_val, mxordn_val, mxords_val;

    /* Check inputs */
    if (nrhs < 5 || nrhs > 8)
        mexErrMsgIdAndTxt("lsoda:nrhs",
            "Usage: [T,Y] = lsoda_mex(odefun, tspan, y0, rtol, atol [, mxstep [, mxordn, mxords]])");
    if (nlhs > 2)
        mexErrMsgIdAndTxt("lsoda:nlhs", "Too many output arguments.");

    /* Parse odefun */
    if (!mxIsClass(prhs[0], "function_handle"))
        mexErrMsgIdAndTxt("lsoda:odefun", "First argument must be a function handle.");

    /* Parse tspan */
    ntout = (int)mxGetNumberOfElements(prhs[1]);
    if (ntout < 2)
        mexErrMsgIdAndTxt("lsoda:tspan", "tspan must have at least 2 elements.");
    tspan = mxGetPr(prhs[1]);

    /* Parse y0 */
    neq = (int)mxGetNumberOfElements(prhs[2]);
    y0 = mxGetPr(prhs[2]);

    /* Parse rtol */
    int rtol_len = (int)mxGetNumberOfElements(prhs[3]);
    rtol_in = mxGetPr(prhs[3]);

    /* Parse atol */
    int atol_len = (int)mxGetNumberOfElements(prhs[4]);
    atol_in = mxGetPr(prhs[4]);

    /* Parse optional mxstep */
    mxstep_val = 500;
    if (nrhs >= 6) {
        mxstep_val = (int)mxGetScalar(prhs[5]);
        if (mxstep_val <= 0) mxstep_val = 500;
    }

    /* Parse optional mxordn, mxords (0 = use library defaults: 12, 5) */
    mxordn_val = 0;
    mxords_val = 0;
    if (nrhs >= 8) {
        mxordn_val = (int)mxGetScalar(prhs[6]);
        mxords_val = (int)mxGetScalar(prhs[7]);
    }

    /* Build rtol/atol arrays (LSODA expects one per equation) */
    double *rtol_arr = (double *)mxCalloc(neq, sizeof(double));
    double *atol_arr = (double *)mxCalloc(neq, sizeof(double));

    for (i = 0; i < neq; i++) {
        rtol_arr[i] = (rtol_len == 1) ? rtol_in[0] : rtol_in[i];
        atol_arr[i] = (atol_len == 1) ? atol_in[0] : atol_in[i];
    }

    /* Set up LSODA options */
    struct lsoda_opt_t opt;
    memset(&opt, 0, sizeof(opt));
    opt.ixpr = 0;
    opt.rtol = rtol_arr;
    opt.atol = atol_arr;
    opt.itask = 1;
    opt.mxstep = mxstep_val;
    opt.mxordn = mxordn_val;
    opt.mxords = mxords_val;

    /* Set up user data */
    mex_userdata_t ud;
    ud.odefun = (mxArray *)prhs[0];
    ud.neq = neq;
    ud.eval_failed = 0;

    /* Set up LSODA context */
    struct lsoda_context_t ctx;
    memset(&ctx, 0, sizeof(ctx));
    ctx.function = mex_odefun;
    ctx.neq = neq;
    ctx.data = &ud;
    ctx.state = 1;

    if (!lsoda_prepare(&ctx, &opt)) {
        mxFree(rtol_arr);
        mxFree(atol_arr);
        mexErrMsgIdAndTxt("lsoda:prepare", "lsoda_prepare failed: %s",
                          ctx.error ? ctx.error : "unknown error");
    }

    /* Working copy of y (LSODA uses 0-based indexing in this C version) */
    double *y = (double *)mxCalloc(neq, sizeof(double));
    for (i = 0; i < neq; i++)
        y[i] = y0[i];

    double t = tspan[0];

    if (ntout == 2) {
        /*
         * Adaptive stepping mode: when tspan has exactly 2 elements [t0, tf],
         * use itask=2 (one-step mode) to collect all internal adaptive steps.
         * This mimics MATLAB's native ODE solver output behavior.
         */
        opt.itask = 2;
        double tf = tspan[1];

        /* Pre-allocate buffers; grow as needed */
        int buf_cap = 1024;
        int npts = 0;
        double *T_buf = (double *)mxMalloc(buf_cap * sizeof(double));
        double *Y_buf = (double *)mxMalloc(buf_cap * neq * sizeof(double));

        /* Store initial condition */
        T_buf[0] = t;
        for (i = 0; i < neq; i++)
            Y_buf[0 * neq + i] = y[i];
        npts = 1;

        /* Step until we reach or pass tf using itask=2 (one-step mode) */
        double dir = (tf > tspan[0]) ? 1.0 : -1.0;
        while ((tf - t) * dir > 0.0) {
            lsoda(&ctx, y, &t, tf);

            if (ctx.state <= 0) {
                raise_lsoda_failure(&ctx, t, rtol_arr, atol_arr, y, T_buf, Y_buf);
            }

            /* Grow buffers if needed */
            if (npts >= buf_cap) {
                buf_cap *= 2;
                T_buf = (double *)mxRealloc(T_buf, buf_cap * sizeof(double));
                Y_buf = (double *)mxRealloc(Y_buf, buf_cap * neq * sizeof(double));
            }

            /* Only store points that haven't overshot tf */
            if ((tf - t) * dir >= 0.0) {
                T_buf[npts] = t;
                for (i = 0; i < neq; i++)
                    Y_buf[npts * neq + i] = y[i];
                npts++;
            }
        }

        /* Final point: interpolate exactly at tf using itask=1 */
        if (npts == 0 || T_buf[npts-1] != tf) {
            opt.itask = 1;
            lsoda(&ctx, y, &t, tf);
            if (ctx.state > 0) {
                if (npts >= buf_cap) {
                    buf_cap *= 2;
                    T_buf = (double *)mxRealloc(T_buf, buf_cap * sizeof(double));
                    Y_buf = (double *)mxRealloc(Y_buf, buf_cap * neq * sizeof(double));
                }
                T_buf[npts] = tf;
                for (i = 0; i < neq; i++)
                    Y_buf[npts * neq + i] = y[i];
                npts++;
            }
        }

        /* Create output arrays of exact size (column-major for MATLAB) */
        plhs[0] = mxCreateDoubleMatrix(npts, 1, mxREAL);
        plhs[1] = mxCreateDoubleMatrix(npts, neq, mxREAL);
        double *T_out = mxGetPr(plhs[0]);
        double *Y_out = mxGetPr(plhs[1]);

        for (iout = 0; iout < npts; iout++) {
            T_out[iout] = T_buf[iout];
            for (i = 0; i < neq; i++)
                Y_out[iout + i * npts] = Y_buf[iout * neq + i];
        }

        mxFree(T_buf);
        mxFree(Y_buf);

    } else {
        /*
         * Fixed output mode: tspan has >2 elements, integrate to each
         * prescribed output time (itask=1).
         */
        plhs[0] = mxCreateDoubleMatrix(ntout, 1, mxREAL);      /* T */
        plhs[1] = mxCreateDoubleMatrix(ntout, neq, mxREAL);    /* Y */
        double *T_out = mxGetPr(plhs[0]);
        double *Y_out = mxGetPr(plhs[1]);

        /* Store initial condition */
        T_out[0] = t;
        for (i = 0; i < neq; i++)
            Y_out[0 + i * ntout] = y[i]; /* column-major */

        /* Integrate to each output time */
        for (iout = 1; iout < ntout; iout++) {
            double tout = tspan[iout];
            lsoda(&ctx, y, &t, tout);

            if (ctx.state <= 0) {
                raise_lsoda_failure(&ctx, t, rtol_arr, atol_arr, y, NULL, NULL);
            }

            T_out[iout] = t;
            for (i = 0; i < neq; i++)
                Y_out[iout + i * ntout] = y[i];
        }
    }

    /* Clean up */
    /* Prevent lsoda_free from printing "unhandled error" for success */
    if (ctx.error) {
        free(ctx.error);
        ctx.error = NULL;
    }
    lsoda_free(&ctx);
    mxFree(rtol_arr);
    mxFree(atol_arr);
    mxFree(y);
}
