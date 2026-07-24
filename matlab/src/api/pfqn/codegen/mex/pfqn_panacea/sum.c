/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * sum.c
 *
 * Code generation for function 'sum'
 *
 */

/* Include files */
#include "sum.h"
#include "eml_int_forloop_overflow_check.h"
#include "pfqn_panacea_data.h"
#include "pfqn_panacea_emxutil.h"
#include "pfqn_panacea_types.h"
#include "rt_nonfinite.h"
#include "sumMatrixIncludeNaN.h"
#include "omp.h"

/* Variable Definitions */
static emlrtRSInfo n_emlrtRSI = {
    20,    /* lineNo */
    "sum", /* fcnName */
    "/usr/local/MATLAB/R2025a/toolbox/eml/lib/matlab/datafun/sum.m" /* pathName
                                                                     */
};

static emlrtRSInfo p_emlrtRSI = {
    86,                      /* lineNo */
    "combineVectorElements", /* fcnName */
    "/usr/local/MATLAB/R2025a/toolbox/eml/lib/matlab/datafun/private/"
    "combineVectorElements.m" /* pathName */
};

static emlrtRSInfo q_emlrtRSI = {
    107,                /* lineNo */
    "blockedSummation", /* fcnName */
    "/usr/local/MATLAB/R2025a/toolbox/eml/lib/matlab/datafun/private/"
    "blockedSummation.m" /* pathName */
};

static emlrtRSInfo r_emlrtRSI = {
    22,                    /* lineNo */
    "sumMatrixIncludeNaN", /* fcnName */
    "/usr/local/MATLAB/R2025a/toolbox/eml/lib/matlab/datafun/private/"
    "sumMatrixIncludeNaN.m" /* pathName */
};

static emlrtRSInfo s_emlrtRSI = {
    41,                 /* lineNo */
    "sumMatrixColumns", /* fcnName */
    "/usr/local/MATLAB/R2025a/toolbox/eml/lib/matlab/datafun/private/"
    "sumMatrixIncludeNaN.m" /* pathName */
};

static emlrtRSInfo t_emlrtRSI = {
    42,                 /* lineNo */
    "sumMatrixColumns", /* fcnName */
    "/usr/local/MATLAB/R2025a/toolbox/eml/lib/matlab/datafun/private/"
    "sumMatrixIncludeNaN.m" /* pathName */
};

static emlrtRSInfo u_emlrtRSI = {
    50,                 /* lineNo */
    "sumMatrixColumns", /* fcnName */
    "/usr/local/MATLAB/R2025a/toolbox/eml/lib/matlab/datafun/private/"
    "sumMatrixIncludeNaN.m" /* pathName */
};

static emlrtRSInfo v_emlrtRSI = {
    53,                 /* lineNo */
    "sumMatrixColumns", /* fcnName */
    "/usr/local/MATLAB/R2025a/toolbox/eml/lib/matlab/datafun/private/"
    "sumMatrixIncludeNaN.m" /* pathName */
};

static emlrtRSInfo w_emlrtRSI = {
    57,                 /* lineNo */
    "sumMatrixColumns", /* fcnName */
    "/usr/local/MATLAB/R2025a/toolbox/eml/lib/matlab/datafun/private/"
    "sumMatrixIncludeNaN.m" /* pathName */
};

static emlrtRSInfo nb_emlrtRSI = {
    99,                 /* lineNo */
    "blockedSummation", /* fcnName */
    "/usr/local/MATLAB/R2025a/toolbox/eml/lib/matlab/datafun/private/"
    "blockedSummation.m" /* pathName */
};

static emlrtRTEInfo kb_emlrtRTEI = {
    20,                                                             /* lineNo */
    1,                                                              /* colNo */
    "sum",                                                          /* fName */
    "/usr/local/MATLAB/R2025a/toolbox/eml/lib/matlab/datafun/sum.m" /* pName */
};

static emlrtRTEInfo lb_emlrtRTEI = {
    35,                    /* lineNo */
    20,                    /* colNo */
    "sumMatrixIncludeNaN", /* fName */
    "/usr/local/MATLAB/R2025a/toolbox/eml/lib/matlab/datafun/private/"
    "sumMatrixIncludeNaN.m" /* pName */
};

/* Function Definitions */
real_T b_sum(const emlrtStack *sp, const emxArray_real_T *x)
{
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack d_st;
  emlrtStack e_st;
  emlrtStack f_st;
  emlrtStack st;
  emxArray_real_T b_x;
  real_T y;
  int32_T c_x;
  int32_T d_x;
  int32_T e_x;
  int32_T f_x;
  int32_T ib;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  d_st.prev = &c_st;
  d_st.tls = c_st.tls;
  e_st.prev = &d_st;
  e_st.tls = d_st.tls;
  f_st.prev = &e_st;
  f_st.tls = e_st.tls;
  st.site = &n_emlrtRSI;
  b_st.site = &o_emlrtRSI;
  c_st.site = &p_emlrtRSI;
  if (x->size[1] == 0) {
    y = 0.0;
  } else {
    d_st.site = &nb_emlrtRSI;
    e_st.site = &r_emlrtRSI;
    if (x->size[1] < 4096) {
      b_x = *x;
      c_x = x->size[1];
      b_x.size = &c_x;
      b_x.numDimensions = 1;
      f_st.site = &t_emlrtRSI;
      y = c_sumColumnB(&f_st, &b_x, x->size[1]);
    } else {
      int32_T inb;
      int32_T nfb;
      int32_T nleft;
      nfb = (int32_T)((uint32_T)x->size[1] >> 12);
      inb = nfb << 12;
      nleft = x->size[1] - inb;
      b_x = *x;
      d_x = x->size[1];
      b_x.size = &d_x;
      b_x.numDimensions = 1;
      y = b_sumColumnB4(&b_x, 1);
      for (ib = 2; ib <= nfb; ib++) {
        b_x = *x;
        e_x = x->size[1];
        b_x.size = &e_x;
        b_x.numDimensions = 1;
        y += b_sumColumnB4(&b_x, ((ib - 1) << 12) + 1);
      }
      if (nleft > 0) {
        b_x = *x;
        f_x = x->size[1];
        b_x.size = &f_x;
        b_x.numDimensions = 1;
        f_st.site = &w_emlrtRSI;
        y += d_sumColumnB(&f_st, &b_x, nleft, inb + 1);
      }
    }
  }
  return y;
}

void c_sum(const emlrtStack *sp, const emxArray_real_T *x, emxArray_real_T *y)
{
  jmp_buf *volatile emlrtJBStack;
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack d_st;
  emlrtStack e_st;
  emlrtStack f_st;
  emlrtStack g_st;
  emlrtStack st;
  real_T *y_data;
  int32_T c_sum_numThreads;
  int32_T col;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  d_st.prev = &c_st;
  d_st.tls = c_st.tls;
  e_st.prev = &d_st;
  e_st.tls = d_st.tls;
  f_st.prev = &e_st;
  f_st.tls = e_st.tls;
  g_st.prev = &f_st;
  g_st.tls = f_st.tls;
  st.site = &n_emlrtRSI;
  b_st.site = &o_emlrtRSI;
  c_st.site = &p_emlrtRSI;
  if (x->size[1] == 0) {
    y->size[0] = 1;
    y->size[1] = 0;
  } else {
    int32_T i;
    int32_T ncols;
    d_st.site = &q_emlrtRSI;
    e_st.site = &r_emlrtRSI;
    ncols = y->size[0] * y->size[1];
    y->size[0] = 1;
    i = x->size[1];
    y->size[1] = x->size[1];
    emxEnsureCapacity_real_T(&e_st, y, ncols, &lb_emlrtRTEI);
    y_data = y->data;
    ncols = x->size[1];
    f_st.site = &s_emlrtRSI;
    if (x->size[1] > 2147483646) {
      g_st.site = &x_emlrtRSI;
      check_forloop_overflow_error(&g_st);
    }
    if (x->size[1] < 1600) {
      for (col = 0; col < i; col++) {
        y_data[col] = e_sumColumnB(x, col + 1);
      }
    } else {
      emlrtEnterParallelRegion(&e_st, omp_in_parallel());
      emlrtPushJmpBuf(&e_st, &emlrtJBStack);
      c_sum_numThreads =
          emlrtAllocRegionTLSs(e_st.tls, omp_in_parallel(),
                               omp_get_max_threads(), omp_get_num_procs());
#pragma omp parallel for num_threads(c_sum_numThreads)

      for (col = 0; col < ncols; col++) {
        y_data[col] = e_sumColumnB(x, col + 1);
      }
      emlrtPopJmpBuf(&e_st, &emlrtJBStack);
      emlrtExitParallelRegion(&e_st, omp_in_parallel());
    }
  }
}

void sum(const emlrtStack *sp, const emxArray_real_T *x, emxArray_real_T *y)
{
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack d_st;
  emlrtStack e_st;
  emlrtStack f_st;
  emlrtStack g_st;
  emlrtStack st;
  real_T *y_data;
  int32_T col;
  int32_T ib;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  d_st.prev = &c_st;
  d_st.tls = c_st.tls;
  e_st.prev = &d_st;
  e_st.tls = d_st.tls;
  f_st.prev = &e_st;
  f_st.tls = e_st.tls;
  g_st.prev = &f_st;
  g_st.tls = f_st.tls;
  st.site = &n_emlrtRSI;
  b_st.site = &o_emlrtRSI;
  c_st.site = &p_emlrtRSI;
  if ((x->size[0] == 0) || (x->size[1] == 0)) {
    int32_T nfb;
    nfb = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = x->size[1];
    emxEnsureCapacity_real_T(&c_st, y, nfb, &kb_emlrtRTEI);
    y_data = y->data;
    nfb = x->size[1];
    for (col = 0; col < nfb; col++) {
      y_data[col] = 0.0;
    }
  } else {
    int32_T i;
    int32_T nfb;
    d_st.site = &q_emlrtRSI;
    e_st.site = &r_emlrtRSI;
    nfb = y->size[0] * y->size[1];
    y->size[0] = 1;
    i = x->size[1];
    y->size[1] = x->size[1];
    emxEnsureCapacity_real_T(&e_st, y, nfb, &lb_emlrtRTEI);
    y_data = y->data;
    if (x->size[0] < 4096) {
      f_st.site = &s_emlrtRSI;
      if (x->size[1] > 2147483646) {
        g_st.site = &x_emlrtRSI;
        check_forloop_overflow_error(&g_st);
      }
      for (col = 0; col < i; col++) {
        f_st.site = &t_emlrtRSI;
        y_data[col] = sumColumnB(&f_st, x, col + 1, x->size[0]);
      }
    } else {
      int32_T inb;
      int32_T nleft;
      nfb = (int32_T)((uint32_T)x->size[0] >> 12);
      inb = nfb << 12;
      nleft = x->size[0] - inb;
      f_st.site = &u_emlrtRSI;
      if (x->size[1] > 2147483646) {
        g_st.site = &x_emlrtRSI;
        check_forloop_overflow_error(&g_st);
      }
      for (col = 0; col < i; col++) {
        real_T s;
        s = sumColumnB4(x, col + 1, 1);
        f_st.site = &v_emlrtRSI;
        for (ib = 2; ib <= nfb; ib++) {
          s += sumColumnB4(x, col + 1, ((ib - 1) << 12) + 1);
        }
        if (nleft > 0) {
          f_st.site = &w_emlrtRSI;
          s += b_sumColumnB(&f_st, x, col + 1, nleft, inb + 1);
        }
        y_data[col] = s;
      }
    }
  }
}

/* End of code generation (sum.c) */
