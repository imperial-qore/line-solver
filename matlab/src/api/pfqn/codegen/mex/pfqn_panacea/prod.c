/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * prod.c
 *
 * Code generation for function 'prod'
 *
 */

/* Include files */
#include "prod.h"
#include "eml_int_forloop_overflow_check.h"
#include "pfqn_panacea_data.h"
#include "pfqn_panacea_types.h"
#include "rt_nonfinite.h"

/* Variable Definitions */
static emlrtRSInfo bd_emlrtRSI = {
    11,     /* lineNo */
    "prod", /* fcnName */
    "/usr/local/MATLAB/R2025a/toolbox/eml/lib/matlab/datafun/prod.m" /* pathName
                                                                      */
};

/* Function Definitions */
real_T prod(const emlrtStack *sp, const emxArray_real_T *x)
{
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack d_st;
  emlrtStack e_st;
  emlrtStack st;
  const real_T *x_data;
  real_T y;
  int32_T k;
  int32_T vlen;
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
  x_data = x->data;
  st.site = &bd_emlrtRSI;
  b_st.site = &o_emlrtRSI;
  vlen = x->size[1];
  if (x->size[1] == 0) {
    y = 1.0;
  } else {
    c_st.site = &cd_emlrtRSI;
    y = x_data[0];
    d_st.site = &dd_emlrtRSI;
    if (x->size[1] > 2147483646) {
      e_st.site = &x_emlrtRSI;
      check_forloop_overflow_error(&e_st);
    }
    for (k = 2; k <= vlen; k++) {
      y *= x_data[k - 1];
    }
  }
  return y;
}

/* End of code generation (prod.c) */
