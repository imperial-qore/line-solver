/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * pfqn_panacea.c
 *
 * Code generation for function 'pfqn_panacea'
 *
 */

/* Include files */
#include "pfqn_panacea.h"
#include "applyScalarFunctionInPlace.h"
#include "assertCompatibleDims.h"
#include "div.h"
#include "eml_int_forloop_overflow_check.h"
#include "log.h"
#include "pfqn_ca.h"
#include "pfqn_panacea_data.h"
#include "pfqn_panacea_emxutil.h"
#include "pfqn_panacea_types.h"
#include "repmat.h"
#include "rt_nonfinite.h"
#include "sum.h"
#include "sumMatrixIncludeNaN.h"
#include "blas.h"
#include "mwmathutil.h"
#include "omp.h"
#include <emmintrin.h>
#include <stddef.h>

/* Variable Definitions */
static emlrtRSInfo emlrtRSI =
    {
        27,             /* lineNo */
        "pfqn_panacea", /* fcnName */
        "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
        "pfqn_panacea.m" /* pathName */
};

static emlrtRSInfo b_emlrtRSI =
    {
        28,             /* lineNo */
        "pfqn_panacea", /* fcnName */
        "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
        "pfqn_panacea.m" /* pathName */
};

static emlrtRSInfo c_emlrtRSI =
    {
        32,             /* lineNo */
        "pfqn_panacea", /* fcnName */
        "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
        "pfqn_panacea.m" /* pathName */
};

static emlrtRSInfo d_emlrtRSI =
    {
        33,             /* lineNo */
        "pfqn_panacea", /* fcnName */
        "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
        "pfqn_panacea.m" /* pathName */
};

static emlrtRSInfo e_emlrtRSI =
    {
        36,             /* lineNo */
        "pfqn_panacea", /* fcnName */
        "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
        "pfqn_panacea.m" /* pathName */
};

static emlrtRSInfo f_emlrtRSI =
    {
        37,             /* lineNo */
        "pfqn_panacea", /* fcnName */
        "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
        "pfqn_panacea.m" /* pathName */
};

static emlrtRSInfo g_emlrtRSI =
    {
        38,             /* lineNo */
        "pfqn_panacea", /* fcnName */
        "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
        "pfqn_panacea.m" /* pathName */
};

static emlrtRSInfo h_emlrtRSI =
    {
        49,             /* lineNo */
        "pfqn_panacea", /* fcnName */
        "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
        "pfqn_panacea.m" /* pathName */
};

static emlrtRSInfo i_emlrtRSI =
    {
        55,             /* lineNo */
        "pfqn_panacea", /* fcnName */
        "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
        "pfqn_panacea.m" /* pathName */
};

static emlrtRSInfo j_emlrtRSI =
    {
        57,             /* lineNo */
        "pfqn_panacea", /* fcnName */
        "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
        "pfqn_panacea.m" /* pathName */
};

static emlrtRSInfo k_emlrtRSI =
    {
        61,             /* lineNo */
        "pfqn_panacea", /* fcnName */
        "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
        "pfqn_panacea.m" /* pathName */
};

static emlrtRSInfo l_emlrtRSI =
    {
        87,             /* lineNo */
        "pfqn_panacea", /* fcnName */
        "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
        "pfqn_panacea.m" /* pathName */
};

static emlrtRSInfo m_emlrtRSI =
    {
        90,             /* lineNo */
        "pfqn_panacea", /* fcnName */
        "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
        "pfqn_panacea.m" /* pathName */
};

static emlrtRSInfo db_emlrtRSI =
    {
        18,            /* lineNo */
        "ifWhileCond", /* fcnName */
        "/usr/local/MATLAB/R2025a/toolbox/eml/eml/+coder/+internal/"
        "ifWhileCond.m" /* pathName */
};

static emlrtRSInfo eb_emlrtRSI =
    {
        31,            /* lineNo */
        "checkNoNaNs", /* fcnName */
        "/usr/local/MATLAB/R2025a/toolbox/eml/eml/+coder/+internal/"
        "ifWhileCond.m" /* pathName */
};

static emlrtRSInfo tb_emlrtRSI =
    {
        34,               /* lineNo */
        "rdivide_helper", /* fcnName */
        "/usr/local/MATLAB/R2025a/toolbox/eml/eml/+coder/+internal/"
        "rdivide_helper.m" /* pathName */
};

static emlrtRSInfo ub_emlrtRSI = {
    53,    /* lineNo */
    "div", /* fcnName */
    "/usr/local/MATLAB/R2025a/toolbox/eml/eml/+coder/+internal/div.m" /* pathName
                                                                       */
};

static emlrtRSInfo vb_emlrtRSI = {
    15,    /* lineNo */
    "max", /* fcnName */
    "/usr/local/MATLAB/R2025a/toolbox/eml/lib/matlab/datafun/max.m" /* pathName
                                                                     */
};

static emlrtRSInfo wb_emlrtRSI = {
    73,         /* lineNo */
    "minOrMax", /* fcnName */
    "/usr/local/MATLAB/R2025a/toolbox/eml/eml/+coder/+internal/minOrMax.m" /* pathName
                                                                            */
};

static emlrtRSInfo xb_emlrtRSI = {
    108,       /* lineNo */
    "maximum", /* fcnName */
    "/usr/local/MATLAB/R2025a/toolbox/eml/eml/+coder/+internal/minOrMax.m" /* pathName
                                                                            */
};

static emlrtRSInfo yb_emlrtRSI =
    {
        255,             /* lineNo */
        "unaryMinOrMax", /* fcnName */
        "/usr/local/MATLAB/R2025a/toolbox/eml/eml/+coder/+internal/"
        "unaryMinOrMax.m" /* pathName */
};

static emlrtRSInfo ac_emlrtRSI =
    {
        966,                    /* lineNo */
        "maxRealVectorOmitNaN", /* fcnName */
        "/usr/local/MATLAB/R2025a/toolbox/eml/eml/+coder/+internal/"
        "unaryMinOrMax.m" /* pathName */
};

static emlrtRSInfo fc_emlrtRSI =
    {
        94,                  /* lineNo */
        "eml_mtimes_helper", /* fcnName */
        "/usr/local/MATLAB/R2025a/toolbox/eml/lib/matlab/ops/"
        "eml_mtimes_helper.m" /* pathName */
};

static emlrtRSInfo gc_emlrtRSI =
    {
        69,                  /* lineNo */
        "eml_mtimes_helper", /* fcnName */
        "/usr/local/MATLAB/R2025a/toolbox/eml/lib/matlab/ops/"
        "eml_mtimes_helper.m" /* pathName */
};

static emlrtRSInfo hc_emlrtRSI =
    {
        142,      /* lineNo */
        "mtimes", /* fcnName */
        "/usr/local/MATLAB/R2025a/toolbox/eml/eml/+coder/+internal/+blas/"
        "mtimes.m" /* pathName */
};

static emlrtRSInfo ic_emlrtRSI =
    {
        177,           /* lineNo */
        "mtimes_blas", /* fcnName */
        "/usr/local/MATLAB/R2025a/toolbox/eml/eml/+coder/+internal/+blas/"
        "mtimes.m" /* pathName */
};

static emlrtRSInfo id_emlrtRSI = {
    44,       /* lineNo */
    "mpower", /* fcnName */
    "/usr/local/MATLAB/R2025a/toolbox/eml/lib/matlab/matfun/mpower.m" /* pathName
                                                                       */
};

static emlrtRSInfo jd_emlrtRSI = {
    71,                                                           /* lineNo */
    "power",                                                      /* fcnName */
    "/usr/local/MATLAB/R2025a/toolbox/eml/lib/matlab/ops/power.m" /* pathName */
};

static emlrtRTEInfo b_emlrtRTEI =
    {
        133,                   /* lineNo */
        23,                    /* colNo */
        "dynamic_size_checks", /* fName */
        "/usr/local/MATLAB/R2025a/toolbox/eml/lib/matlab/ops/"
        "eml_mtimes_helper.m" /* pName */
};

static emlrtRTEInfo c_emlrtRTEI =
    {
        138,                   /* lineNo */
        23,                    /* colNo */
        "dynamic_size_checks", /* fName */
        "/usr/local/MATLAB/R2025a/toolbox/eml/lib/matlab/ops/"
        "eml_mtimes_helper.m" /* pName */
};

static emlrtECInfo emlrtECI =
    {
        2,              /* nDims */
        90,             /* lineNo */
        30,             /* colNo */
        "pfqn_panacea", /* fName */
        "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
        "pfqn_panacea.m" /* pName */
};

static emlrtBCInfo emlrtBCI = {
    -1,             /* iFirst */
    -1,             /* iLast */
    57,             /* lineNo */
    24,             /* colNo */
    "beta",         /* aName */
    "pfqn_panacea", /* fName */
    "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
    "pfqn_panacea.m", /* pName */
    0                 /* checkKind */
};

static emlrtECInfo b_emlrtECI =
    {
        2,              /* nDims */
        28,             /* lineNo */
        34,             /* colNo */
        "pfqn_panacea", /* fName */
        "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
        "pfqn_panacea.m" /* pName */
};

static emlrtECInfo c_emlrtECI =
    {
        2,              /* nDims */
        27,             /* lineNo */
        17,             /* colNo */
        "pfqn_panacea", /* fName */
        "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
        "pfqn_panacea.m" /* pName */
};

static emlrtBCInfo b_emlrtBCI = {
    -1,             /* iFirst */
    -1,             /* iLast */
    48,             /* lineNo */
    23,             /* colNo */
    "m",            /* aName */
    "pfqn_panacea", /* fName */
    "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
    "pfqn_panacea.m", /* pName */
    0                 /* checkKind */
};

static emlrtBCInfo c_emlrtBCI = {
    -1,             /* iFirst */
    -1,             /* iLast */
    49,             /* lineNo */
    19,             /* colNo */
    "beta",         /* aName */
    "pfqn_panacea", /* fName */
    "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
    "pfqn_panacea.m", /* pName */
    0                 /* checkKind */
};

static emlrtBCInfo d_emlrtBCI = {
    -1,             /* iFirst */
    -1,             /* iLast */
    54,             /* lineNo */
    23,             /* colNo */
    "m",            /* aName */
    "pfqn_panacea", /* fName */
    "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
    "pfqn_panacea.m", /* pName */
    0                 /* checkKind */
};

static emlrtBCInfo e_emlrtBCI = {
    -1,             /* iFirst */
    -1,             /* iLast */
    55,             /* lineNo */
    24,             /* colNo */
    "beta",         /* aName */
    "pfqn_panacea", /* fName */
    "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
    "pfqn_panacea.m", /* pName */
    0                 /* checkKind */
};

static emlrtBCInfo f_emlrtBCI = {
    -1,             /* iFirst */
    -1,             /* iLast */
    56,             /* lineNo */
    23,             /* colNo */
    "m",            /* aName */
    "pfqn_panacea", /* fName */
    "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
    "pfqn_panacea.m", /* pName */
    0                 /* checkKind */
};

static emlrtBCInfo g_emlrtBCI = {
    -1,             /* iFirst */
    -1,             /* iLast */
    60,             /* lineNo */
    31,             /* colNo */
    "m",            /* aName */
    "pfqn_panacea", /* fName */
    "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
    "pfqn_panacea.m", /* pName */
    0                 /* checkKind */
};

static emlrtBCInfo h_emlrtBCI = {
    -1,             /* iFirst */
    -1,             /* iLast */
    60,             /* lineNo */
    39,             /* colNo */
    "m",            /* aName */
    "pfqn_panacea", /* fName */
    "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
    "pfqn_panacea.m", /* pName */
    0                 /* checkKind */
};

static emlrtBCInfo i_emlrtBCI = {
    -1,             /* iFirst */
    -1,             /* iLast */
    61,             /* lineNo */
    34,             /* colNo */
    "beta",         /* aName */
    "pfqn_panacea", /* fName */
    "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
    "pfqn_panacea.m", /* pName */
    0                 /* checkKind */
};

static emlrtBCInfo j_emlrtBCI = {
    -1,             /* iFirst */
    -1,             /* iLast */
    61,             /* lineNo */
    44,             /* colNo */
    "beta",         /* aName */
    "pfqn_panacea", /* fName */
    "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
    "pfqn_panacea.m", /* pName */
    0                 /* checkKind */
};

static emlrtRTEInfo m_emlrtRTEI =
    {
        25,             /* lineNo */
        5,              /* colNo */
        "pfqn_panacea", /* fName */
        "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
        "pfqn_panacea.m" /* pName */
};

static emlrtRTEInfo n_emlrtRTEI =
    {
        27,             /* lineNo */
        17,             /* colNo */
        "pfqn_panacea", /* fName */
        "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
        "pfqn_panacea.m" /* pName */
};

static emlrtRTEInfo p_emlrtRTEI =
    {
        32,             /* lineNo */
        1,              /* colNo */
        "pfqn_panacea", /* fName */
        "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
        "pfqn_panacea.m" /* pName */
};

static emlrtRTEInfo q_emlrtRTEI =
    {
        34,               /* lineNo */
        1,                /* colNo */
        "rdivide_helper", /* fName */
        "/usr/local/MATLAB/R2025a/toolbox/eml/eml/+coder/+internal/"
        "rdivide_helper.m" /* pName */
};

static emlrtRTEInfo r_emlrtRTEI =
    {
        28,             /* lineNo */
        34,             /* colNo */
        "pfqn_panacea", /* fName */
        "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
        "pfqn_panacea.m" /* pName */
};

static emlrtRTEInfo s_emlrtRTEI =
    {
        34,             /* lineNo */
        1,              /* colNo */
        "pfqn_panacea", /* fName */
        "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
        "pfqn_panacea.m" /* pName */
};

static emlrtRTEInfo t_emlrtRTEI =
    {
        35,             /* lineNo */
        1,              /* colNo */
        "pfqn_panacea", /* fName */
        "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
        "pfqn_panacea.m" /* pName */
};

static emlrtRTEInfo u_emlrtRTEI =
    {
        94,                  /* lineNo */
        5,                   /* colNo */
        "eml_mtimes_helper", /* fName */
        "/usr/local/MATLAB/R2025a/toolbox/eml/lib/matlab/ops/"
        "eml_mtimes_helper.m" /* pName */
};

static emlrtRTEInfo v_emlrtRTEI =
    {
        36,             /* lineNo */
        1,              /* colNo */
        "pfqn_panacea", /* fName */
        "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
        "pfqn_panacea.m" /* pName */
};

static emlrtRTEInfo w_emlrtRTEI =
    {
        218,      /* lineNo */
        20,       /* colNo */
        "mtimes", /* fName */
        "/usr/local/MATLAB/R2025a/toolbox/eml/eml/+coder/+internal/+blas/"
        "mtimes.m" /* pName */
};

static emlrtRTEInfo x_emlrtRTEI =
    {
        37,             /* lineNo */
        30,             /* colNo */
        "pfqn_panacea", /* fName */
        "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
        "pfqn_panacea.m" /* pName */
};

static emlrtRTEInfo y_emlrtRTEI =
    {
        48,             /* lineNo */
        5,              /* colNo */
        "pfqn_panacea", /* fName */
        "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
        "pfqn_panacea.m" /* pName */
};

static emlrtRTEInfo ab_emlrtRTEI =
    {
        54,             /* lineNo */
        5,              /* colNo */
        "pfqn_panacea", /* fName */
        "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
        "pfqn_panacea.m" /* pName */
};

static emlrtRTEInfo bb_emlrtRTEI =
    {
        56,             /* lineNo */
        5,              /* colNo */
        "pfqn_panacea", /* fName */
        "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
        "pfqn_panacea.m" /* pName */
};

static emlrtRTEInfo cb_emlrtRTEI =
    {
        60,             /* lineNo */
        13,             /* colNo */
        "pfqn_panacea", /* fName */
        "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
        "pfqn_panacea.m" /* pName */
};

static emlrtRTEInfo db_emlrtRTEI =
    {
        90,             /* lineNo */
        30,             /* colNo */
        "pfqn_panacea", /* fName */
        "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
        "pfqn_panacea.m" /* pName */
};

static emlrtRTEInfo eb_emlrtRTEI =
    {
        19,             /* lineNo */
        19,             /* colNo */
        "pfqn_panacea", /* fName */
        "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
        "pfqn_panacea.m" /* pName */
};

static emlrtRTEInfo fb_emlrtRTEI =
    {
        33,             /* lineNo */
        10,             /* colNo */
        "pfqn_panacea", /* fName */
        "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
        "pfqn_panacea.m" /* pName */
};

static emlrtRSInfo kd_emlrtRSI = {
    54,    /* lineNo */
    "div", /* fcnName */
    "/usr/local/MATLAB/R2025a/toolbox/eml/eml/+coder/+internal/div.m" /* pathName
                                                                       */
};

/* Function Declarations */
static real_T binary_expand_op_1(const emlrtStack *sp, const emlrtRSInfo in1,
                                 const emxArray_real_T *in2,
                                 const emxArray_real_T *in3,
                                 const emxArray_real_T *in4, real_T in5,
                                 const emxArray_real_T *in6);

static void binary_expand_op_2(const emlrtStack *sp, emxArray_boolean_T *in1,
                               const emxArray_real_T *in2,
                               const emxArray_real_T *in3);

/* Function Definitions */
static real_T binary_expand_op_1(const emlrtStack *sp, const emlrtRSInfo in1,
                                 const emxArray_real_T *in2,
                                 const emxArray_real_T *in3,
                                 const emxArray_real_T *in4, real_T in5,
                                 const emxArray_real_T *in6)
{
  jmp_buf *volatile emlrtJBStack;
  emlrtStack st;
  emxArray_real_T *b_in3;
  const real_T *in3_data;
  const real_T *in4_data;
  real_T out1;
  real_T *b_in3_data;
  int32_T binary_expand_op_1_numThreads;
  int32_T i;
  int32_T loop_ub;
  int32_T stride_0_1;
  int32_T stride_1_1;
  st.prev = sp;
  st.tls = sp->tls;
  in4_data = in4->data;
  in3_data = in3->data;
  emlrtHeapReferenceStackEnterFcnR2012b((emlrtConstCTX)sp);
  emxInit_real_T(sp, &b_in3, 2, &db_emlrtRTEI);
  stride_0_1 = b_in3->size[0] * b_in3->size[1];
  b_in3->size[0] = 1;
  if (in4->size[1] == 1) {
    loop_ub = in3->size[1];
  } else {
    loop_ub = in4->size[1];
  }
  b_in3->size[1] = loop_ub;
  emxEnsureCapacity_real_T(sp, b_in3, stride_0_1, &db_emlrtRTEI);
  b_in3_data = b_in3->data;
  stride_0_1 = (in3->size[1] != 1);
  stride_1_1 = (in4->size[1] != 1);
  if (loop_ub < 1600) {
    for (i = 0; i < loop_ub; i++) {
      b_in3_data[i] = in3_data[i * stride_0_1] * in4_data[i * stride_1_1];
    }
  } else {
    emlrtEnterParallelRegion((emlrtCTX)sp, omp_in_parallel());
    emlrtPushJmpBuf((emlrtCTX)sp, &emlrtJBStack);
    binary_expand_op_1_numThreads = emlrtAllocRegionTLSs(
        sp->tls, omp_in_parallel(), omp_get_max_threads(), omp_get_num_procs());
#pragma omp parallel for num_threads(binary_expand_op_1_numThreads)

    for (i = 0; i < loop_ub; i++) {
      b_in3_data[i] = in3_data[i * stride_0_1] * in4_data[i * stride_1_1];
    }
    emlrtPopJmpBuf((emlrtCTX)sp, &emlrtJBStack);
    emlrtExitParallelRegion((emlrtCTX)sp, omp_in_parallel());
  }
  st.site = (emlrtRSInfo *)&in1;
  out1 = ((-b_sum(&st, in2) + b_sum(&st, b_in3)) + in5) - b_sum(&st, in6);
  emxFree_real_T(sp, &b_in3);
  emlrtHeapReferenceStackLeaveFcnR2012b((emlrtConstCTX)sp);
  return out1;
}

static void binary_expand_op_2(const emlrtStack *sp, emxArray_boolean_T *in1,
                               const emxArray_real_T *in2,
                               const emxArray_real_T *in3)
{
  jmp_buf *volatile emlrtJBStack;
  const real_T *in2_data;
  int32_T binary_expand_op_2_numThreads;
  int32_T i;
  int32_T stride_0_1;
  int32_T unnamed_idx_1;
  boolean_T *in1_data;
  in2_data = in2->data;
  unnamed_idx_1 = in3->size[1];
  stride_0_1 = in1->size[0] * in1->size[1];
  in1->size[0] = 1;
  emxEnsureCapacity_boolean_T(sp, in1, stride_0_1, &n_emlrtRTEI);
  if (unnamed_idx_1 == 1) {
    unnamed_idx_1 = in2->size[1];
  }
  stride_0_1 = in1->size[0] * in1->size[1];
  in1->size[1] = unnamed_idx_1;
  emxEnsureCapacity_boolean_T(sp, in1, stride_0_1, &n_emlrtRTEI);
  in1_data = in1->data;
  stride_0_1 = (in2->size[1] != 1);
  if (unnamed_idx_1 < 1600) {
    for (i = 0; i < unnamed_idx_1; i++) {
      in1_data[i] = (in2_data[i * stride_0_1] == 0.0);
    }
  } else {
    emlrtEnterParallelRegion((emlrtCTX)sp, omp_in_parallel());
    emlrtPushJmpBuf((emlrtCTX)sp, &emlrtJBStack);
    binary_expand_op_2_numThreads = emlrtAllocRegionTLSs(
        sp->tls, omp_in_parallel(), omp_get_max_threads(), omp_get_num_procs());
#pragma omp parallel for num_threads(binary_expand_op_2_numThreads)

    for (i = 0; i < unnamed_idx_1; i++) {
      in1_data[i] = (in2_data[i * stride_0_1] == 0.0);
    }
    emlrtPopJmpBuf((emlrtCTX)sp, &emlrtJBStack);
    emlrtExitParallelRegion((emlrtCTX)sp, omp_in_parallel());
  }
}

real_T binary_expand_op(const emlrtStack *sp, const emlrtRSInfo in1,
                        const emxArray_real_T *in2, const emxArray_real_T *in3,
                        const emxArray_real_T *in4)
{
  jmp_buf *volatile emlrtJBStack;
  emlrtStack st;
  emxArray_real_T *b_in3;
  const real_T *in3_data;
  const real_T *in4_data;
  real_T out1;
  real_T *b_in3_data;
  int32_T binary_expand_op_numThreads;
  int32_T i;
  int32_T loop_ub;
  int32_T stride_0_1;
  int32_T stride_1_1;
  st.prev = sp;
  st.tls = sp->tls;
  in4_data = in4->data;
  in3_data = in3->data;
  emlrtHeapReferenceStackEnterFcnR2012b((emlrtConstCTX)sp);
  emxInit_real_T(sp, &b_in3, 2, &r_emlrtRTEI);
  stride_0_1 = b_in3->size[0] * b_in3->size[1];
  b_in3->size[0] = 1;
  if (in4->size[1] == 1) {
    loop_ub = in3->size[1];
  } else {
    loop_ub = in4->size[1];
  }
  b_in3->size[1] = loop_ub;
  emxEnsureCapacity_real_T(sp, b_in3, stride_0_1, &r_emlrtRTEI);
  b_in3_data = b_in3->data;
  stride_0_1 = (in3->size[1] != 1);
  stride_1_1 = (in4->size[1] != 1);
  if (loop_ub < 1600) {
    for (i = 0; i < loop_ub; i++) {
      b_in3_data[i] = in3_data[i * stride_0_1] * in4_data[i * stride_1_1];
    }
  } else {
    emlrtEnterParallelRegion((emlrtCTX)sp, omp_in_parallel());
    emlrtPushJmpBuf((emlrtCTX)sp, &emlrtJBStack);
    binary_expand_op_numThreads = emlrtAllocRegionTLSs(
        sp->tls, omp_in_parallel(), omp_get_max_threads(), omp_get_num_procs());
#pragma omp parallel for num_threads(binary_expand_op_numThreads)

    for (i = 0; i < loop_ub; i++) {
      b_in3_data[i] = in3_data[i * stride_0_1] * in4_data[i * stride_1_1];
    }
    emlrtPopJmpBuf((emlrtCTX)sp, &emlrtJBStack);
    emlrtExitParallelRegion((emlrtCTX)sp, omp_in_parallel());
  }
  st.site = (emlrtRSInfo *)&in1;
  out1 = -b_sum(&st, in2) + b_sum(&st, b_in3);
  emxFree_real_T(sp, &b_in3);
  emlrtHeapReferenceStackLeaveFcnR2012b((emlrtConstCTX)sp);
  return out1;
}

emlrtCTX emlrtGetRootTLSGlobal(void)
{
  return emlrtRootTLSGlobal;
}

void emlrtLockerFunction(EmlrtLockeeFunction aLockee, emlrtConstCTX aTLS,
                         void *aData)
{
  omp_set_lock(&emlrtLockGlobal);
  emlrtCallLockeeFunction(aLockee, aTLS, aData);
  omp_unset_lock(&emlrtLockGlobal);
}

void pfqn_panacea(const emlrtStack *sp, const emxArray_real_T *L,
                  const emxArray_real_T *N, emxArray_real_T *Z, real_T *Gn,
                  real_T *lGn)
{
  __m128d r2;
  ptrdiff_t k_t;
  ptrdiff_t lda_t;
  ptrdiff_t ldb_t;
  ptrdiff_t ldc_t;
  ptrdiff_t m_t;
  ptrdiff_t n_t;
  jmp_buf *volatile emlrtJBStack;
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack d_st;
  emlrtStack e_st;
  emlrtStack f_st;
  emlrtStack g_st;
  emlrtStack h_st;
  emlrtStack i_st;
  emlrtStack st;
  emxArray_boolean_T *b_r;
  emxArray_real_T *alpha;
  emxArray_real_T *b_N;
  emxArray_real_T *b_gamma;
  emxArray_real_T *beta;
  emxArray_real_T *m;
  emxArray_real_T *r;
  emxArray_real_T *z;
  const real_T *L_data;
  const real_T *N_data;
  real_T A1;
  real_T A2;
  real_T *Z_data;
  real_T *alpha_data;
  real_T *beta_data;
  real_T *gamma_data;
  real_T *z_data;
  int32_T b_j;
  int32_T i;
  int32_T idx;
  int32_T j;
  int32_T k;
  int32_T last;
  int32_T p;
  int32_T pfqn_panacea_numThreads;
  int32_T scalarLB;
  char_T TRANSA1;
  char_T TRANSB1;
  boolean_T guard1;
  boolean_T *r1;
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
  h_st.prev = &g_st;
  h_st.tls = g_st.tls;
  i_st.prev = &h_st;
  i_st.tls = h_st.tls;
  N_data = N->data;
  L_data = L->data;
  emlrtHeapReferenceStackEnterFcnR2012b((emlrtConstCTX)sp);
  /* { */
  /*  % @file pfqn_panacea.m */
  /*  % @brief PANACEA (PAth-based Normal Approximation for Closed networks
   * Estimation Algorithm). */
  /* } */
  /* { */
  /*  % @brief PANACEA (PAth-based Normal Approximation for Closed networks
   * Estimation Algorithm). */
  /*  % @fn pfqn_panacea(L, N, Z) */
  /*  % @param L Service demand matrix. */
  /*  % @param N Population vector. */
  /*  % @param Z Think time vector. */
  /*  % @return Gn Normalizing constant. */
  /*  % @return lGn Logarithm of normalizing constant. */
  /* } */
  /*  [GN,LGN]=PFQN_PANACEA(L,N,Z) */
  /*  K = population vector */
  p = L->size[1];
  if (Z->size[1] == 0) {
    scalarLB = Z->size[0] * Z->size[1];
    Z->size[0] = 1;
    idx = N->size[1];
    Z->size[1] = N->size[1];
    emxEnsureCapacity_real_T(sp, Z, scalarLB, &m_emlrtRTEI);
    Z_data = Z->data;
    scalarLB = (N->size[1] / 2) << 1;
    last = scalarLB - 2;
    for (k = 0; k <= last; k += 2) {
      _mm_storeu_pd(&Z_data[k], _mm_add_pd(_mm_mul_pd(_mm_loadu_pd(&N_data[k]),
                                                      _mm_set1_pd(0.0)),
                                           _mm_set1_pd(1.0E-8)));
    }
    for (k = scalarLB; k < idx; k++) {
      Z_data[k] = N_data[k] * 0.0 + 1.0E-8;
    }
  }
  emxInit_real_T(sp, &r, 2, &p_emlrtRTEI);
  emxInit_real_T(sp, &beta, 2, &s_emlrtRTEI);
  emxInit_real_T(sp, &b_gamma, 2, &t_emlrtRTEI);
  emxInit_real_T(sp, &alpha, 2, &v_emlrtRTEI);
  emxInit_real_T(sp, &m, 2, &y_emlrtRTEI);
  emxInit_boolean_T(sp, &b_r, &eb_emlrtRTEI);
  emxInit_real_T(sp, &z, 1, &fb_emlrtRTEI);
  emxInit_real_T(sp, &b_N, 2, &db_emlrtRTEI);
  guard1 = false;
  if ((L->size[0] == 0) || (L->size[1] == 0)) {
    guard1 = true;
  } else {
    boolean_T exitg1;
    boolean_T y;
    st.site = &emlrtRSI;
    sum(&st, L, beta);
    beta_data = beta->data;
    last = beta->size[1];
    if ((beta->size[1] != L->size[1]) &&
        ((beta->size[1] != 1) && (L->size[1] != 1))) {
      emlrtDimSizeImpxCheckR2021b(beta->size[1], L->size[1], &c_emlrtECI,
                                  (emlrtConstCTX)sp);
    }
    if (beta->size[1] == L->size[1]) {
      scalarLB = b_r->size[0] * b_r->size[1];
      b_r->size[0] = 1;
      b_r->size[1] = beta->size[1];
      emxEnsureCapacity_boolean_T(sp, b_r, scalarLB, &n_emlrtRTEI);
      r1 = b_r->data;
      scalarLB = beta->size[1];
      if (beta->size[1] < 1600) {
        for (i = 0; i < last; i++) {
          r1[i] = (beta_data[i] == 0.0);
        }
      } else {
        emlrtEnterParallelRegion((emlrtCTX)sp, omp_in_parallel());
        emlrtPushJmpBuf((emlrtCTX)sp, &emlrtJBStack);
        pfqn_panacea_numThreads =
            emlrtAllocRegionTLSs(sp->tls, omp_in_parallel(),
                                 omp_get_max_threads(), omp_get_num_procs());
#pragma omp parallel for num_threads(pfqn_panacea_numThreads)

        for (i = 0; i < scalarLB; i++) {
          r1[i] = (beta_data[i] == 0.0);
        }
        emlrtPopJmpBuf((emlrtCTX)sp, &emlrtJBStack);
        emlrtExitParallelRegion((emlrtCTX)sp, omp_in_parallel());
      }
    } else {
      st.site = &emlrtRSI;
      binary_expand_op_2(&st, b_r, beta, L);
      r1 = b_r->data;
    }
    st.site = &emlrtRSI;
    y = (b_r->size[1] != 0);
    if (y) {
      b_st.site = &db_emlrtRSI;
      c_st.site = &eb_emlrtRSI;
      if (b_r->size[1] > 2147483646) {
        d_st.site = &x_emlrtRSI;
        check_forloop_overflow_error(&d_st);
      }
      scalarLB = 0;
      exitg1 = false;
      while ((!exitg1) && (scalarLB <= b_r->size[1] - 1)) {
        if (!r1[scalarLB]) {
          y = false;
          exitg1 = true;
        } else {
          scalarLB++;
        }
      }
    }
    if (y) {
      guard1 = true;
    } else {
      real_T Nt;
      int32_T b_loop_ub;
      int32_T b_scalarLB;
      int32_T c_scalarLB;
      int32_T loop_ub;
      st.site = &c_emlrtRSI;
      b_st.site = &c_emlrtRSI;
      repmat(&b_st, Z, L->size[0], r);
      b_st.site = &tb_emlrtRSI;
      c_st.site = &ub_emlrtRSI;
      assertCompatibleDims(&c_st, L, r);
      if ((L->size[0] == r->size[0]) && (L->size[1] == r->size[1])) {
        idx = L->size[0] * L->size[1];
        scalarLB = r->size[0] * r->size[1];
        r->size[0] = L->size[0];
        r->size[1] = p;
        emxEnsureCapacity_real_T(&b_st, r, scalarLB, &p_emlrtRTEI);
        Z_data = r->data;
        scalarLB = (idx / 2) << 1;
        last = scalarLB - 2;
        for (k = 0; k <= last; k += 2) {
          r2 = _mm_loadu_pd(&Z_data[k]);
          _mm_storeu_pd(&Z_data[k], _mm_div_pd(_mm_loadu_pd(&L_data[k]), r2));
        }
        for (k = scalarLB; k < idx; k++) {
          Z_data[k] = L_data[k] / Z_data[k];
        }
      } else {
        c_st.site = &kd_emlrtRSI;
        b_rdivide(&c_st, r, L);
        Z_data = r->data;
      }
      loop_ub = r->size[0] * r->size[1];
      scalarLB = z->size[0];
      z->size[0] = loop_ub;
      emxEnsureCapacity_real_T(sp, z, scalarLB, &q_emlrtRTEI);
      z_data = z->data;
      b_scalarLB = (loop_ub / 2) << 1;
      scalarLB = b_scalarLB - 2;
      for (k = 0; k <= scalarLB; k += 2) {
        r2 = _mm_loadu_pd(&Z_data[k]);
        _mm_storeu_pd(&z_data[k], _mm_div_pd(_mm_set1_pd(1.0), r2));
      }
      for (k = b_scalarLB; k < loop_ub; k++) {
        z_data[k] = 1.0 / Z_data[k];
      }
      st.site = &d_emlrtRSI;
      b_st.site = &vb_emlrtRSI;
      c_st.site = &wb_emlrtRSI;
      d_st.site = &xb_emlrtRSI;
      if (z->size[0] < 1) {
        emlrtErrorWithMessageIdR2018a(
            &d_st, &d_emlrtRTEI, "Coder:toolbox:eml_min_or_max_varDimZero",
            "Coder:toolbox:eml_min_or_max_varDimZero", 0);
      }
      e_st.site = &yb_emlrtRSI;
      f_st.site = &ac_emlrtRSI;
      last = z->size[0];
      if (z->size[0] <= 2) {
        if (z->size[0] == 1) {
          Nt = z_data[0];
        } else if ((z_data[0] < z_data[1]) ||
                   (muDoubleScalarIsNaN(z_data[0]) &&
                    (!muDoubleScalarIsNaN(z_data[1])))) {
          Nt = z_data[1];
        } else {
          Nt = z_data[0];
        }
      } else {
        g_st.site = &cc_emlrtRSI;
        if (!muDoubleScalarIsNaN(z_data[0])) {
          idx = 1;
        } else {
          idx = 0;
          h_st.site = &dc_emlrtRSI;
          if (z->size[0] > 2147483646) {
            i_st.site = &x_emlrtRSI;
            check_forloop_overflow_error(&i_st);
          }
          scalarLB = 2;
          exitg1 = false;
          while ((!exitg1) && (scalarLB <= last)) {
            if (!muDoubleScalarIsNaN(z_data[scalarLB - 1])) {
              idx = scalarLB;
              exitg1 = true;
            } else {
              scalarLB++;
            }
          }
        }
        if (idx == 0) {
          Nt = z_data[0];
        } else {
          g_st.site = &bc_emlrtRSI;
          Nt = z_data[idx - 1];
          scalarLB = idx + 1;
          h_st.site = &ec_emlrtRSI;
          if ((idx + 1 <= z->size[0]) && (z->size[0] > 2147483646)) {
            i_st.site = &x_emlrtRSI;
            check_forloop_overflow_error(&i_st);
          }
          for (j = scalarLB; j <= last; j++) {
            A1 = z_data[j - 1];
            if (Nt < A1) {
              Nt = A1;
            }
          }
        }
      }
      scalarLB = beta->size[0] * beta->size[1];
      beta->size[0] = 1;
      b_loop_ub = N->size[1];
      beta->size[1] = N->size[1];
      emxEnsureCapacity_real_T(sp, beta, scalarLB, &s_emlrtRTEI);
      beta_data = beta->data;
      c_scalarLB = (N->size[1] / 2) << 1;
      scalarLB = c_scalarLB - 2;
      for (j = 0; j <= scalarLB; j += 2) {
        _mm_storeu_pd(&beta_data[j],
                      _mm_div_pd(_mm_loadu_pd(&N_data[j]), _mm_set1_pd(Nt)));
      }
      for (j = c_scalarLB; j < b_loop_ub; j++) {
        beta_data[j] = N_data[j] / Nt;
      }
      idx = r->size[0];
      scalarLB = b_gamma->size[0] * b_gamma->size[1];
      b_gamma->size[0] = r->size[0];
      b_gamma->size[1] = r->size[1];
      emxEnsureCapacity_real_T(sp, b_gamma, scalarLB, &t_emlrtRTEI);
      gamma_data = b_gamma->data;
      scalarLB = b_scalarLB - 2;
      for (j = 0; j <= scalarLB; j += 2) {
        r2 = _mm_loadu_pd(&Z_data[j]);
        _mm_storeu_pd(&gamma_data[j], _mm_mul_pd(r2, _mm_set1_pd(Nt)));
      }
      for (j = b_scalarLB; j < loop_ub; j++) {
        gamma_data[j] = Z_data[j] * Nt;
      }
      st.site = &e_emlrtRSI;
      b_st.site = &gc_emlrtRSI;
      if (N->size[1] != r->size[1]) {
        if ((N->size[1] == 1) || ((r->size[0] == 1) && (r->size[1] == 1))) {
          emlrtErrorWithMessageIdR2018a(
              &b_st, &b_emlrtRTEI,
              "Coder:toolbox:mtimes_noDynamicScalarExpansion",
              "Coder:toolbox:mtimes_noDynamicScalarExpansion", 0);
        } else {
          emlrtErrorWithMessageIdR2018a(&b_st, &c_emlrtRTEI, "MATLAB:innerdim",
                                        "MATLAB:innerdim", 0);
        }
      }
      b_st.site = &fc_emlrtRSI;
      if ((N->size[1] == 0) || (r->size[0] == 0) || (r->size[1] == 0)) {
        scalarLB = alpha->size[0] * alpha->size[1];
        alpha->size[0] = 1;
        alpha->size[1] = r->size[0];
        emxEnsureCapacity_real_T(&b_st, alpha, scalarLB, &u_emlrtRTEI);
        alpha_data = alpha->data;
        for (k = 0; k < idx; k++) {
          alpha_data[k] = 0.0;
        }
      } else {
        c_st.site = &hc_emlrtRSI;
        d_st.site = &ic_emlrtRSI;
        TRANSB1 = 'T';
        TRANSA1 = 'N';
        A1 = 1.0;
        A2 = 0.0;
        m_t = (ptrdiff_t)1;
        n_t = (ptrdiff_t)r->size[0];
        k_t = (ptrdiff_t)N->size[1];
        lda_t = (ptrdiff_t)1;
        ldb_t = (ptrdiff_t)r->size[0];
        ldc_t = (ptrdiff_t)1;
        scalarLB = alpha->size[0] * alpha->size[1];
        alpha->size[0] = 1;
        alpha->size[1] = r->size[0];
        emxEnsureCapacity_real_T(&d_st, alpha, scalarLB, &w_emlrtRTEI);
        alpha_data = alpha->data;
        dgemm(&TRANSA1, &TRANSB1, &m_t, &n_t, &k_t, &A1, (real_T *)&N_data[0],
              &lda_t, &Z_data[0], &ldb_t, &A2, &alpha_data[0], &ldc_t);
      }
      scalarLB = alpha->size[0] * alpha->size[1];
      alpha->size[0] = 1;
      emxEnsureCapacity_real_T(sp, alpha, scalarLB, &v_emlrtRTEI);
      alpha_data = alpha->data;
      scalarLB = alpha->size[1] - 1;
      last = (alpha->size[1] / 2) << 1;
      idx = last - 2;
      for (j = 0; j <= idx; j += 2) {
        r2 = _mm_loadu_pd(&alpha_data[j]);
        _mm_storeu_pd(&alpha_data[j], _mm_sub_pd(_mm_set1_pd(1.0), r2));
      }
      for (j = last; j <= scalarLB; j++) {
        alpha_data[j] = 1.0 - alpha_data[j];
      }
      st.site = &f_emlrtRSI;
      b_scalarLB = alpha->size[1];
      scalarLB = z->size[0];
      z->size[0] = alpha->size[1];
      emxEnsureCapacity_real_T(&st, z, scalarLB, &x_emlrtRTEI);
      z_data = z->data;
      for (k = 0; k < b_scalarLB; k++) {
        z_data[k] = alpha_data[k];
      }
      b_st.site = &f_emlrtRSI;
      b_repmat(&b_st, z, L->size[1], r);
      Z_data = r->data;
      b_st.site = &tb_emlrtRSI;
      c_st.site = &ub_emlrtRSI;
      assertCompatibleDims(&c_st, b_gamma, r);
      if ((b_gamma->size[0] == r->size[0]) &&
          (b_gamma->size[1] == r->size[1])) {
        scalarLB = b_gamma->size[0] * b_gamma->size[1];
        idx = (scalarLB / 2) << 1;
        last = idx - 2;
        for (j = 0; j <= last; j += 2) {
          __m128d r3;
          r2 = _mm_loadu_pd(&gamma_data[j]);
          r3 = _mm_loadu_pd(&Z_data[j]);
          _mm_storeu_pd(&gamma_data[j], _mm_div_pd(r2, r3));
        }
        for (j = idx; j < scalarLB; j++) {
          gamma_data[j] /= Z_data[j];
        }
      } else {
        c_st.site = &kd_emlrtRSI;
        rdivide(&c_st, b_gamma, r);
      }
      st.site = &g_emlrtRSI;
      b_st.site = &lc_emlrtRSI;
      c_st.site = &mc_emlrtRSI;
      d_st.site = &nc_emlrtRSI;
      if (alpha->size[1] < 1) {
        emlrtErrorWithMessageIdR2018a(
            &d_st, &d_emlrtRTEI, "Coder:toolbox:eml_min_or_max_varDimZero",
            "Coder:toolbox:eml_min_or_max_varDimZero", 0);
      }
      e_st.site = &oc_emlrtRSI;
      f_st.site = &pc_emlrtRSI;
      if (alpha->size[1] <= 2) {
        if (alpha->size[1] == 1) {
          A1 = alpha_data[0];
        } else if ((alpha_data[0] > alpha_data[1]) ||
                   (muDoubleScalarIsNaN(alpha_data[0]) &&
                    (!muDoubleScalarIsNaN(alpha_data[1])))) {
          A1 = alpha_data[1];
        } else {
          A1 = alpha_data[0];
        }
      } else {
        g_st.site = &cc_emlrtRSI;
        if (!muDoubleScalarIsNaN(alpha_data[0])) {
          last = 1;
        } else {
          last = 0;
          h_st.site = &dc_emlrtRSI;
          if (alpha->size[1] > 2147483646) {
            i_st.site = &x_emlrtRSI;
            check_forloop_overflow_error(&i_st);
          }
          scalarLB = 2;
          exitg1 = false;
          while ((!exitg1) && (scalarLB <= b_scalarLB)) {
            if (!muDoubleScalarIsNaN(alpha_data[scalarLB - 1])) {
              last = scalarLB;
              exitg1 = true;
            } else {
              scalarLB++;
            }
          }
        }
        if (last == 0) {
          A1 = alpha_data[0];
        } else {
          g_st.site = &bc_emlrtRSI;
          A1 = alpha_data[last - 1];
          scalarLB = last + 1;
          h_st.site = &ec_emlrtRSI;
          if ((last + 1 <= alpha->size[1]) && (alpha->size[1] > 2147483646)) {
            i_st.site = &x_emlrtRSI;
            check_forloop_overflow_error(&i_st);
          }
          for (k = scalarLB; k <= b_scalarLB; k++) {
            A2 = alpha_data[k - 1];
            if (A1 > A2) {
              A1 = A2;
            }
          }
        }
      }
      if (A1 < 0.0) {
        /*     line_warning(mfilename,'Model is not in normal usage'); */
        *Gn = rtNaN;
        *lGn = rtNaN;
      } else {
        A1 = 0.0;
        for (j = 0; j < p; j++) {
          scalarLB = m->size[0] * m->size[1];
          m->size[0] = 1;
          m->size[1] = p;
          emxEnsureCapacity_real_T(sp, m, scalarLB, &y_emlrtRTEI);
          Z_data = m->data;
          for (k = 0; k < p; k++) {
            Z_data[k] = 0.0;
          }
          if (j + 1 > p) {
            emlrtDynamicBoundsCheckR2012b(j + 1, 1, p, &b_emlrtBCI,
                                          (emlrtConstCTX)sp);
          }
          Z_data[j] = 2.0;
          if (j + 1 > beta->size[1]) {
            emlrtDynamicBoundsCheckR2012b(j + 1, 1, beta->size[1], &c_emlrtBCI,
                                          (emlrtConstCTX)sp);
          }
          st.site = &h_emlrtRSI;
          A1 -= beta_data[j] * pfqn_ca(&st, b_gamma, m);
          if (*emlrtBreakCheckR2012bFlagVar != 0) {
            emlrtBreakCheckR2012b((emlrtConstCTX)sp);
          }
        }
        A2 = 0.0;
        for (b_j = 0; b_j < p; b_j++) {
          real_T A2_tmp;
          scalarLB = m->size[0] * m->size[1];
          m->size[0] = 1;
          m->size[1] = p;
          emxEnsureCapacity_real_T(sp, m, scalarLB, &ab_emlrtRTEI);
          Z_data = m->data;
          for (k = 0; k < p; k++) {
            Z_data[k] = 0.0;
          }
          if (b_j + 1 > p) {
            emlrtDynamicBoundsCheckR2012b(b_j + 1, 1, p, &d_emlrtBCI,
                                          (emlrtConstCTX)sp);
          }
          Z_data[b_j] = 3.0;
          idx = beta->size[1];
          if (b_j + 1 > beta->size[1]) {
            emlrtDynamicBoundsCheckR2012b(b_j + 1, 1, beta->size[1],
                                          &e_emlrtBCI, (emlrtConstCTX)sp);
          }
          A2_tmp = beta_data[b_j];
          st.site = &i_emlrtRSI;
          A2 += 2.0 * A2_tmp * pfqn_ca(&st, b_gamma, m);
          scalarLB = m->size[0] * m->size[1];
          m->size[0] = 1;
          m->size[1] = p;
          emxEnsureCapacity_real_T(sp, m, scalarLB, &bb_emlrtRTEI);
          Z_data = m->data;
          for (k = 0; k < p; k++) {
            Z_data[k] = 0.0;
          }
          if (b_j + 1 > p) {
            emlrtDynamicBoundsCheckR2012b(b_j + 1, 1, p, &f_emlrtBCI,
                                          (emlrtConstCTX)sp);
          }
          Z_data[b_j] = 4.0;
          st.site = &j_emlrtRSI;
          if (b_j + 1 > beta->size[1]) {
            emlrtDynamicBoundsCheckR2012b(b_j + 1, 1, beta->size[1], &emlrtBCI,
                                          &st);
          }
          b_st.site = &id_emlrtRSI;
          c_st.site = &jd_emlrtRSI;
          st.site = &j_emlrtRSI;
          A2 += 3.0 * (beta_data[b_j] * beta_data[b_j]) *
                pfqn_ca(&st, b_gamma, m);
          for (k = 0; k < p; k++) {
            if (k != b_j) {
              scalarLB = m->size[0] * m->size[1];
              m->size[0] = 1;
              m->size[1] = p;
              emxEnsureCapacity_real_T(sp, m, scalarLB, &cb_emlrtRTEI);
              Z_data = m->data;
              for (j = 0; j < p; j++) {
                Z_data[j] = 0.0;
              }
              if (b_j + 1 > p) {
                emlrtDynamicBoundsCheckR2012b(b_j + 1, 1, p, &g_emlrtBCI,
                                              (emlrtConstCTX)sp);
              }
              Z_data[b_j] = 2.0;
              if (k + 1 > m->size[1]) {
                emlrtDynamicBoundsCheckR2012b(k + 1, 1, m->size[1], &h_emlrtBCI,
                                              (emlrtConstCTX)sp);
              }
              Z_data[k] = 2.0;
              if (b_j + 1 > idx) {
                emlrtDynamicBoundsCheckR2012b(b_j + 1, 1, idx, &i_emlrtBCI,
                                              (emlrtConstCTX)sp);
              }
              if (k + 1 > idx) {
                emlrtDynamicBoundsCheckR2012b(k + 1, 1, idx, &j_emlrtBCI,
                                              (emlrtConstCTX)sp);
              }
              st.site = &k_emlrtRSI;
              A2 += 0.5 * A2_tmp * beta_data[k] * pfqn_ca(&st, b_gamma, m);
            }
            if (*emlrtBreakCheckR2012bFlagVar != 0) {
              emlrtBreakCheckR2012b((emlrtConstCTX)sp);
            }
          }
          if (*emlrtBreakCheckR2012bFlagVar != 0) {
            emlrtBreakCheckR2012b((emlrtConstCTX)sp);
          }
        }
        /*  if false */
        /*      A3 = 0; */
        /*      for j=1:p */
        /*          m = zeros(1,p); m(j)=4; */
        /*          A3 = A3 - 6 * beta(j) * pfqn_ca(gammatilde,m); */
        /*          m = zeros(1,p); m(j)=5; */
        /*          A3 = A3 - 20 * beta(j)^2 * pfqn_ca(gammatilde,m); */
        /*          m = zeros(1,p); m(j)=6; */
        /*          A3 = A3 - 15 * beta(j)^3 * pfqn_ca(gammatilde,m); */
        /*          for k=setdiff(1:p,j) */
        /*              m = zeros(1,p); m(j)=4; m(k)=2; */
        /*              A3 = A3 - 2 * beta(j) * beta(k) * pfqn_ca(gammatilde,m);
         */
        /*              m = zeros(1,p); m(j)=2; m(k)=3; */
        /*              A3 = A3 - 3 * beta(j)^2 * beta(k) *
         * pfqn_ca(gammatilde,m); */
        /*              for l=setdiff(1:p,[j,k]) */
        /*                  m = zeros(1,p); m(j)=2; m(k)=2; m(l)=2; */
        /*                  A3 = A3 - (1/6) * beta(j) * beta(k) * beta(l) *
         * pfqn_ca(gammatilde,m); */
        /*              end */
        /*          end */
        /*      end */
        /*  end */
        st.site = &l_emlrtRSI;
        b_st.site = &id_emlrtRSI;
        c_st.site = &jd_emlrtRSI;
        /* , A3/N^3*0]; */
        st.site = &m_emlrtRSI;
        c_sum(&st, Z, beta);
        st.site = &m_emlrtRSI;
        b_log(&st, beta);
        beta_data = beta->data;
        if ((N->size[1] != beta->size[1]) &&
            ((N->size[1] != 1) && (beta->size[1] != 1))) {
          emlrtDimSizeImpxCheckR2021b(N->size[1], beta->size[1], &emlrtECI,
                                      (emlrtConstCTX)sp);
        }
        st.site = &m_emlrtRSI;
        /*  lf=FACTLN(n) */
        /*  Compure the logarithm of n!        */
        /*  */
        /*  Copyright (c) 2012-2026, Imperial College London */
        /*  All rights reserved.   */
        b_st.site = &fb_emlrtRSI;
        scalarLB = m->size[0] * m->size[1];
        m->size[0] = 1;
        m->size[1] = N->size[1];
        emxEnsureCapacity_real_T(&b_st, m, scalarLB, &o_emlrtRTEI);
        Z_data = m->data;
        scalarLB = c_scalarLB - 2;
        for (k = 0; k <= scalarLB; k += 2) {
          _mm_storeu_pd(&Z_data[k],
                        _mm_add_pd(_mm_loadu_pd(&N_data[k]), _mm_set1_pd(1.0)));
        }
        for (k = c_scalarLB; k < b_loop_ub; k++) {
          Z_data[k] = N_data[k] + 1.0;
        }
        real_T dv[3];
        c_st.site = &gb_emlrtRSI;
        applyScalarFunctionInPlace(&c_st, m);
        dv[0] = 1.0;
        dv[1] = A1 / Nt;
        dv[2] = A2 / (Nt * Nt);
        A1 = f_sumColumnB(dv);
        st.site = &m_emlrtRSI;
        if (A1 < 0.0) {
          emlrtErrorWithMessageIdR2018a(
              &st, &emlrtRTEI, "Coder:toolbox:ElFunDomainError",
              "Coder:toolbox:ElFunDomainError", 3, 4, 3, "log");
        }
        A1 = muDoubleScalarLog(A1);
        st.site = &m_emlrtRSI;
        b_log(&st, alpha);
        if (N->size[1] == beta->size[1]) {
          scalarLB = b_N->size[0] * b_N->size[1];
          b_N->size[0] = 1;
          b_N->size[1] = N->size[1];
          emxEnsureCapacity_real_T(sp, b_N, scalarLB, &db_emlrtRTEI);
          Z_data = b_N->data;
          scalarLB = c_scalarLB - 2;
          for (k = 0; k <= scalarLB; k += 2) {
            r2 = _mm_loadu_pd(&beta_data[k]);
            _mm_storeu_pd(&Z_data[k], _mm_mul_pd(_mm_loadu_pd(&N_data[k]), r2));
          }
          for (k = c_scalarLB; k < b_loop_ub; k++) {
            Z_data[k] = N_data[k] * beta_data[k];
          }
          st.site = &m_emlrtRSI;
          *lGn = ((-b_sum(&st, m) + b_sum(&st, b_N)) + A1) - b_sum(&st, alpha);
        } else {
          st.site = &m_emlrtRSI;
          *lGn = binary_expand_op_1(&st, m_emlrtRSI, m, N, beta, A1, alpha);
        }
        *Gn = muDoubleScalarExp(*lGn);
        if (muDoubleScalarIsInf(*lGn) || muDoubleScalarIsNaN(*lGn)) {
          *Gn = rtNaN;
          *lGn = rtNaN;
        }
      }
    }
  }
  if (guard1) {
    st.site = &b_emlrtRSI;
    c_sum(&st, Z, beta);
    st.site = &b_emlrtRSI;
    b_log(&st, beta);
    beta_data = beta->data;
    idx = N->size[1];
    if ((N->size[1] != beta->size[1]) &&
        ((N->size[1] != 1) && (beta->size[1] != 1))) {
      emlrtDimSizeImpxCheckR2021b(N->size[1], beta->size[1], &b_emlrtECI,
                                  (emlrtConstCTX)sp);
    }
    st.site = &b_emlrtRSI;
    /*  lf=FACTLN(n) */
    /*  Compure the logarithm of n!        */
    /*  */
    /*  Copyright (c) 2012-2026, Imperial College London */
    /*  All rights reserved.   */
    b_st.site = &fb_emlrtRSI;
    scalarLB = m->size[0] * m->size[1];
    m->size[0] = 1;
    m->size[1] = N->size[1];
    emxEnsureCapacity_real_T(&b_st, m, scalarLB, &o_emlrtRTEI);
    Z_data = m->data;
    last = (N->size[1] / 2) << 1;
    scalarLB = last - 2;
    for (k = 0; k <= scalarLB; k += 2) {
      _mm_storeu_pd(&Z_data[k],
                    _mm_add_pd(_mm_loadu_pd(&N_data[k]), _mm_set1_pd(1.0)));
    }
    for (k = last; k < idx; k++) {
      Z_data[k] = N_data[k] + 1.0;
    }
    c_st.site = &gb_emlrtRSI;
    applyScalarFunctionInPlace(&c_st, m);
    if (N->size[1] == beta->size[1]) {
      scalarLB = b_N->size[0] * b_N->size[1];
      b_N->size[0] = 1;
      b_N->size[1] = N->size[1];
      emxEnsureCapacity_real_T(sp, b_N, scalarLB, &r_emlrtRTEI);
      Z_data = b_N->data;
      scalarLB = last - 2;
      for (k = 0; k <= scalarLB; k += 2) {
        r2 = _mm_loadu_pd(&beta_data[k]);
        _mm_storeu_pd(&Z_data[k], _mm_mul_pd(_mm_loadu_pd(&N_data[k]), r2));
      }
      for (k = last; k < idx; k++) {
        Z_data[k] = N_data[k] * beta_data[k];
      }
      st.site = &b_emlrtRSI;
      *lGn = -b_sum(&st, m) + b_sum(&st, b_N);
    } else {
      st.site = &b_emlrtRSI;
      *lGn = binary_expand_op(&st, b_emlrtRSI, m, N, beta);
    }
    *Gn = muDoubleScalarExp(*lGn);
  }
  emxFree_real_T(sp, &b_N);
  emxFree_real_T(sp, &z);
  emxFree_boolean_T(sp, &b_r);
  emxFree_real_T(sp, &m);
  emxFree_real_T(sp, &alpha);
  emxFree_real_T(sp, &b_gamma);
  emxFree_real_T(sp, &beta);
  emxFree_real_T(sp, &r);
  emlrtHeapReferenceStackLeaveFcnR2012b((emlrtConstCTX)sp);
}

/* End of code generation (pfqn_panacea.c) */
