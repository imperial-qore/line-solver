/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * pfqn_ca.c
 *
 * Code generation for function 'pfqn_ca'
 *
 */

/* Include files */
#include "pfqn_ca.h"
#include "applyScalarFunctionInPlace.h"
#include "eml_int_forloop_overflow_check.h"
#include "log.h"
#include "pfqn_panacea.h"
#include "pfqn_panacea_data.h"
#include "pfqn_panacea_emxutil.h"
#include "pfqn_panacea_types.h"
#include "prod.h"
#include "rt_nonfinite.h"
#include "sum.h"
#include "mwmathutil.h"
#include "omp.h"
#include <emmintrin.h>

/* Variable Definitions */
static emlrtRSInfo
    qc_emlrtRSI =
        {
            62,        /* lineNo */
            "pfqn_ca", /* fcnName */
            "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
            "pfqn_ca.m" /* pathName */
};

static emlrtRSInfo
    rc_emlrtRSI =
        {
            59,        /* lineNo */
            "pfqn_ca", /* fcnName */
            "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
            "pfqn_ca.m" /* pathName */
};

static emlrtRSInfo
    sc_emlrtRSI =
        {
            53,        /* lineNo */
            "pfqn_ca", /* fcnName */
            "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
            "pfqn_ca.m" /* pathName */
};

static emlrtRSInfo
    tc_emlrtRSI =
        {
            47,        /* lineNo */
            "pfqn_ca", /* fcnName */
            "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
            "pfqn_ca.m" /* pathName */
};

static emlrtRSInfo
    uc_emlrtRSI =
        {
            46,        /* lineNo */
            "pfqn_ca", /* fcnName */
            "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
            "pfqn_ca.m" /* pathName */
};

static emlrtRSInfo
    vc_emlrtRSI =
        {
            45,        /* lineNo */
            "pfqn_ca", /* fcnName */
            "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
            "pfqn_ca.m" /* pathName */
};

static emlrtRSInfo
    wc_emlrtRSI =
        {
            43,        /* lineNo */
            "pfqn_ca", /* fcnName */
            "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
            "pfqn_ca.m" /* pathName */
};

static emlrtRSInfo
    xc_emlrtRSI =
        {
            36,        /* lineNo */
            "pfqn_ca", /* fcnName */
            "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
            "pfqn_ca.m" /* pathName */
};

static emlrtRSInfo
    yc_emlrtRSI =
        {
            30,        /* lineNo */
            "pfqn_ca", /* fcnName */
            "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
            "pfqn_ca.m" /* pathName */
};

static emlrtRSInfo
    ad_emlrtRSI =
        {
            25,        /* lineNo */
            "pfqn_ca", /* fcnName */
            "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
            "pfqn_ca.m" /* pathName */
};

static emlrtRSInfo
    ed_emlrtRSI =
        {
            73,        /* lineNo */
            "hashpop", /* fcnName */
            "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
            "pfqn_ca.m" /* pathName */
};

static emlrtRSInfo
    fd_emlrtRSI =
        {
            118,  /* lineNo */
            "Fz", /* fcnName */
            "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
            "pfqn_ca.m" /* pathName */
};

static emlrtRSInfo
    gd_emlrtRSI =
        {
            96,      /* lineNo */
            "pprod", /* fcnName */
            "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
            "pfqn_ca.m" /* pathName */
};

static emlrtRSInfo hd_emlrtRSI = {
    15,    /* lineNo */
    "sum", /* fcnName */
    "/usr/local/MATLAB/R2025a/toolbox/eml/lib/matlab/datafun/sum.m" /* pathName
                                                                     */
};

static emlrtBCInfo
    k_emlrtBCI =
        {
            -1,        /* iFirst */
            -1,        /* iLast */
            73,        /* lineNo */
            29,        /* colNo */
            "N",       /* aName */
            "hashpop", /* fName */
            "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
            "pfqn_ca.m", /* pName */
            0            /* checkKind */
};

static emlrtECInfo
    d_emlrtECI =
        {
            2,       /* nDims */
            96,      /* lineNo */
            8,       /* colNo */
            "pprod", /* fName */
            "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
            "pfqn_ca.m" /* pName */
};

static emlrtBCInfo
    l_emlrtBCI =
        {
            -1,   /* iFirst */
            -1,   /* iLast */
            124,  /* lineNo */
            10,   /* colNo */
            "Z",  /* aName */
            "Fz", /* fName */
            "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
            "pfqn_ca.m", /* pName */
            0            /* checkKind */
};

static emlrtECInfo
    e_emlrtECI =
        {
            2,         /* nDims */
            25,        /* lineNo */
            34,        /* colNo */
            "pfqn_ca", /* fName */
            "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
            "pfqn_ca.m" /* pName */
};

static emlrtDCInfo
    emlrtDCI =
        {
            43,        /* lineNo */
            10,        /* colNo */
            "pfqn_ca", /* fName */
            "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
            "pfqn_ca.m", /* pName */
            1            /* checkKind */
};

static emlrtDCInfo
    b_emlrtDCI =
        {
            43,        /* lineNo */
            14,        /* colNo */
            "pfqn_ca", /* fName */
            "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
            "pfqn_ca.m", /* pName */
            1            /* checkKind */
};

static emlrtRTEInfo
    j_emlrtRTEI =
        {
            48,        /* lineNo */
            11,        /* colNo */
            "pfqn_ca", /* fName */
            "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
            "pfqn_ca.m" /* pName */
};

static emlrtBCInfo
    m_emlrtBCI =
        {
            -1,        /* iFirst */
            -1,        /* iLast */
            61,        /* lineNo */
            6,         /* colNo */
            "G",       /* aName */
            "pfqn_ca", /* fName */
            "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
            "pfqn_ca.m", /* pName */
            0            /* checkKind */
};

static emlrtBCInfo
    n_emlrtBCI =
        {
            -1,        /* iFirst */
            -1,        /* iLast */
            61,        /* lineNo */
            10,        /* colNo */
            "G",       /* aName */
            "pfqn_ca", /* fName */
            "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
            "pfqn_ca.m", /* pName */
            0            /* checkKind */
};

static emlrtBCInfo
    o_emlrtBCI =
        {
            -1,        /* iFirst */
            -1,        /* iLast */
            47,        /* lineNo */
            9,         /* colNo */
            "G",       /* aName */
            "pfqn_ca", /* fName */
            "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
            "pfqn_ca.m", /* pName */
            0            /* checkKind */
};

static emlrtDCInfo
    c_emlrtDCI =
        {
            47,        /* lineNo */
            9,         /* colNo */
            "pfqn_ca", /* fName */
            "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
            "pfqn_ca.m", /* pName */
            1            /* checkKind */
};

static emlrtBCInfo
    p_emlrtBCI =
        {
            -1,        /* iFirst */
            -1,        /* iLast */
            73,        /* lineNo */
            39,        /* colNo */
            "n",       /* aName */
            "hashpop", /* fName */
            "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
            "pfqn_ca.m", /* pName */
            0            /* checkKind */
};

static emlrtBCInfo
    q_emlrtBCI =
        {
            -1,   /* iFirst */
            -1,   /* iLast */
            127,  /* lineNo */
            14,   /* colNo */
            "n",  /* aName */
            "Fz", /* fName */
            "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
            "pfqn_ca.m", /* pName */
            0            /* checkKind */
};

static emlrtBCInfo
    r_emlrtBCI =
        {
            -1,        /* iFirst */
            -1,        /* iLast */
            49,        /* lineNo */
            23,        /* colNo */
            "G",       /* aName */
            "pfqn_ca", /* fName */
            "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
            "pfqn_ca.m", /* pName */
            0            /* checkKind */
};

static emlrtBCInfo
    s_emlrtBCI =
        {
            -1,        /* iFirst */
            -1,        /* iLast */
            49,        /* lineNo */
            27,        /* colNo */
            "G",       /* aName */
            "pfqn_ca", /* fName */
            "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
            "pfqn_ca.m", /* pName */
            0            /* checkKind */
};

static emlrtBCInfo
    t_emlrtBCI =
        {
            -1,        /* iFirst */
            -1,        /* iLast */
            49,        /* lineNo */
            11,        /* colNo */
            "G",       /* aName */
            "pfqn_ca", /* fName */
            "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
            "pfqn_ca.m", /* pName */
            0            /* checkKind */
};

static emlrtBCInfo
    u_emlrtBCI =
        {
            -1,        /* iFirst */
            -1,        /* iLast */
            49,        /* lineNo */
            13,        /* colNo */
            "G",       /* aName */
            "pfqn_ca", /* fName */
            "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
            "pfqn_ca.m", /* pName */
            0            /* checkKind */
};

static emlrtBCInfo
    v_emlrtBCI =
        {
            -1,        /* iFirst */
            -1,        /* iLast */
            51,        /* lineNo */
            18,        /* colNo */
            "n",       /* aName */
            "pfqn_ca", /* fName */
            "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
            "pfqn_ca.m", /* pName */
            0            /* checkKind */
};

static emlrtBCInfo
    w_emlrtBCI =
        {
            -1,        /* iFirst */
            -1,        /* iLast */
            52,        /* lineNo */
            26,        /* colNo */
            "n",       /* aName */
            "pfqn_ca", /* fName */
            "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
            "pfqn_ca.m", /* pName */
            0            /* checkKind */
};

static emlrtBCInfo
    x_emlrtBCI =
        {
            -1,        /* iFirst */
            -1,        /* iLast */
            54,        /* lineNo */
            26,        /* colNo */
            "n",       /* aName */
            "pfqn_ca", /* fName */
            "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
            "pfqn_ca.m", /* pName */
            0            /* checkKind */
};

static emlrtBCInfo
    y_emlrtBCI =
        {
            -1,        /* iFirst */
            -1,        /* iLast */
            55,        /* lineNo */
            31,        /* colNo */
            "G",       /* aName */
            "pfqn_ca", /* fName */
            "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
            "pfqn_ca.m", /* pName */
            0            /* checkKind */
};

static emlrtBCInfo
    ab_emlrtBCI =
        {
            -1,        /* iFirst */
            -1,        /* iLast */
            55,        /* lineNo */
            33,        /* colNo */
            "G",       /* aName */
            "pfqn_ca", /* fName */
            "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
            "pfqn_ca.m", /* pName */
            0            /* checkKind */
};

static emlrtBCInfo
    bb_emlrtBCI =
        {
            -1,        /* iFirst */
            -1,        /* iLast */
            55,        /* lineNo */
            43,        /* colNo */
            "L",       /* aName */
            "pfqn_ca", /* fName */
            "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
            "pfqn_ca.m", /* pName */
            0            /* checkKind */
};

static emlrtBCInfo
    cb_emlrtBCI =
        {
            -1,        /* iFirst */
            -1,        /* iLast */
            55,        /* lineNo */
            47,        /* colNo */
            "L",       /* aName */
            "pfqn_ca", /* fName */
            "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
            "pfqn_ca.m", /* pName */
            0            /* checkKind */
};

static emlrtBCInfo
    db_emlrtBCI =
        {
            -1,        /* iFirst */
            -1,        /* iLast */
            55,        /* lineNo */
            52,        /* colNo */
            "G",       /* aName */
            "pfqn_ca", /* fName */
            "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
            "pfqn_ca.m", /* pName */
            0            /* checkKind */
};

static emlrtBCInfo
    eb_emlrtBCI =
        {
            -1,        /* iFirst */
            -1,        /* iLast */
            55,        /* lineNo */
            54,        /* colNo */
            "G",       /* aName */
            "pfqn_ca", /* fName */
            "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
            "pfqn_ca.m", /* pName */
            0            /* checkKind */
};

static emlrtDCInfo
    d_emlrtDCI =
        {
            55,        /* lineNo */
            54,        /* colNo */
            "pfqn_ca", /* fName */
            "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
            "pfqn_ca.m", /* pName */
            1            /* checkKind */
};

static emlrtBCInfo
    fb_emlrtBCI =
        {
            -1,        /* iFirst */
            -1,        /* iLast */
            55,        /* lineNo */
            19,        /* colNo */
            "G",       /* aName */
            "pfqn_ca", /* fName */
            "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
            "pfqn_ca.m", /* pName */
            0            /* checkKind */
};

static emlrtBCInfo
    gb_emlrtBCI =
        {
            -1,        /* iFirst */
            -1,        /* iLast */
            55,        /* lineNo */
            21,        /* colNo */
            "G",       /* aName */
            "pfqn_ca", /* fName */
            "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
            "pfqn_ca.m", /* pName */
            0            /* checkKind */
};

static emlrtBCInfo
    hb_emlrtBCI =
        {
            -1,      /* iFirst */
            -1,      /* iLast */
            102,     /* lineNo */
            16,      /* colNo */
            "n",     /* aName */
            "pprod", /* fName */
            "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
            "pfqn_ca.m", /* pName */
            0            /* checkKind */
};

static emlrtBCInfo
    ib_emlrtBCI =
        {
            -1,      /* iFirst */
            -1,      /* iLast */
            102,     /* lineNo */
            22,      /* colNo */
            "N",     /* aName */
            "pprod", /* fName */
            "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
            "pfqn_ca.m", /* pName */
            0            /* checkKind */
};

static emlrtBCInfo
    jb_emlrtBCI =
        {
            -1,      /* iFirst */
            -1,      /* iLast */
            103,     /* lineNo */
            7,       /* colNo */
            "n",     /* aName */
            "pprod", /* fName */
            "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
            "pfqn_ca.m", /* pName */
            0            /* checkKind */
};

static emlrtBCInfo
    kb_emlrtBCI =
        {
            -1,      /* iFirst */
            -1,      /* iLast */
            110,     /* lineNo */
            8,       /* colNo */
            "n",     /* aName */
            "pprod", /* fName */
            "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
            "pfqn_ca.m", /* pName */
            0            /* checkKind */
};

static emlrtRTEInfo
    jb_emlrtRTEI =
        {
            22,        /* lineNo */
            7,         /* colNo */
            "pfqn_ca", /* fName */
            "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
            "pfqn_ca.m" /* pName */
};

static emlrtRTEInfo
    kb_emlrtRTEI =
        {
            43,        /* lineNo */
            19,        /* colNo */
            "pfqn_ca", /* fName */
            "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
            "pfqn_ca.m" /* pName */
};

static emlrtRTEInfo
    lb_emlrtRTEI =
        {
            25,        /* lineNo */
            34,        /* colNo */
            "pfqn_ca", /* fName */
            "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
            "pfqn_ca.m" /* pName */
};

static emlrtRTEInfo
    mb_emlrtRTEI =
        {
            43,        /* lineNo */
            1,         /* colNo */
            "pfqn_ca", /* fName */
            "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
            "pfqn_ca.m" /* pName */
};

static emlrtRTEInfo
    nb_emlrtRTEI =
        {
            44,        /* lineNo */
            1,         /* colNo */
            "pfqn_ca", /* fName */
            "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
            "pfqn_ca.m" /* pName */
};

static emlrtRTEInfo
    ob_emlrtRTEI =
        {
            73,        /* lineNo */
            25,        /* colNo */
            "pfqn_ca", /* fName */
            "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
            "pfqn_ca.m" /* pName */
};

static emlrtRTEInfo
    pb_emlrtRTEI =
        {
            96,        /* lineNo */
            8,         /* colNo */
            "pfqn_ca", /* fName */
            "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
            "pfqn_ca.m" /* pName */
};

static emlrtRTEInfo
    qb_emlrtRTEI =
        {
            59,        /* lineNo */
            5,         /* colNo */
            "pfqn_ca", /* fName */
            "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
            "pfqn_ca.m" /* pName */
};

static emlrtRTEInfo rb_emlrtRTEI = {
    8,                                                             /* lineNo */
    14,                                                            /* colNo */
    "factln",                                                      /* fName */
    "/home/gcasale/Dropbox/code/line-dev.git/matlab/util/factln.m" /* pName */
};

/* Function Declarations */
static void eq(const emlrtStack *sp, emxArray_boolean_T *in1,
               const emxArray_real_T *in2, const emxArray_real_T *in3);

/* Function Definitions */
static void eq(const emlrtStack *sp, emxArray_boolean_T *in1,
               const emxArray_real_T *in2, const emxArray_real_T *in3)
{
  jmp_buf *volatile emlrtJBStack;
  const real_T *in2_data;
  const real_T *in3_data;
  int32_T eq_numThreads;
  int32_T i;
  int32_T loop_ub;
  int32_T stride_0_1;
  int32_T stride_1_1;
  boolean_T *in1_data;
  in3_data = in3->data;
  in2_data = in2->data;
  stride_0_1 = in1->size[0] * in1->size[1];
  in1->size[0] = 1;
  emxEnsureCapacity_boolean_T(sp, in1, stride_0_1, &pb_emlrtRTEI);
  if (in3->size[1] == 1) {
    loop_ub = in2->size[1];
  } else {
    loop_ub = in3->size[1];
  }
  stride_0_1 = in1->size[0] * in1->size[1];
  in1->size[1] = loop_ub;
  emxEnsureCapacity_boolean_T(sp, in1, stride_0_1, &pb_emlrtRTEI);
  in1_data = in1->data;
  stride_0_1 = (in2->size[1] != 1);
  stride_1_1 = (in3->size[1] != 1);
  if (loop_ub < 1600) {
    for (i = 0; i < loop_ub; i++) {
      in1_data[i] = (in2_data[i * stride_0_1] == in3_data[i * stride_1_1]);
    }
  } else {
    emlrtEnterParallelRegion((emlrtCTX)sp, omp_in_parallel());
    emlrtPushJmpBuf((emlrtCTX)sp, &emlrtJBStack);
    eq_numThreads = emlrtAllocRegionTLSs(
        sp->tls, omp_in_parallel(), omp_get_max_threads(), omp_get_num_procs());
#pragma omp parallel for num_threads(eq_numThreads)

    for (i = 0; i < loop_ub; i++) {
      in1_data[i] = (in2_data[i * stride_0_1] == in3_data[i * stride_1_1]);
    }
    emlrtPopJmpBuf((emlrtCTX)sp, &emlrtJBStack);
    emlrtExitParallelRegion((emlrtCTX)sp, omp_in_parallel());
  }
}

real_T pfqn_ca(const emlrtStack *sp, const emxArray_real_T *L,
               const emxArray_real_T *N)
{
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack d_st;
  emlrtStack e_st;
  emlrtStack f_st;
  emlrtStack g_st;
  emlrtStack h_st;
  emlrtStack i_st;
  emlrtStack st;
  emxArray_boolean_T *x;
  emxArray_real_T *G;
  emxArray_real_T *b_N;
  emxArray_real_T *b_x;
  emxArray_real_T *n;
  const real_T *L_data;
  const real_T *N_data;
  real_T Gn;
  real_T *G_data;
  real_T *b_N_data;
  real_T *n_data;
  int32_T M;
  int32_T R;
  int32_T b_r;
  int32_T c_r;
  int32_T d_r;
  int32_T k;
  boolean_T *x_data;
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
  /*  % @file pfqn_ca.m */
  /*  % @brief Convolution Algorithm for exact normalizing constant computation.
   */
  /* } */
  /* { */
  /*  % @brief Convolution Algorithm for exact normalizing constant computation.
   */
  /*  % @fn pfqn_ca(L, N, Z) */
  /*  % @param L Service demand matrix. */
  /*  % @param N Population vector. */
  /*  % @param Z Think time vector. */
  /*  % @return Gn Normalizing constant. */
  /*  % @return lGn Logarithm of the normalizing constant. */
  /* } */
  R = L->size[1];
  M = L->size[0];
  emxInit_real_T(sp, &G, 2, &mb_emlrtRTEI);
  emxInit_real_T(sp, &n, 2, &nb_emlrtRTEI);
  emxInit_boolean_T(sp, &x, &pb_emlrtRTEI);
  emxInit_real_T(sp, &b_x, 2, &rb_emlrtRTEI);
  emxInit_real_T(sp, &b_N, 2, &kb_emlrtRTEI);
  if (L->size[0] == 0) {
    int32_T idx;
    int32_T loop_ub;
    int32_T nz;
    idx = b_N->size[0] * b_N->size[1];
    b_N->size[0] = 1;
    b_N->size[1] = L->size[1];
    emxEnsureCapacity_real_T(sp, b_N, idx, &jb_emlrtRTEI);
    b_N_data = b_N->data;
    for (k = 0; k < R; k++) {
      b_N_data[k] = 0.0;
    }
    st.site = &ad_emlrtRSI;
    c_sum(&st, b_N, n);
    st.site = &ad_emlrtRSI;
    b_log(&st, n);
    n_data = n->data;
    loop_ub = N->size[1];
    if ((N->size[1] != n->size[1]) &&
        ((N->size[1] != 1) && (n->size[1] != 1))) {
      emlrtDimSizeImpxCheckR2021b(N->size[1], n->size[1], &e_emlrtECI,
                                  (emlrtConstCTX)sp);
    }
    st.site = &ad_emlrtRSI;
    /*  lf=FACTLN(n) */
    /*  Compure the logarithm of n!        */
    /*  */
    /*  Copyright (c) 2012-2026, Imperial College London */
    /*  All rights reserved.   */
    b_st.site = &fb_emlrtRSI;
    idx = b_x->size[0] * b_x->size[1];
    b_x->size[0] = 1;
    b_x->size[1] = N->size[1];
    emxEnsureCapacity_real_T(&b_st, b_x, idx, &o_emlrtRTEI);
    b_N_data = b_x->data;
    nz = (N->size[1] / 2) << 1;
    idx = nz - 2;
    for (k = 0; k <= idx; k += 2) {
      _mm_storeu_pd(&b_N_data[k],
                    _mm_add_pd(_mm_loadu_pd(&N_data[k]), _mm_set1_pd(1.0)));
    }
    for (k = nz; k < loop_ub; k++) {
      b_N_data[k] = N_data[k] + 1.0;
    }
    c_st.site = &gb_emlrtRSI;
    applyScalarFunctionInPlace(&c_st, b_x);
    if (N->size[1] == n->size[1]) {
      idx = b_N->size[0] * b_N->size[1];
      b_N->size[0] = 1;
      b_N->size[1] = N->size[1];
      emxEnsureCapacity_real_T(sp, b_N, idx, &lb_emlrtRTEI);
      b_N_data = b_N->data;
      idx = nz - 2;
      for (k = 0; k <= idx; k += 2) {
        __m128d r;
        r = _mm_loadu_pd(&n_data[k]);
        _mm_storeu_pd(&b_N_data[k], _mm_mul_pd(_mm_loadu_pd(&N_data[k]), r));
      }
      for (k = nz; k < loop_ub; k++) {
        b_N_data[k] = N_data[k] * n_data[k];
      }
      st.site = &ad_emlrtRSI;
      Gn = -b_sum(&st, b_x) + b_sum(&st, b_N);
    } else {
      st.site = &ad_emlrtRSI;
      Gn = binary_expand_op(&st, ad_emlrtRSI, b_x, N, n);
    }
    Gn = muDoubleScalarExp(Gn);
  } else {
    int32_T idx;
    int32_T nz;
    boolean_T exitg1;
    st.site = &yc_emlrtRSI;
    b_st.site = &lc_emlrtRSI;
    c_st.site = &mc_emlrtRSI;
    d_st.site = &nc_emlrtRSI;
    if (N->size[1] < 1) {
      emlrtErrorWithMessageIdR2018a(
          &d_st, &d_emlrtRTEI, "Coder:toolbox:eml_min_or_max_varDimZero",
          "Coder:toolbox:eml_min_or_max_varDimZero", 0);
    }
    e_st.site = &oc_emlrtRSI;
    f_st.site = &pc_emlrtRSI;
    if (N->size[1] > 2) {
      g_st.site = &cc_emlrtRSI;
      if (!muDoubleScalarIsNaN(N_data[0])) {
        idx = 1;
      } else {
        idx = 0;
        h_st.site = &dc_emlrtRSI;
        if (N->size[1] > 2147483646) {
          i_st.site = &x_emlrtRSI;
          check_forloop_overflow_error(&i_st);
        }
        nz = 2;
        exitg1 = false;
        while ((!exitg1) && (nz <= N->size[1])) {
          if (!muDoubleScalarIsNaN(N_data[nz - 1])) {
            idx = nz;
            exitg1 = true;
          } else {
            nz++;
          }
        }
      }
      if (idx != 0) {
        g_st.site = &bc_emlrtRSI;
        h_st.site = &ec_emlrtRSI;
        if ((idx + 1 <= N->size[1]) && (N->size[1] > 2147483646)) {
          i_st.site = &x_emlrtRSI;
          check_forloop_overflow_error(&i_st);
        }
      }
    }
    st.site = &xc_emlrtRSI;
    if (b_sum(&st, N) == 0.0) {
      Gn = 1.0;
    } else {
      int32_T b_loop_ub;
      if ((real_T)L->size[0] + 1.0 != L->size[0] + 1) {
        emlrtIntegerCheckR2012b((real_T)L->size[0] + 1.0, &emlrtDCI,
                                (emlrtConstCTX)sp);
      }
      idx = b_N->size[0] * b_N->size[1];
      b_N->size[0] = 1;
      b_loop_ub = N->size[1];
      b_N->size[1] = N->size[1];
      emxEnsureCapacity_real_T(sp, b_N, idx, &kb_emlrtRTEI);
      b_N_data = b_N->data;
      idx = (N->size[1] / 2) << 1;
      nz = idx - 2;
      for (k = 0; k <= nz; k += 2) {
        _mm_storeu_pd(&b_N_data[k],
                      _mm_add_pd(_mm_loadu_pd(&N_data[k]), _mm_set1_pd(1.0)));
      }
      for (k = idx; k < b_loop_ub; k++) {
        b_N_data[k] = N_data[k] + 1.0;
      }
      st.site = &wc_emlrtRSI;
      Gn = prod(&st, b_N);
      if (Gn != (int32_T)muDoubleScalarFloor(Gn)) {
        emlrtIntegerCheckR2012b(Gn, &b_emlrtDCI, (emlrtConstCTX)sp);
      }
      idx = G->size[0] * G->size[1];
      G->size[0] = L->size[0] + 1;
      G->size[1] = (int32_T)Gn;
      emxEnsureCapacity_real_T(sp, G, idx, &mb_emlrtRTEI);
      G_data = G->data;
      idx = (L->size[0] + 1) * (int32_T)Gn;
      for (k = 0; k < idx; k++) {
        G_data[k] = 1.0;
      }
      /*  stores G across recursion */
      /*  [N]=PPROD(N,N) */
      /*  sequentially generate all vectors n: 0<=n<=N */
      /*  n=pprod(N) - init */
      /*  n=pprod(n,N) - next state */
      idx = n->size[0] * n->size[1];
      n->size[0] = 1;
      n->size[1] = N->size[1];
      emxEnsureCapacity_real_T(sp, n, idx, &nb_emlrtRTEI);
      n_data = n->data;
      idx = N->size[1];
      for (k = 0; k < idx; k++) {
        n_data[k] = 0.0;
      }
      int32_T exitg11;
      do {
        exitg11 = 0;
        st.site = &vc_emlrtRSI;
        if (b_sum(&st, n) != -1.0) {
          real_T idxn;
          int32_T loop_ub;
          st.site = &uc_emlrtRSI;
          /*  IDX=HASHPOP(N,N,R,PRODS) */
          /*  hash a population vector in n: 0<=n<=N */
          idxn = 1.0;
          for (b_r = 0; b_r < b_loop_ub; b_r++) {
            if (b_r < 1) {
              loop_ub = 0;
            } else {
              if (b_r > b_loop_ub) {
                emlrtDynamicBoundsCheckR2012b(b_r, 1, b_loop_ub, &k_emlrtBCI,
                                              &st);
              }
              loop_ub = b_r;
            }
            idx = b_N->size[0] * b_N->size[1];
            b_N->size[0] = 1;
            b_N->size[1] = loop_ub;
            emxEnsureCapacity_real_T(&st, b_N, idx, &ob_emlrtRTEI);
            b_N_data = b_N->data;
            idx = (loop_ub / 2) << 1;
            nz = idx - 2;
            for (k = 0; k <= nz; k += 2) {
              _mm_storeu_pd(&b_N_data[k], _mm_add_pd(_mm_loadu_pd(&N_data[k]),
                                                     _mm_set1_pd(1.0)));
            }
            for (k = idx; k < loop_ub; k++) {
              b_N_data[k] = N_data[k] + 1.0;
            }
            if (b_r + 1 > n->size[1]) {
              emlrtDynamicBoundsCheckR2012b(b_r + 1, 1, n->size[1], &p_emlrtBCI,
                                            &st);
            }
            b_st.site = &ed_emlrtRSI;
            idxn += prod(&b_st, b_N) * n_data[b_r];
            if (*emlrtBreakCheckR2012bFlagVar != 0) {
              emlrtBreakCheckR2012b(&st);
            }
          }
          st.site = &tc_emlrtRSI;
          /*  F=FZ(Z,N) */
          b_st.site = &fd_emlrtRSI;
          if (b_sum(&b_st, n) == 0.0) {
            if (idxn != (int32_T)muDoubleScalarFloor(idxn)) {
              emlrtIntegerCheckR2012b(idxn, &c_emlrtDCI, &st);
            }
            if (((int32_T)idxn < 1) || ((int32_T)idxn > G->size[1])) {
              emlrtDynamicBoundsCheckR2012b((int32_T)idxn, 1, G->size[1],
                                            &o_emlrtBCI, &st);
            }
            G_data[G->size[0] * ((int32_T)idxn - 1)] = 1.0;
          } else {
            idx = 0;
            int32_T exitg2;
            do {
              exitg2 = 0;
              if (idx <= n->size[1] - 1) {
                if (idx + 1 > R) {
                  emlrtDynamicBoundsCheckR2012b(idx + 1, 1, R, &l_emlrtBCI,
                                                &st);
                }
                if (idx + 1 > n->size[1]) {
                  emlrtDynamicBoundsCheckR2012b(idx + 1, 1, n->size[1],
                                                &q_emlrtBCI, &st);
                }
                if (n_data[idx] > 0.0) {
                  if (idxn != (int32_T)muDoubleScalarFloor(idxn)) {
                    emlrtIntegerCheckR2012b(idxn, &c_emlrtDCI, &st);
                  }
                  if (((int32_T)idxn < 1) || ((int32_T)idxn > G->size[1])) {
                    emlrtDynamicBoundsCheckR2012b((int32_T)idxn, 1, G->size[1],
                                                  &o_emlrtBCI, &st);
                  }
                  G_data[G->size[0] * ((int32_T)idxn - 1)] = 0.0;
                  exitg2 = 1;
                } else {
                  idx++;
                }
              } else {
                if (idxn != (int32_T)muDoubleScalarFloor(idxn)) {
                  emlrtIntegerCheckR2012b(idxn, &c_emlrtDCI, &st);
                }
                if (((int32_T)idxn < 1) || ((int32_T)idxn > G->size[1])) {
                  emlrtDynamicBoundsCheckR2012b((int32_T)idxn, 1, G->size[1],
                                                &o_emlrtBCI, &st);
                }
                G_data[G->size[0] * ((int32_T)idxn - 1)] = 1.0;
                exitg2 = 1;
              }
              if (*emlrtBreakCheckR2012bFlagVar != 0) {
                emlrtBreakCheckR2012b(&st);
              }
            } while (exitg2 == 0);
          }
          emlrtForLoopVectorCheckR2021a(2.0, 1.0, (real_T)M + 1.0,
                                        mxDOUBLE_CLASS, M, &j_emlrtRTEI,
                                        (emlrtConstCTX)sp);
          for (b_r = 0; b_r < M; b_r++) {
            if ((int32_T)((uint32_T)b_r + 1U) > G->size[0]) {
              emlrtDynamicBoundsCheckR2012b((int32_T)((uint32_T)b_r + 1U), 1,
                                            G->size[0], &r_emlrtBCI,
                                            (emlrtConstCTX)sp);
            }
            if (((int32_T)idxn < 1) || ((int32_T)idxn > G->size[1])) {
              emlrtDynamicBoundsCheckR2012b((int32_T)idxn, 1, G->size[1],
                                            &s_emlrtBCI, (emlrtConstCTX)sp);
            }
            if (((int32_T)((uint32_T)b_r + 2U) < 1) ||
                ((int32_T)((uint32_T)b_r + 2U) > G->size[0])) {
              emlrtDynamicBoundsCheckR2012b((int32_T)((uint32_T)b_r + 2U), 1,
                                            G->size[0], &t_emlrtBCI,
                                            (emlrtConstCTX)sp);
            }
            if (((int32_T)idxn < 1) || ((int32_T)idxn > G->size[1])) {
              emlrtDynamicBoundsCheckR2012b((int32_T)idxn, 1, G->size[1],
                                            &u_emlrtBCI, (emlrtConstCTX)sp);
            }
            G_data[(b_r + G->size[0] * ((int32_T)idxn - 1)) + 1] =
                G_data[b_r + G->size[0] * ((int32_T)idxn - 1)];
            /*  norm constant with m-1 queues */
            for (c_r = 0; c_r < R; c_r++) {
              if (c_r + 1 > n->size[1]) {
                emlrtDynamicBoundsCheckR2012b(c_r + 1, 1, n->size[1],
                                              &v_emlrtBCI, (emlrtConstCTX)sp);
              }
              Gn = n_data[c_r];
              if (Gn >= 1.0) {
                if (c_r + 1 > n->size[1]) {
                  emlrtDynamicBoundsCheckR2012b(c_r + 1, 1, n->size[1],
                                                &w_emlrtBCI, (emlrtConstCTX)sp);
                }
                n_data[c_r] = Gn - 1.0;
                st.site = &sc_emlrtRSI;
                /*  IDX=HASHPOP(N,N,R,PRODS) */
                /*  hash a population vector in n: 0<=n<=N */
                Gn = 1.0;
                for (d_r = 0; d_r < b_loop_ub; d_r++) {
                  if (d_r < 1) {
                    loop_ub = 0;
                  } else {
                    if (d_r > b_loop_ub) {
                      emlrtDynamicBoundsCheckR2012b(d_r, 1, b_loop_ub,
                                                    &k_emlrtBCI, &st);
                    }
                    loop_ub = d_r;
                  }
                  idx = b_N->size[0] * b_N->size[1];
                  b_N->size[0] = 1;
                  b_N->size[1] = loop_ub;
                  emxEnsureCapacity_real_T(&st, b_N, idx, &ob_emlrtRTEI);
                  b_N_data = b_N->data;
                  idx = (loop_ub / 2) << 1;
                  nz = idx - 2;
                  for (k = 0; k <= nz; k += 2) {
                    _mm_storeu_pd(
                        &b_N_data[k],
                        _mm_add_pd(_mm_loadu_pd(&N_data[k]), _mm_set1_pd(1.0)));
                  }
                  for (k = idx; k < loop_ub; k++) {
                    b_N_data[k] = N_data[k] + 1.0;
                  }
                  if (d_r + 1 > n->size[1]) {
                    emlrtDynamicBoundsCheckR2012b(d_r + 1, 1, n->size[1],
                                                  &p_emlrtBCI, &st);
                  }
                  b_st.site = &ed_emlrtRSI;
                  Gn += prod(&b_st, b_N) * n_data[d_r];
                  if (*emlrtBreakCheckR2012bFlagVar != 0) {
                    emlrtBreakCheckR2012b(&st);
                  }
                }
                if (c_r + 1 > n->size[1]) {
                  emlrtDynamicBoundsCheckR2012b(c_r + 1, 1, n->size[1],
                                                &x_emlrtBCI, (emlrtConstCTX)sp);
                }
                n_data[c_r]++;
                if (((int32_T)((uint32_T)b_r + 2U) < 1) ||
                    ((int32_T)((uint32_T)b_r + 2U) > G->size[0])) {
                  emlrtDynamicBoundsCheckR2012b((int32_T)((uint32_T)b_r + 2U),
                                                1, G->size[0], &y_emlrtBCI,
                                                (emlrtConstCTX)sp);
                }
                if (((int32_T)idxn < 1) || ((int32_T)idxn > G->size[1])) {
                  emlrtDynamicBoundsCheckR2012b((int32_T)idxn, 1, G->size[1],
                                                &ab_emlrtBCI,
                                                (emlrtConstCTX)sp);
                }
                if ((int32_T)((uint32_T)b_r + 1U) > M) {
                  emlrtDynamicBoundsCheckR2012b((int32_T)((uint32_T)b_r + 1U),
                                                1, M, &bb_emlrtBCI,
                                                (emlrtConstCTX)sp);
                }
                if (c_r + 1 > R) {
                  emlrtDynamicBoundsCheckR2012b(c_r + 1, 1, R, &cb_emlrtBCI,
                                                (emlrtConstCTX)sp);
                }
                if (((int32_T)((uint32_T)b_r + 2U) < 1) ||
                    ((int32_T)((uint32_T)b_r + 2U) > G->size[0])) {
                  emlrtDynamicBoundsCheckR2012b((int32_T)((uint32_T)b_r + 2U),
                                                1, G->size[0], &db_emlrtBCI,
                                                (emlrtConstCTX)sp);
                }
                if (Gn != (int32_T)muDoubleScalarFloor(Gn)) {
                  emlrtIntegerCheckR2012b(Gn, &d_emlrtDCI, (emlrtConstCTX)sp);
                }
                if (((int32_T)Gn < 1) || ((int32_T)Gn > G->size[1])) {
                  emlrtDynamicBoundsCheckR2012b((int32_T)Gn, 1, G->size[1],
                                                &eb_emlrtBCI,
                                                (emlrtConstCTX)sp);
                }
                if (((int32_T)((uint32_T)b_r + 2U) < 1) ||
                    ((int32_T)((uint32_T)b_r + 2U) > G->size[0])) {
                  emlrtDynamicBoundsCheckR2012b((int32_T)((uint32_T)b_r + 2U),
                                                1, G->size[0], &fb_emlrtBCI,
                                                (emlrtConstCTX)sp);
                }
                if (((int32_T)idxn < 1) || ((int32_T)idxn > G->size[1])) {
                  emlrtDynamicBoundsCheckR2012b((int32_T)idxn, 1, G->size[1],
                                                &gb_emlrtBCI,
                                                (emlrtConstCTX)sp);
                }
                G_data[(b_r + G->size[0] * ((int32_T)idxn - 1)) + 1] +=
                    L_data[b_r + L->size[0] * c_r] *
                    G_data[(b_r + G->size[0] * ((int32_T)Gn - 1)) + 1];
              }
              if (*emlrtBreakCheckR2012bFlagVar != 0) {
                emlrtBreakCheckR2012b((emlrtConstCTX)sp);
              }
            }
            if (*emlrtBreakCheckR2012bFlagVar != 0) {
              emlrtBreakCheckR2012b((emlrtConstCTX)sp);
            }
          }
          st.site = &rc_emlrtRSI;
          /*  [N]=PPROD(N,N) */
          /*  sequentially generate all vectors n: 0<=n<=N */
          /*  n=pprod(N) - init */
          /*  n=pprod(n,N) - next state */
          nz = n->size[1];
          if ((n->size[1] != b_loop_ub) &&
              ((n->size[1] != 1) && (b_loop_ub != 1))) {
            emlrtDimSizeImpxCheckR2021b(n->size[1], b_loop_ub, &d_emlrtECI,
                                        &st);
          }
          b_st.site = &gd_emlrtRSI;
          if (n->size[1] == N->size[1]) {
            idx = x->size[0] * x->size[1];
            x->size[0] = 1;
            x->size[1] = n->size[1];
            emxEnsureCapacity_boolean_T(&b_st, x, idx, &pb_emlrtRTEI);
            x_data = x->data;
            for (k = 0; k < nz; k++) {
              x_data[k] = (n_data[k] == N_data[k]);
            }
          } else {
            c_st.site = &gd_emlrtRSI;
            eq(&c_st, x, n, N);
            x_data = x->data;
          }
          c_st.site = &hd_emlrtRSI;
          d_st.site = &o_emlrtRSI;
          idx = x->size[1];
          if (x->size[1] == 0) {
            nz = 0;
          } else {
            e_st.site = &cd_emlrtRSI;
            nz = x_data[0];
            f_st.site = &dd_emlrtRSI;
            if (x->size[1] > 2147483646) {
              g_st.site = &x_emlrtRSI;
              check_forloop_overflow_error(&g_st);
            }
            for (k = 2; k <= idx; k++) {
              nz += x_data[k - 1];
            }
          }
          if (nz == N->size[1]) {
            idx = n->size[0] * n->size[1];
            n->size[0] = 1;
            n->size[1] = 1;
            emxEnsureCapacity_real_T(&st, n, idx, &qb_emlrtRTEI);
            n_data = n->data;
            n_data[0] = -1.0;
          } else {
            idx = N->size[1];
            exitg1 = false;
            while ((!exitg1) && (idx > 0)) {
              if (idx > n->size[1]) {
                emlrtDynamicBoundsCheckR2012b(idx, 1, n->size[1], &hb_emlrtBCI,
                                              &st);
              }
              if (idx > b_loop_ub) {
                emlrtDynamicBoundsCheckR2012b(idx, 1, b_loop_ub, &ib_emlrtBCI,
                                              &st);
              }
              if (n_data[idx - 1] == N_data[idx - 1]) {
                if (idx > n->size[1]) {
                  emlrtDynamicBoundsCheckR2012b(idx, 1, n->size[1],
                                                &jb_emlrtBCI, &st);
                }
                n_data[idx - 1] = 0.0;
                idx--;
              } else {
                exitg1 = true;
              }
              if (*emlrtBreakCheckR2012bFlagVar != 0) {
                emlrtBreakCheckR2012b(&st);
              }
            }
            if (idx != 0) {
              if (idx > n->size[1]) {
                emlrtDynamicBoundsCheckR2012b(idx, 1, n->size[1], &kb_emlrtBCI,
                                              &st);
              }
              n_data[idx - 1]++;
            } else {
              /* n=-1*ones(1,R); */
            }
          }
          if (*emlrtBreakCheckR2012bFlagVar != 0) {
            emlrtBreakCheckR2012b((emlrtConstCTX)sp);
          }
        } else {
          exitg11 = 1;
        }
      } while (exitg11 == 0);
      if (((int32_T)((uint32_T)L->size[0] + 1U) < 1) ||
          ((int32_T)((uint32_T)L->size[0] + 1U) > G->size[0])) {
        emlrtDynamicBoundsCheckR2012b((int32_T)((uint32_T)L->size[0] + 1U), 1,
                                      G->size[0], &m_emlrtBCI,
                                      (emlrtConstCTX)sp);
      }
      if (G->size[1] < 1) {
        emlrtDynamicBoundsCheckR2012b(G->size[1], 1, G->size[1], &n_emlrtBCI,
                                      (emlrtConstCTX)sp);
      }
      Gn = G_data[L->size[0] + G->size[0] * (G->size[1] - 1)];
      st.site = &qc_emlrtRSI;
      if (Gn < 0.0) {
        emlrtErrorWithMessageIdR2018a(
            &st, &emlrtRTEI, "Coder:toolbox:ElFunDomainError",
            "Coder:toolbox:ElFunDomainError", 3, 4, 3, "log");
      }
    }
  }
  emxFree_real_T(sp, &b_N);
  emxFree_real_T(sp, &b_x);
  emxFree_boolean_T(sp, &x);
  emxFree_real_T(sp, &n);
  emxFree_real_T(sp, &G);
  emlrtHeapReferenceStackLeaveFcnR2012b((emlrtConstCTX)sp);
  return Gn;
}

/* End of code generation (pfqn_ca.c) */
