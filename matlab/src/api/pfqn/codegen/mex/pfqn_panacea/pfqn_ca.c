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
#include <math.h>

/* Variable Definitions */
static emlrtRSInfo
    sc_emlrtRSI =
        {
            142,       /* lineNo */
            "pfqn_ca", /* fcnName */
            "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
            "pfqn_ca.m" /* pathName */
};

static emlrtRSInfo
    tc_emlrtRSI =
        {
            136,       /* lineNo */
            "pfqn_ca", /* fcnName */
            "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
            "pfqn_ca.m" /* pathName */
};

static emlrtRSInfo
    uc_emlrtRSI =
        {
            131,       /* lineNo */
            "pfqn_ca", /* fcnName */
            "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
            "pfqn_ca.m" /* pathName */
};

static emlrtRSInfo
    vc_emlrtRSI =
        {
            125,       /* lineNo */
            "pfqn_ca", /* fcnName */
            "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
            "pfqn_ca.m" /* pathName */
};

static emlrtRSInfo
    wc_emlrtRSI =
        {
            119,       /* lineNo */
            "pfqn_ca", /* fcnName */
            "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
            "pfqn_ca.m" /* pathName */
};

static emlrtRSInfo
    xc_emlrtRSI =
        {
            118,       /* lineNo */
            "pfqn_ca", /* fcnName */
            "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
            "pfqn_ca.m" /* pathName */
};

static emlrtRSInfo
    yc_emlrtRSI =
        {
            117,       /* lineNo */
            "pfqn_ca", /* fcnName */
            "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
            "pfqn_ca.m" /* pathName */
};

static emlrtRSInfo
    ad_emlrtRSI =
        {
            115,       /* lineNo */
            "pfqn_ca", /* fcnName */
            "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
            "pfqn_ca.m" /* pathName */
};

static emlrtRSInfo
    bd_emlrtRSI =
        {
            111,       /* lineNo */
            "pfqn_ca", /* fcnName */
            "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
            "pfqn_ca.m" /* pathName */
};

static emlrtRSInfo
    cd_emlrtRSI =
        {
            109,       /* lineNo */
            "pfqn_ca", /* fcnName */
            "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
            "pfqn_ca.m" /* pathName */
};

static emlrtRSInfo
    dd_emlrtRSI =
        {
            93,        /* lineNo */
            "pfqn_ca", /* fcnName */
            "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
            "pfqn_ca.m" /* pathName */
};

static emlrtRSInfo
    ed_emlrtRSI =
        {
            85,        /* lineNo */
            "pfqn_ca", /* fcnName */
            "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
            "pfqn_ca.m" /* pathName */
};

static emlrtRSInfo
    fd_emlrtRSI =
        {
            78,        /* lineNo */
            "pfqn_ca", /* fcnName */
            "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
            "pfqn_ca.m" /* pathName */
};

static emlrtRSInfo
    gd_emlrtRSI =
        {
            36,        /* lineNo */
            "pfqn_ca", /* fcnName */
            "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
            "pfqn_ca.m" /* pathName */
};

static emlrtRSInfo
    hd_emlrtRSI =
        {
            30,        /* lineNo */
            "pfqn_ca", /* fcnName */
            "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
            "pfqn_ca.m" /* pathName */
};

static emlrtRSInfo
    id_emlrtRSI =
        {
            25,        /* lineNo */
            "pfqn_ca", /* fcnName */
            "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
            "pfqn_ca.m" /* pathName */
};

static emlrtRSInfo jd_emlrtRSI = {
    13,                                                         /* lineNo */
    "any",                                                      /* fcnName */
    "/usr/local/MATLAB/R2025a/toolbox/eml/lib/matlab/ops/any.m" /* pathName */
};

static emlrtRSInfo kd_emlrtRSI = {
    143,        /* lineNo */
    "allOrAny", /* fcnName */
    "/usr/local/MATLAB/R2025a/toolbox/eml/eml/+coder/+internal/allOrAny.m" /* pathName
                                                                            */
};

static emlrtRSInfo ld_emlrtRSI = {
    12,                                                            /* lineNo */
    "pow2",                                                        /* fcnName */
    "/usr/local/MATLAB/R2025a/toolbox/eml/lib/matlab/elfun/pow2.m" /* pathName
                                                                    */
};

static emlrtRSInfo md_emlrtRSI = {
    48,                    /* lineNo */
    "applyScalarFunction", /* fcnName */
    "/usr/local/MATLAB/R2025a/toolbox/eml/eml/+coder/+internal/"
    "applyScalarFunction.m" /* pathName */
};

static emlrtRSInfo nd_emlrtRSI =
    {
        12,     /* lineNo */
        "pow2", /* fcnName */
        "/usr/local/MATLAB/R2025a/toolbox/eml/eml/+coder/+internal/+scalar/"
        "pow2.m" /* pathName */
};

static emlrtRSInfo
    sd_emlrtRSI =
        {
            153,       /* lineNo */
            "hashpop", /* fcnName */
            "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
            "pfqn_ca.m" /* pathName */
};

static emlrtRSInfo
    td_emlrtRSI =
        {
            198,  /* lineNo */
            "Fz", /* fcnName */
            "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
            "pfqn_ca.m" /* pathName */
};

static emlrtRSInfo
    ud_emlrtRSI =
        {
            176,     /* lineNo */
            "pprod", /* fcnName */
            "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
            "pfqn_ca.m" /* pathName */
};

static emlrtRSInfo vd_emlrtRSI = {
    15,    /* lineNo */
    "sum", /* fcnName */
    "/usr/local/MATLAB/R2025a/toolbox/eml/lib/matlab/datafun/sum.m" /* pathName
                                                                     */
};

static emlrtBCInfo
    l_emlrtBCI =
        {
            -1,        /* iFirst */
            -1,        /* iLast */
            142,       /* lineNo */
            17,        /* colNo */
            "G",       /* aName */
            "pfqn_ca", /* fName */
            "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
            "pfqn_ca.m", /* pName */
            0            /* checkKind */
};

static emlrtBCInfo
    m_emlrtBCI =
        {
            -1,        /* iFirst */
            -1,        /* iLast */
            136,       /* lineNo */
            17,        /* colNo */
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
            136,       /* lineNo */
            13,        /* colNo */
            "G",       /* aName */
            "pfqn_ca", /* fName */
            "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
            "pfqn_ca.m", /* pName */
            0            /* checkKind */
};

static emlrtRTEInfo
    k_emlrtRTEI =
        {
            120,       /* lineNo */
            11,        /* colNo */
            "pfqn_ca", /* fName */
            "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
            "pfqn_ca.m" /* pName */
};

static emlrtDCInfo
    emlrtDCI =
        {
            115,       /* lineNo */
            14,        /* colNo */
            "pfqn_ca", /* fName */
            "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
            "pfqn_ca.m", /* pName */
            1            /* checkKind */
};

static emlrtDCInfo
    b_emlrtDCI =
        {
            115,       /* lineNo */
            10,        /* colNo */
            "pfqn_ca", /* fName */
            "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
            "pfqn_ca.m", /* pName */
            1            /* checkKind */
};

static emlrtECInfo
    d_emlrtECI =
        {
            2,         /* nDims */
            25,        /* lineNo */
            34,        /* colNo */
            "pfqn_ca", /* fName */
            "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
            "pfqn_ca.m" /* pName */
};

static emlrtBCInfo
    o_emlrtBCI =
        {
            -1,   /* iFirst */
            -1,   /* iLast */
            204,  /* lineNo */
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
            2,       /* nDims */
            176,     /* lineNo */
            8,       /* colNo */
            "pprod", /* fName */
            "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
            "pfqn_ca.m" /* pName */
};

static emlrtBCInfo
    p_emlrtBCI =
        {
            -1,        /* iFirst */
            -1,        /* iLast */
            153,       /* lineNo */
            29,        /* colNo */
            "N",       /* aName */
            "hashpop", /* fName */
            "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
            "pfqn_ca.m", /* pName */
            0            /* checkKind */
};

static emlrtBCInfo
    q_emlrtBCI =
        {
            -1,        /* iFirst */
            -1,        /* iLast */
            83,        /* lineNo */
            14,        /* colNo */
            "N",       /* aName */
            "pfqn_ca", /* fName */
            "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
            "pfqn_ca.m", /* pName */
            0            /* checkKind */
};

static emlrtBCInfo
    r_emlrtBCI =
        {
            -1,        /* iFirst */
            -1,        /* iLast */
            84,        /* lineNo */
            18,        /* colNo */
            "L",       /* aName */
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
            84,        /* lineNo */
            20,        /* colNo */
            "L",       /* aName */
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
            85,        /* lineNo */
            36,        /* colNo */
            "L",       /* aName */
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
            85,        /* lineNo */
            38,        /* colNo */
            "L",       /* aName */
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
            85,        /* lineNo */
            27,        /* colNo */
            "N",       /* aName */
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
            119,       /* lineNo */
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
            119,       /* lineNo */
            9,         /* colNo */
            "pfqn_ca", /* fName */
            "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
            "pfqn_ca.m", /* pName */
            1            /* checkKind */
};

static emlrtBCInfo
    x_emlrtBCI =
        {
            -1,        /* iFirst */
            -1,        /* iLast */
            153,       /* lineNo */
            39,        /* colNo */
            "n",       /* aName */
            "hashpop", /* fName */
            "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
            "pfqn_ca.m", /* pName */
            0            /* checkKind */
};

static emlrtBCInfo
    y_emlrtBCI =
        {
            -1,   /* iFirst */
            -1,   /* iLast */
            207,  /* lineNo */
            14,   /* colNo */
            "n",  /* aName */
            "Fz", /* fName */
            "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
            "pfqn_ca.m", /* pName */
            0            /* checkKind */
};

static emlrtBCInfo
    ab_emlrtBCI =
        {
            -1,        /* iFirst */
            -1,        /* iLast */
            121,       /* lineNo */
            23,        /* colNo */
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
            121,       /* lineNo */
            27,        /* colNo */
            "G",       /* aName */
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
            121,       /* lineNo */
            11,        /* colNo */
            "G",       /* aName */
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
            121,       /* lineNo */
            13,        /* colNo */
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
            123,       /* lineNo */
            18,        /* colNo */
            "n",       /* aName */
            "pfqn_ca", /* fName */
            "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
            "pfqn_ca.m", /* pName */
            0            /* checkKind */
};

static emlrtBCInfo
    fb_emlrtBCI =
        {
            -1,        /* iFirst */
            -1,        /* iLast */
            124,       /* lineNo */
            26,        /* colNo */
            "n",       /* aName */
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
            126,       /* lineNo */
            26,        /* colNo */
            "n",       /* aName */
            "pfqn_ca", /* fName */
            "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
            "pfqn_ca.m", /* pName */
            0            /* checkKind */
};

static emlrtBCInfo
    hb_emlrtBCI =
        {
            -1,        /* iFirst */
            -1,        /* iLast */
            127,       /* lineNo */
            31,        /* colNo */
            "G",       /* aName */
            "pfqn_ca", /* fName */
            "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
            "pfqn_ca.m", /* pName */
            0            /* checkKind */
};

static emlrtBCInfo
    ib_emlrtBCI =
        {
            -1,        /* iFirst */
            -1,        /* iLast */
            127,       /* lineNo */
            33,        /* colNo */
            "G",       /* aName */
            "pfqn_ca", /* fName */
            "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
            "pfqn_ca.m", /* pName */
            0            /* checkKind */
};

static emlrtBCInfo
    jb_emlrtBCI =
        {
            -1,        /* iFirst */
            -1,        /* iLast */
            127,       /* lineNo */
            43,        /* colNo */
            "L",       /* aName */
            "pfqn_ca", /* fName */
            "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
            "pfqn_ca.m", /* pName */
            0            /* checkKind */
};

static emlrtBCInfo
    kb_emlrtBCI =
        {
            -1,        /* iFirst */
            -1,        /* iLast */
            127,       /* lineNo */
            47,        /* colNo */
            "L",       /* aName */
            "pfqn_ca", /* fName */
            "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
            "pfqn_ca.m", /* pName */
            0            /* checkKind */
};

static emlrtBCInfo
    lb_emlrtBCI =
        {
            -1,        /* iFirst */
            -1,        /* iLast */
            127,       /* lineNo */
            52,        /* colNo */
            "G",       /* aName */
            "pfqn_ca", /* fName */
            "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
            "pfqn_ca.m", /* pName */
            0            /* checkKind */
};

static emlrtBCInfo
    mb_emlrtBCI =
        {
            -1,        /* iFirst */
            -1,        /* iLast */
            127,       /* lineNo */
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
            127,       /* lineNo */
            54,        /* colNo */
            "pfqn_ca", /* fName */
            "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
            "pfqn_ca.m", /* pName */
            1            /* checkKind */
};

static emlrtBCInfo
    nb_emlrtBCI =
        {
            -1,        /* iFirst */
            -1,        /* iLast */
            127,       /* lineNo */
            19,        /* colNo */
            "G",       /* aName */
            "pfqn_ca", /* fName */
            "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
            "pfqn_ca.m", /* pName */
            0            /* checkKind */
};

static emlrtBCInfo
    ob_emlrtBCI =
        {
            -1,        /* iFirst */
            -1,        /* iLast */
            127,       /* lineNo */
            21,        /* colNo */
            "G",       /* aName */
            "pfqn_ca", /* fName */
            "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
            "pfqn_ca.m", /* pName */
            0            /* checkKind */
};

static emlrtBCInfo
    pb_emlrtBCI =
        {
            -1,      /* iFirst */
            -1,      /* iLast */
            182,     /* lineNo */
            16,      /* colNo */
            "n",     /* aName */
            "pprod", /* fName */
            "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
            "pfqn_ca.m", /* pName */
            0            /* checkKind */
};

static emlrtBCInfo
    qb_emlrtBCI =
        {
            -1,      /* iFirst */
            -1,      /* iLast */
            182,     /* lineNo */
            22,      /* colNo */
            "N",     /* aName */
            "pprod", /* fName */
            "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
            "pfqn_ca.m", /* pName */
            0            /* checkKind */
};

static emlrtBCInfo
    rb_emlrtBCI =
        {
            -1,      /* iFirst */
            -1,      /* iLast */
            183,     /* lineNo */
            7,       /* colNo */
            "n",     /* aName */
            "pprod", /* fName */
            "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
            "pfqn_ca.m", /* pName */
            0            /* checkKind */
};

static emlrtBCInfo
    sb_emlrtBCI =
        {
            -1,      /* iFirst */
            -1,      /* iLast */
            190,     /* lineNo */
            8,       /* colNo */
            "n",     /* aName */
            "pprod", /* fName */
            "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
            "pfqn_ca.m", /* pName */
            0            /* checkKind */
};

static emlrtRTEInfo
    pb_emlrtRTEI =
        {
            22,        /* lineNo */
            7,         /* colNo */
            "pfqn_ca", /* fName */
            "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
            "pfqn_ca.m" /* pName */
};

static emlrtRTEInfo
    qb_emlrtRTEI =
        {
            25,        /* lineNo */
            34,        /* colNo */
            "pfqn_ca", /* fName */
            "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
            "pfqn_ca.m" /* pName */
};

static emlrtRTEInfo
    rb_emlrtRTEI =
        {
            115,       /* lineNo */
            19,        /* colNo */
            "pfqn_ca", /* fName */
            "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
            "pfqn_ca.m" /* pName */
};

static emlrtRTEInfo
    sb_emlrtRTEI =
        {
            115,       /* lineNo */
            1,         /* colNo */
            "pfqn_ca", /* fName */
            "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
            "pfqn_ca.m" /* pName */
};

static emlrtRTEInfo
    tb_emlrtRTEI =
        {
            116,       /* lineNo */
            1,         /* colNo */
            "pfqn_ca", /* fName */
            "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
            "pfqn_ca.m" /* pName */
};

static emlrtRTEInfo
    ub_emlrtRTEI =
        {
            153,       /* lineNo */
            25,        /* colNo */
            "pfqn_ca", /* fName */
            "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
            "pfqn_ca.m" /* pName */
};

static emlrtRTEInfo
    vb_emlrtRTEI =
        {
            176,       /* lineNo */
            8,         /* colNo */
            "pfqn_ca", /* fName */
            "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
            "pfqn_ca.m" /* pName */
};

static emlrtRTEInfo
    wb_emlrtRTEI =
        {
            131,       /* lineNo */
            5,         /* colNo */
            "pfqn_ca", /* fName */
            "/home/gcasale/Dropbox/code/line-dev.git/matlab/src/api/pfqn/"
            "pfqn_ca.m" /* pName */
};

static emlrtRTEInfo xb_emlrtRTEI = {
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
  emxEnsureCapacity_boolean_T(sp, in1, stride_0_1, &vb_emlrtRTEI);
  if (in3->size[1] == 1) {
    loop_ub = in2->size[1];
  } else {
    loop_ub = in3->size[1];
  }
  stride_0_1 = in1->size[0] * in1->size[1];
  in1->size[1] = loop_ub;
  emxEnsureCapacity_boolean_T(sp, in1, stride_0_1, &vb_emlrtRTEI);
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

real_T pfqn_ca(const emlrtStack *sp, emxArray_real_T *L,
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
  const real_T *N_data;
  real_T Gn;
  real_T *G_data;
  real_T *L_data;
  real_T *b_N_data;
  real_T *n_data;
  int32_T M;
  int32_T R;
  int32_T b_r;
  int32_T c_r;
  int32_T d_r;
  int32_T i;
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
  emxInit_real_T(sp, &G, 2, &sb_emlrtRTEI);
  emxInit_real_T(sp, &n, 2, &tb_emlrtRTEI);
  emxInit_boolean_T(sp, &x, &vb_emlrtRTEI);
  emxInit_real_T(sp, &b_x, 2, &xb_emlrtRTEI);
  emxInit_real_T(sp, &b_N, 2, &rb_emlrtRTEI);
  if (L->size[0] == 0) {
    int32_T idx;
    int32_T nz;
    int32_T scalarLB;
    idx = b_N->size[0] * b_N->size[1];
    b_N->size[0] = 1;
    b_N->size[1] = L->size[1];
    emxEnsureCapacity_real_T(sp, b_N, idx, &pb_emlrtRTEI);
    b_N_data = b_N->data;
    for (i = 0; i < R; i++) {
      b_N_data[i] = 0.0;
    }
    st.site = &id_emlrtRSI;
    c_sum(&st, b_N, n);
    st.site = &id_emlrtRSI;
    b_log(&st, n);
    n_data = n->data;
    nz = N->size[1];
    if ((N->size[1] != n->size[1]) &&
        ((N->size[1] != 1) && (n->size[1] != 1))) {
      emlrtDimSizeImpxCheckR2021b(N->size[1], n->size[1], &d_emlrtECI,
                                  (emlrtConstCTX)sp);
    }
    st.site = &id_emlrtRSI;
    /*  lf=FACTLN(n) */
    /*  Compure the logarithm of n!        */
    /*  */
    /*  Copyright (c) 2012-2026, Imperial College London */
    /*  All rights reserved.   */
    b_st.site = &fb_emlrtRSI;
    idx = b_x->size[0] * b_x->size[1];
    b_x->size[0] = 1;
    b_x->size[1] = N->size[1];
    emxEnsureCapacity_real_T(&b_st, b_x, idx, &q_emlrtRTEI);
    L_data = b_x->data;
    scalarLB = (N->size[1] / 2) << 1;
    idx = scalarLB - 2;
    for (i = 0; i <= idx; i += 2) {
      _mm_storeu_pd(&L_data[i],
                    _mm_add_pd(_mm_loadu_pd(&N_data[i]), _mm_set1_pd(1.0)));
    }
    for (i = scalarLB; i < nz; i++) {
      L_data[i] = N_data[i] + 1.0;
    }
    c_st.site = &gb_emlrtRSI;
    applyScalarFunctionInPlace(&c_st, b_x);
    if (N->size[1] == n->size[1]) {
      idx = b_N->size[0] * b_N->size[1];
      b_N->size[0] = 1;
      b_N->size[1] = N->size[1];
      emxEnsureCapacity_real_T(sp, b_N, idx, &qb_emlrtRTEI);
      b_N_data = b_N->data;
      idx = scalarLB - 2;
      for (i = 0; i <= idx; i += 2) {
        __m128d r;
        r = _mm_loadu_pd(&n_data[i]);
        _mm_storeu_pd(&b_N_data[i], _mm_mul_pd(_mm_loadu_pd(&N_data[i]), r));
      }
      for (i = scalarLB; i < nz; i++) {
        b_N_data[i] = N_data[i] * n_data[i];
      }
      st.site = &id_emlrtRSI;
      Gn = -b_sum(&st, b_x) + b_sum(&st, b_N);
    } else {
      st.site = &id_emlrtRSI;
      Gn = binary_expand_op(&st, id_emlrtRSI, b_x, N, n);
    }
    Gn = muDoubleScalarExp(Gn);
  } else {
    int32_T idx;
    int32_T nz;
    boolean_T exitg1;
    st.site = &hd_emlrtRSI;
    b_st.site = &nc_emlrtRSI;
    c_st.site = &oc_emlrtRSI;
    d_st.site = &pc_emlrtRSI;
    if (N->size[1] < 1) {
      emlrtErrorWithMessageIdR2018a(
          &d_st, &b_emlrtRTEI, "Coder:toolbox:eml_min_or_max_varDimZero",
          "Coder:toolbox:eml_min_or_max_varDimZero", 0);
    }
    e_st.site = &qc_emlrtRSI;
    f_st.site = &rc_emlrtRSI;
    if (N->size[1] > 2) {
      g_st.site = &ec_emlrtRSI;
      if (!muDoubleScalarIsNaN(N_data[0])) {
        idx = 1;
      } else {
        idx = 0;
        h_st.site = &fc_emlrtRSI;
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
        g_st.site = &dc_emlrtRSI;
        h_st.site = &gc_emlrtRSI;
        if ((idx + 1 <= N->size[1]) && (N->size[1] > 2147483646)) {
          i_st.site = &x_emlrtRSI;
          check_forloop_overflow_error(&i_st);
        }
      }
    }
    st.site = &gd_emlrtRSI;
    if (b_sum(&st, N) == 0.0) {
      Gn = 1.0;
    } else {
      real_T lGest;
      real_T t;
      int32_T loop_ub;
      int32_T scalarLB;
      boolean_T ok;
      /*  Demand scaling, so that lGn stays computable once G(N) leaves the */
      /*  double-precision range. The recursion below runs in linear space, so
       * it */
      /*  overflows to Inf as soon as G(N) > realmax, i.e. log G > 709.78 -- and
       * lGn */
      /*  was then returned as Inf even though log G is perfectly representable.
       * This */
      /*  is the floating-point range problem Reiser and Lavenberg (1980, JACM
       * 27(2), */
      /*  p.319) report for the convolution algorithm, and the scaling remedy is
       * Lam */
      /*  (1982), "Dynamic scaling and growth behavior of queuing network
       * normalization */
      /*  constants". */
      /*  */
      /*  Every state in G(N) carries the same total population sum(N), so
       * dividing all */
      /*  demands and think times by a constant c divides G(N) by exactly
       * c^sum(N): */
      /*    G(N; L/c, Z/c) = G(N; L, Z) / c^sum(N) */
      /*  hence log G = log G_scaled + sum(N) log c, which is exact, not an */
      /*  approximation. */
      /*  */
      /*  c must be chosen to CENTRE log G_scaled near 0, not merely to shrink
       * the */
      /*  demands. G(N) can leave the double range in EITHER direction: Reiser
       * (1981, */
      /*  Perf. Eval. 1:7-18, Sec. 6.1) reports the unnormalized convolution */
      /*  UNDERFLOWING (< 1e-75, losing all significant digits) on his Fig. 5 */
      /*  central-server model for K > 160, and OVERFLOWING (> 1e75) on the same
       * model */
      /*  if the think time is raised. Scaling by max(L), the textbook choice,
       * only */
      /*  shrinks G and so makes the underflow strictly worse -- on Reiser's
       * model it */
      /*  turns a representable G(200) ~ 1e-87 into an underflow to 0. */
      /*  */
      /*  log G is estimated first, from the largest SINGLE-STATE term, which is
       * a */
      /*  lower bound on G(N) and in practice within O(log #states) of log G: */
      /*    all jobs at one queueing station i : sum_r N_r log L_ir */
      /*    all jobs at the delay              : sum_r [N_r log Z_r - log N_r!]
       */
      /*  On Reiser's model at K=200 this gives max(-460.5, -264.1) = -264.1
       * against a */
      /*  true log G of about -263, i.e. accurate enough to place the scaled
       * value well */
      /*  inside the exponent range. */
      /*  */
      /*  c is then the power of two nearest exp(lGest/sum(N)), so that */
      /*  log G_scaled = log G - sum(N) log c is near 0. A power of two matters:
       * it */
      /*  shifts exponents only, so L/c and Z/c stay exactly representable and
       * the */
      /*  recursion is bit-for-bit the unscaled one with shifted exponents. */
      st.site = &fd_emlrtRSI;
      Gn = b_sum(&st, N);
      lGest = rtMinusInf;
      for (i = 0; i < M; i++) {
        t = 0.0;
        ok = true;
        idx = 1;
        exitg1 = false;
        while ((!exitg1) && (idx - 1 <= R - 1)) {
          real_T d;
          if ((idx < 1) || (idx > N->size[1])) {
            emlrtDynamicBoundsCheckR2012b(idx, 1, N->size[1], &q_emlrtBCI,
                                          (emlrtConstCTX)sp);
          }
          d = N_data[idx - 1];
          if (d > 0.0) {
            real_T c_x;
            if (i + 1 > M) {
              emlrtDynamicBoundsCheckR2012b(i + 1, 1, M, &r_emlrtBCI,
                                            (emlrtConstCTX)sp);
            }
            if (idx > R) {
              emlrtDynamicBoundsCheckR2012b(idx, 1, R, &s_emlrtBCI,
                                            (emlrtConstCTX)sp);
            }
            c_x = L_data[i + L->size[0] * (idx - 1)];
            if (c_x > 0.0) {
              st.site = &ed_emlrtRSI;
              if (i + 1 > M) {
                emlrtDynamicBoundsCheckR2012b(i + 1, 1, M, &t_emlrtBCI, &st);
              }
              if (idx > R) {
                emlrtDynamicBoundsCheckR2012b(idx, 1, R, &u_emlrtBCI, &st);
              }
              if (idx > N->size[1]) {
                emlrtDynamicBoundsCheckR2012b(idx, 1, N->size[1], &v_emlrtBCI,
                                              (emlrtConstCTX)sp);
              }
              t += d * muDoubleScalarLog(c_x);
              idx++;
            } else {
              ok = false;
              exitg1 = true;
            }
          } else {
            idx++;
          }
          if (*emlrtBreakCheckR2012bFlagVar != 0) {
            emlrtBreakCheckR2012b((emlrtConstCTX)sp);
          }
        }
        if (ok) {
          lGest = muDoubleScalarMax(lGest, t);
        }
        if (*emlrtBreakCheckR2012bFlagVar != 0) {
          emlrtBreakCheckR2012b((emlrtConstCTX)sp);
        }
      }
      st.site = &dd_emlrtRSI;
      b_st.site = &jd_emlrtRSI;
      c_st.site = &kd_emlrtRSI;
      if (L->size[1] > 2147483646) {
        d_st.site = &x_emlrtRSI;
        check_forloop_overflow_error(&d_st);
      }
      if (muDoubleScalarIsInf(lGest)) {
        t = 0.0;
      } else {
        st.site = &cd_emlrtRSI;
        t = muDoubleScalarRound(lGest / (Gn * 0.69314718055994529));
      }
      st.site = &bd_emlrtRSI;
      b_st.site = &ld_emlrtRSI;
      c_st.site = &md_emlrtRSI;
      d_st.site = &nd_emlrtRSI;
      e_st.site = &od_emlrtRSI;
      Gn = muDoubleScalarPower(2.0, t);
      idx = L->size[0] * L->size[1];
      nz = (idx / 2) << 1;
      scalarLB = nz - 2;
      for (i = 0; i <= scalarLB; i += 2) {
        __m128d r;
        r = _mm_loadu_pd(&L_data[i]);
        _mm_storeu_pd(&L_data[i], _mm_div_pd(r, _mm_set1_pd(Gn)));
      }
      for (i = nz; i < idx; i++) {
        L_data[i] /= Gn;
      }
      if ((real_T)L->size[0] + 1.0 != L->size[0] + 1) {
        emlrtIntegerCheckR2012b((real_T)L->size[0] + 1.0, &b_emlrtDCI,
                                (emlrtConstCTX)sp);
      }
      idx = b_N->size[0] * b_N->size[1];
      b_N->size[0] = 1;
      loop_ub = N->size[1];
      b_N->size[1] = N->size[1];
      emxEnsureCapacity_real_T(sp, b_N, idx, &rb_emlrtRTEI);
      b_N_data = b_N->data;
      idx = (N->size[1] / 2) << 1;
      nz = idx - 2;
      for (i = 0; i <= nz; i += 2) {
        _mm_storeu_pd(&b_N_data[i],
                      _mm_add_pd(_mm_loadu_pd(&N_data[i]), _mm_set1_pd(1.0)));
      }
      for (i = idx; i < loop_ub; i++) {
        b_N_data[i] = N_data[i] + 1.0;
      }
      st.site = &ad_emlrtRSI;
      Gn = prod(&st, b_N);
      if (Gn != (int32_T)muDoubleScalarFloor(Gn)) {
        emlrtIntegerCheckR2012b(Gn, &emlrtDCI, (emlrtConstCTX)sp);
      }
      idx = G->size[0] * G->size[1];
      G->size[0] = L->size[0] + 1;
      G->size[1] = (int32_T)Gn;
      emxEnsureCapacity_real_T(sp, G, idx, &sb_emlrtRTEI);
      G_data = G->data;
      idx = (L->size[0] + 1) * (int32_T)Gn;
      for (i = 0; i < idx; i++) {
        G_data[i] = 1.0;
      }
      /*  stores G across recursion */
      /*  [N]=PPROD(N,N) */
      /*  sequentially generate all vectors n: 0<=n<=N */
      /*  n=pprod(N) - init */
      /*  n=pprod(n,N) - next state */
      idx = n->size[0] * n->size[1];
      n->size[0] = 1;
      n->size[1] = N->size[1];
      emxEnsureCapacity_real_T(sp, n, idx, &tb_emlrtRTEI);
      n_data = n->data;
      idx = N->size[1];
      for (i = 0; i < idx; i++) {
        n_data[i] = 0.0;
      }
      int32_T exitg11;
      do {
        exitg11 = 0;
        st.site = &yc_emlrtRSI;
        if (b_sum(&st, n) != -1.0) {
          st.site = &xc_emlrtRSI;
          /*  IDX=HASHPOP(N,N,R,PRODS) */
          /*  hash a population vector in n: 0<=n<=N */
          lGest = 1.0;
          for (b_r = 0; b_r < loop_ub; b_r++) {
            if (b_r < 1) {
              scalarLB = 0;
            } else {
              if (b_r > loop_ub) {
                emlrtDynamicBoundsCheckR2012b(b_r, 1, loop_ub, &p_emlrtBCI,
                                              &st);
              }
              scalarLB = b_r;
            }
            idx = b_N->size[0] * b_N->size[1];
            b_N->size[0] = 1;
            b_N->size[1] = scalarLB;
            emxEnsureCapacity_real_T(&st, b_N, idx, &ub_emlrtRTEI);
            b_N_data = b_N->data;
            idx = (scalarLB / 2) << 1;
            nz = idx - 2;
            for (i = 0; i <= nz; i += 2) {
              _mm_storeu_pd(&b_N_data[i], _mm_add_pd(_mm_loadu_pd(&N_data[i]),
                                                     _mm_set1_pd(1.0)));
            }
            for (i = idx; i < scalarLB; i++) {
              b_N_data[i] = N_data[i] + 1.0;
            }
            if (b_r + 1 > n->size[1]) {
              emlrtDynamicBoundsCheckR2012b(b_r + 1, 1, n->size[1], &x_emlrtBCI,
                                            &st);
            }
            b_st.site = &sd_emlrtRSI;
            lGest += prod(&b_st, b_N) * n_data[b_r];
            if (*emlrtBreakCheckR2012bFlagVar != 0) {
              emlrtBreakCheckR2012b(&st);
            }
          }
          st.site = &wc_emlrtRSI;
          /*  F=FZ(Z,N) */
          b_st.site = &td_emlrtRSI;
          if (b_sum(&b_st, n) == 0.0) {
            if (lGest != (int32_T)muDoubleScalarFloor(lGest)) {
              emlrtIntegerCheckR2012b(lGest, &c_emlrtDCI, &st);
            }
            if (((int32_T)lGest < 1) || ((int32_T)lGest > G->size[1])) {
              emlrtDynamicBoundsCheckR2012b((int32_T)lGest, 1, G->size[1],
                                            &w_emlrtBCI, &st);
            }
            G_data[G->size[0] * ((int32_T)lGest - 1)] = 1.0;
          } else {
            idx = 0;
            int32_T exitg2;
            do {
              exitg2 = 0;
              if (idx <= n->size[1] - 1) {
                if (idx + 1 > R) {
                  emlrtDynamicBoundsCheckR2012b(idx + 1, 1, R, &o_emlrtBCI,
                                                &st);
                }
                if (idx + 1 > n->size[1]) {
                  emlrtDynamicBoundsCheckR2012b(idx + 1, 1, n->size[1],
                                                &y_emlrtBCI, &st);
                }
                if (n_data[idx] > 0.0) {
                  if (lGest != (int32_T)muDoubleScalarFloor(lGest)) {
                    emlrtIntegerCheckR2012b(lGest, &c_emlrtDCI, &st);
                  }
                  if (((int32_T)lGest < 1) || ((int32_T)lGest > G->size[1])) {
                    emlrtDynamicBoundsCheckR2012b((int32_T)lGest, 1, G->size[1],
                                                  &w_emlrtBCI, &st);
                  }
                  G_data[G->size[0] * ((int32_T)lGest - 1)] = 0.0;
                  exitg2 = 1;
                } else {
                  idx++;
                }
              } else {
                if (lGest != (int32_T)muDoubleScalarFloor(lGest)) {
                  emlrtIntegerCheckR2012b(lGest, &c_emlrtDCI, &st);
                }
                if (((int32_T)lGest < 1) || ((int32_T)lGest > G->size[1])) {
                  emlrtDynamicBoundsCheckR2012b((int32_T)lGest, 1, G->size[1],
                                                &w_emlrtBCI, &st);
                }
                G_data[G->size[0] * ((int32_T)lGest - 1)] = 1.0;
                exitg2 = 1;
              }
              if (*emlrtBreakCheckR2012bFlagVar != 0) {
                emlrtBreakCheckR2012b(&st);
              }
            } while (exitg2 == 0);
          }
          emlrtForLoopVectorCheckR2021a(2.0, 1.0, (real_T)M + 1.0,
                                        mxDOUBLE_CLASS, M, &k_emlrtRTEI,
                                        (emlrtConstCTX)sp);
          for (b_r = 0; b_r < M; b_r++) {
            if ((int32_T)((uint32_T)b_r + 1U) > G->size[0]) {
              emlrtDynamicBoundsCheckR2012b((int32_T)((uint32_T)b_r + 1U), 1,
                                            G->size[0], &ab_emlrtBCI,
                                            (emlrtConstCTX)sp);
            }
            if (((int32_T)lGest < 1) || ((int32_T)lGest > G->size[1])) {
              emlrtDynamicBoundsCheckR2012b((int32_T)lGest, 1, G->size[1],
                                            &bb_emlrtBCI, (emlrtConstCTX)sp);
            }
            if (((int32_T)((uint32_T)b_r + 2U) < 1) ||
                ((int32_T)((uint32_T)b_r + 2U) > G->size[0])) {
              emlrtDynamicBoundsCheckR2012b((int32_T)((uint32_T)b_r + 2U), 1,
                                            G->size[0], &cb_emlrtBCI,
                                            (emlrtConstCTX)sp);
            }
            if (((int32_T)lGest < 1) || ((int32_T)lGest > G->size[1])) {
              emlrtDynamicBoundsCheckR2012b((int32_T)lGest, 1, G->size[1],
                                            &db_emlrtBCI, (emlrtConstCTX)sp);
            }
            G_data[(b_r + G->size[0] * ((int32_T)lGest - 1)) + 1] =
                G_data[b_r + G->size[0] * ((int32_T)lGest - 1)];
            /*  norm constant with m-1 queues */
            for (c_r = 0; c_r < R; c_r++) {
              if (c_r + 1 > n->size[1]) {
                emlrtDynamicBoundsCheckR2012b(c_r + 1, 1, n->size[1],
                                              &eb_emlrtBCI, (emlrtConstCTX)sp);
              }
              Gn = n_data[c_r];
              if (Gn >= 1.0) {
                if (c_r + 1 > n->size[1]) {
                  emlrtDynamicBoundsCheckR2012b(
                      c_r + 1, 1, n->size[1], &fb_emlrtBCI, (emlrtConstCTX)sp);
                }
                n_data[c_r] = Gn - 1.0;
                st.site = &vc_emlrtRSI;
                /*  IDX=HASHPOP(N,N,R,PRODS) */
                /*  hash a population vector in n: 0<=n<=N */
                Gn = 1.0;
                for (d_r = 0; d_r < loop_ub; d_r++) {
                  if (d_r < 1) {
                    scalarLB = 0;
                  } else {
                    if (d_r > loop_ub) {
                      emlrtDynamicBoundsCheckR2012b(d_r, 1, loop_ub,
                                                    &p_emlrtBCI, &st);
                    }
                    scalarLB = d_r;
                  }
                  idx = b_N->size[0] * b_N->size[1];
                  b_N->size[0] = 1;
                  b_N->size[1] = scalarLB;
                  emxEnsureCapacity_real_T(&st, b_N, idx, &ub_emlrtRTEI);
                  b_N_data = b_N->data;
                  idx = (scalarLB / 2) << 1;
                  nz = idx - 2;
                  for (i = 0; i <= nz; i += 2) {
                    _mm_storeu_pd(
                        &b_N_data[i],
                        _mm_add_pd(_mm_loadu_pd(&N_data[i]), _mm_set1_pd(1.0)));
                  }
                  for (i = idx; i < scalarLB; i++) {
                    b_N_data[i] = N_data[i] + 1.0;
                  }
                  if (d_r + 1 > n->size[1]) {
                    emlrtDynamicBoundsCheckR2012b(d_r + 1, 1, n->size[1],
                                                  &x_emlrtBCI, &st);
                  }
                  b_st.site = &sd_emlrtRSI;
                  Gn += prod(&b_st, b_N) * n_data[d_r];
                  if (*emlrtBreakCheckR2012bFlagVar != 0) {
                    emlrtBreakCheckR2012b(&st);
                  }
                }
                if (c_r + 1 > n->size[1]) {
                  emlrtDynamicBoundsCheckR2012b(
                      c_r + 1, 1, n->size[1], &gb_emlrtBCI, (emlrtConstCTX)sp);
                }
                n_data[c_r]++;
                if (((int32_T)((uint32_T)b_r + 2U) < 1) ||
                    ((int32_T)((uint32_T)b_r + 2U) > G->size[0])) {
                  emlrtDynamicBoundsCheckR2012b((int32_T)((uint32_T)b_r + 2U),
                                                1, G->size[0], &hb_emlrtBCI,
                                                (emlrtConstCTX)sp);
                }
                if (((int32_T)lGest < 1) || ((int32_T)lGest > G->size[1])) {
                  emlrtDynamicBoundsCheckR2012b((int32_T)lGest, 1, G->size[1],
                                                &ib_emlrtBCI,
                                                (emlrtConstCTX)sp);
                }
                if ((int32_T)((uint32_T)b_r + 1U) > L->size[0]) {
                  emlrtDynamicBoundsCheckR2012b((int32_T)((uint32_T)b_r + 1U),
                                                1, L->size[0], &jb_emlrtBCI,
                                                (emlrtConstCTX)sp);
                }
                if (c_r + 1 > L->size[1]) {
                  emlrtDynamicBoundsCheckR2012b(
                      c_r + 1, 1, L->size[1], &kb_emlrtBCI, (emlrtConstCTX)sp);
                }
                if (((int32_T)((uint32_T)b_r + 2U) < 1) ||
                    ((int32_T)((uint32_T)b_r + 2U) > G->size[0])) {
                  emlrtDynamicBoundsCheckR2012b((int32_T)((uint32_T)b_r + 2U),
                                                1, G->size[0], &lb_emlrtBCI,
                                                (emlrtConstCTX)sp);
                }
                if (Gn != (int32_T)muDoubleScalarFloor(Gn)) {
                  emlrtIntegerCheckR2012b(Gn, &d_emlrtDCI, (emlrtConstCTX)sp);
                }
                if (((int32_T)Gn < 1) || ((int32_T)Gn > G->size[1])) {
                  emlrtDynamicBoundsCheckR2012b((int32_T)Gn, 1, G->size[1],
                                                &mb_emlrtBCI,
                                                (emlrtConstCTX)sp);
                }
                if (((int32_T)((uint32_T)b_r + 2U) < 1) ||
                    ((int32_T)((uint32_T)b_r + 2U) > G->size[0])) {
                  emlrtDynamicBoundsCheckR2012b((int32_T)((uint32_T)b_r + 2U),
                                                1, G->size[0], &nb_emlrtBCI,
                                                (emlrtConstCTX)sp);
                }
                if (((int32_T)lGest < 1) || ((int32_T)lGest > G->size[1])) {
                  emlrtDynamicBoundsCheckR2012b((int32_T)lGest, 1, G->size[1],
                                                &ob_emlrtBCI,
                                                (emlrtConstCTX)sp);
                }
                G_data[(b_r + G->size[0] * ((int32_T)lGest - 1)) + 1] +=
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
          st.site = &uc_emlrtRSI;
          /*  [N]=PPROD(N,N) */
          /*  sequentially generate all vectors n: 0<=n<=N */
          /*  n=pprod(N) - init */
          /*  n=pprod(n,N) - next state */
          nz = n->size[1];
          if ((n->size[1] != loop_ub) &&
              ((n->size[1] != 1) && (loop_ub != 1))) {
            emlrtDimSizeImpxCheckR2021b(n->size[1], loop_ub, &e_emlrtECI, &st);
          }
          b_st.site = &ud_emlrtRSI;
          if (n->size[1] == N->size[1]) {
            idx = x->size[0] * x->size[1];
            x->size[0] = 1;
            x->size[1] = n->size[1];
            emxEnsureCapacity_boolean_T(&b_st, x, idx, &vb_emlrtRTEI);
            x_data = x->data;
            for (i = 0; i < nz; i++) {
              x_data[i] = (n_data[i] == N_data[i]);
            }
          } else {
            c_st.site = &ud_emlrtRSI;
            eq(&c_st, x, n, N);
            x_data = x->data;
          }
          c_st.site = &vd_emlrtRSI;
          d_st.site = &o_emlrtRSI;
          idx = x->size[1];
          if (x->size[1] == 0) {
            nz = 0;
          } else {
            e_st.site = &qd_emlrtRSI;
            nz = x_data[0];
            f_st.site = &rd_emlrtRSI;
            if (x->size[1] > 2147483646) {
              g_st.site = &x_emlrtRSI;
              check_forloop_overflow_error(&g_st);
            }
            for (i = 2; i <= idx; i++) {
              nz += x_data[i - 1];
            }
          }
          if (nz == N->size[1]) {
            idx = n->size[0] * n->size[1];
            n->size[0] = 1;
            n->size[1] = 1;
            emxEnsureCapacity_real_T(&st, n, idx, &wb_emlrtRTEI);
            n_data = n->data;
            n_data[0] = -1.0;
          } else {
            idx = N->size[1];
            exitg1 = false;
            while ((!exitg1) && (idx > 0)) {
              if (idx > n->size[1]) {
                emlrtDynamicBoundsCheckR2012b(idx, 1, n->size[1], &pb_emlrtBCI,
                                              &st);
              }
              if (idx > loop_ub) {
                emlrtDynamicBoundsCheckR2012b(idx, 1, loop_ub, &qb_emlrtBCI,
                                              &st);
              }
              if (n_data[idx - 1] == N_data[idx - 1]) {
                if (idx > n->size[1]) {
                  emlrtDynamicBoundsCheckR2012b(idx, 1, n->size[1],
                                                &rb_emlrtBCI, &st);
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
                emlrtDynamicBoundsCheckR2012b(idx, 1, n->size[1], &sb_emlrtBCI,
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
      /*  Undo the scaling in log space: log G = log G_scaled + sum(N) log c.
       * lGn is */
      /*  therefore finite whenever log G itself is, even though Gn may
       * legitimately */
      /*  overflow to Inf (the true constant really is outside double range). */
      st.site = &tc_emlrtRSI;
      ok = ((L->size[0] + 1 < 1) || (L->size[0] + 1 > G->size[0]));
      if (ok) {
        emlrtDynamicBoundsCheckR2012b(L->size[0] + 1, 1, G->size[0],
                                      &n_emlrtBCI, &st);
      }
      if (G->size[1] < 1) {
        emlrtDynamicBoundsCheckR2012b(G->size[1], 1, G->size[1], &m_emlrtBCI,
                                      &st);
      }
      lGest = G_data[L->size[0] + G->size[0] * (G->size[1] - 1)];
      if (lGest < 0.0) {
        emlrtErrorWithMessageIdR2018a(
            &st, &e_emlrtRTEI, "Coder:toolbox:ElFunDomainError",
            "Coder:toolbox:ElFunDomainError", 3, 4, 3, "log");
      }
      st.site = &tc_emlrtRSI;
      b_sum(&st, N);
      st.site = &tc_emlrtRSI;
      /*  Gn is recovered by scaling the mantissa back with pow2 (an exact
       * exponent */
      /*  adjustment), not as exp(lGn): on models that never overflowed this
       * returns */
      /*  bit-for-bit the value the unscaled recursion used to return, and it
       * still */
      /*  goes to Inf when the constant genuinely leaves double range -- in
       * which case */
      /*  lGn above remains finite and usable. */
      if (G->size[1] < 1) {
        emlrtDynamicBoundsCheckR2012b(G->size[1], 1, G->size[1], &l_emlrtBCI,
                                      (emlrtConstCTX)sp);
      }
      st.site = &sc_emlrtRSI;
      Gn = b_sum(&st, N) * t;
      if (Gn < 0.0) {
        Gn = muDoubleScalarCeil(Gn);
        if (Gn < -32768.0) {
          Gn = -32768.0;
        }
      } else {
        Gn = muDoubleScalarFloor(Gn);
        if (Gn > 32767.0) {
          Gn = 32767.0;
        }
      }
      Gn = ldexp(lGest, (int32_T)Gn);
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
