/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * pfqn_panacea.h
 *
 * Code generation for function 'pfqn_panacea'
 *
 */

#pragma once

/* Include files */
#include "pfqn_panacea_types.h"
#include "rtwtypes.h"
#include "emlrt.h"
#include "mex.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Function Declarations */
real_T binary_expand_op(const emlrtStack *sp, const emlrtRSInfo in1,
                        const emxArray_real_T *in2, const emxArray_real_T *in3,
                        const emxArray_real_T *in4);

emlrtCTX emlrtGetRootTLSGlobal(void);

void emlrtLockerFunction(EmlrtLockeeFunction aLockee, emlrtConstCTX aTLS,
                         void *aData);

void pfqn_panacea(const emlrtStack *sp, const emxArray_real_T *L,
                  const emxArray_real_T *N, emxArray_real_T *Z, real_T *Gn,
                  real_T *lGn);

/* End of code generation (pfqn_panacea.h) */
