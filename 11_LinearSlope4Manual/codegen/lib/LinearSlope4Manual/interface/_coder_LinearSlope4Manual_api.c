/*
 * File: _coder_LinearSlope4Manual_api.c
 *
 * MATLAB Coder version            : 3.1
 * C/C++ source code generated on  : 20-Dec-2017 18:14:33
 */

/* Include Files */
#include "tmwtypes.h"
#include "_coder_LinearSlope4Manual_api.h"
#include "_coder_LinearSlope4Manual_mex.h"

/* Variable Definitions */
emlrtCTX emlrtRootTLSGlobal = NULL;
emlrtContext emlrtContextGlobal = { true, false, 131434U, NULL,
  "LinearSlope4Manual", NULL, false, { 2045744189U, 2170104910U, 2743257031U,
    4284093946U }, NULL };

/* Function Declarations */
static real_T (*b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId))[1000];
static real_T c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *n, const
  char_T *identifier);
static real_T d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId);
static real_T (*e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[1000];
static real_T (*emlrt_marshallIn(const emlrtStack *sp, const mxArray *f, const
  char_T *identifier))[1000];
static const mxArray *emlrt_marshallOut(const real_T u);
static real_T f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId);

/* Function Definitions */

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *u
 *                const emlrtMsgIdentifier *parentId
 * Return Type  : real_T (*)[1000]
 */
static real_T (*b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId))[1000]
{
  real_T (*y)[1000];
  y = e_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}
/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *n
 *                const char_T *identifier
 * Return Type  : real_T
 */
  static real_T c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *n, const
  char_T *identifier)
{
  real_T y;
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = d_emlrt_marshallIn(sp, emlrtAlias(n), &thisId);
  emlrtDestroyArray(&n);
  return y;
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *u
 *                const emlrtMsgIdentifier *parentId
 * Return Type  : real_T
 */
static real_T d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId)
{
  real_T y;
  y = f_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *src
 *                const emlrtMsgIdentifier *msgId
 * Return Type  : real_T (*)[1000]
 */
static real_T (*e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[1000]
{
  real_T (*ret)[1000];
  static const int32_T dims[1] = { 1000 };

  emlrtCheckBuiltInR2012b(sp, msgId, src, "double", false, 1U, dims);
  ret = (real_T (*)[1000])mxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}
/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *f
 *                const char_T *identifier
 * Return Type  : real_T (*)[1000]
 */
  static real_T (*emlrt_marshallIn(const emlrtStack *sp, const mxArray *f, const
  char_T *identifier))[1000]
{
  real_T (*y)[1000];
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = b_emlrt_marshallIn(sp, emlrtAlias(f), &thisId);
  emlrtDestroyArray(&f);
  return y;
}

/*
 * Arguments    : const real_T u
 * Return Type  : const mxArray *
 */
static const mxArray *emlrt_marshallOut(const real_T u)
{
  const mxArray *y;
  const mxArray *m0;
  y = NULL;
  m0 = emlrtCreateDoubleScalar(u);
  emlrtAssign(&y, m0);
  return y;
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *src
 *                const emlrtMsgIdentifier *msgId
 * Return Type  : real_T
 */
static real_T f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId)
{
  real_T ret;
  static const int32_T dims = 0;
  emlrtCheckBuiltInR2012b(sp, msgId, src, "double", false, 0U, &dims);
  ret = *(real_T *)mxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

/*
 * Arguments    : const mxArray *prhs[5]
 *                const mxArray *plhs[3]
 * Return Type  : void
 */
void LinearSlope4Manual_api(const mxArray *prhs[5], const mxArray *plhs[3])
{
  real_T (*f)[1000];
  real_T (*r)[1000];
  real_T n;
  real_T MinFreqRange;
  real_T MaxFreqRange;
  real_T SlopeResult;
  real_T PlantGain;
  real_T ErrCode;
  emlrtStack st = { NULL, NULL, NULL };

  st.tls = emlrtRootTLSGlobal;
  prhs[0] = emlrtProtectR2012b(prhs[0], 0, false, -1);
  prhs[1] = emlrtProtectR2012b(prhs[1], 1, false, -1);

  /* Marshall function inputs */
  f = emlrt_marshallIn(&st, emlrtAlias(prhs[0]), "f");
  r = emlrt_marshallIn(&st, emlrtAlias(prhs[1]), "r");
  n = c_emlrt_marshallIn(&st, emlrtAliasP(prhs[2]), "n");
  MinFreqRange = c_emlrt_marshallIn(&st, emlrtAliasP(prhs[3]), "MinFreqRange");
  MaxFreqRange = c_emlrt_marshallIn(&st, emlrtAliasP(prhs[4]), "MaxFreqRange");

  /* Invoke the target function */
  LinearSlope4Manual(*f, *r, n, MinFreqRange, MaxFreqRange, &SlopeResult,
                     &PlantGain, &ErrCode);

  /* Marshall function outputs */
  plhs[0] = emlrt_marshallOut(SlopeResult);
  plhs[1] = emlrt_marshallOut(PlantGain);
  plhs[2] = emlrt_marshallOut(ErrCode);
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void LinearSlope4Manual_atexit(void)
{
  emlrtStack st = { NULL, NULL, NULL };

  mexFunctionCreateRootTLS();
  st.tls = emlrtRootTLSGlobal;
  emlrtEnterRtStackR2012b(&st);
  emlrtLeaveRtStackR2012b(&st);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
  LinearSlope4Manual_xil_terminate();
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void LinearSlope4Manual_initialize(void)
{
  emlrtStack st = { NULL, NULL, NULL };

  mexFunctionCreateRootTLS();
  st.tls = emlrtRootTLSGlobal;
  emlrtClearAllocCountR2012b(&st, false, 0U, 0);
  emlrtEnterRtStackR2012b(&st);
  emlrtFirstTimeR2012b(emlrtRootTLSGlobal);
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void LinearSlope4Manual_terminate(void)
{
  emlrtStack st = { NULL, NULL, NULL };

  st.tls = emlrtRootTLSGlobal;
  emlrtLeaveRtStackR2012b(&st);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
}

/*
 * File trailer for _coder_LinearSlope4Manual_api.c
 *
 * [EOF]
 */
