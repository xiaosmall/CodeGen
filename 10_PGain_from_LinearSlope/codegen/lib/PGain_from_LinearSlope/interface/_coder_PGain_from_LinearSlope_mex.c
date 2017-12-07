/*
 * File: _coder_PGain_from_LinearSlope_mex.c
 *
 * MATLAB Coder version            : 3.1
 * C/C++ source code generated on  : 20-Dec-2017 17:38:07
 */

/* Include Files */
#include "_coder_PGain_from_LinearSlope_api.h"
#include "_coder_PGain_from_LinearSlope_mex.h"

/* Function Declarations */
static void c_PGain_from_LinearSlope_mexFun(int32_T nlhs, mxArray *plhs[6],
  int32_T nrhs, const mxArray *prhs[8]);

/* Function Definitions */

/*
 * Arguments    : int32_T nlhs
 *                const mxArray *plhs[6]
 *                int32_T nrhs
 *                const mxArray *prhs[8]
 * Return Type  : void
 */
static void c_PGain_from_LinearSlope_mexFun(int32_T nlhs, mxArray *plhs[6],
  int32_T nrhs, const mxArray *prhs[8])
{
  int32_T n;
  const mxArray *inputs[8];
  const mxArray *outputs[6];
  int32_T b_nlhs;
  emlrtStack st = { NULL, NULL, NULL };

  st.tls = emlrtRootTLSGlobal;

  /* Check for proper number of arguments. */
  if (nrhs != 8) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:WrongNumberOfInputs", 5, 12, 8, 4,
                        22, "PGain_from_LinearSlope");
  }

  if (nlhs > 6) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:TooManyOutputArguments", 3, 4, 22,
                        "PGain_from_LinearSlope");
  }

  /* Temporary copy for mex inputs. */
  for (n = 0; n < nrhs; n++) {
    inputs[n] = prhs[n];
  }

  /* Call the function. */
  PGain_from_LinearSlope_api(inputs, outputs);

  /* Copy over outputs to the caller. */
  if (nlhs < 1) {
    b_nlhs = 1;
  } else {
    b_nlhs = nlhs;
  }

  emlrtReturnArrays(b_nlhs, plhs, outputs);

  /* Module termination. */
  PGain_from_LinearSlope_terminate();
}

/*
 * Arguments    : int32_T nlhs
 *                const mxArray * const plhs[]
 *                int32_T nrhs
 *                const mxArray * const prhs[]
 * Return Type  : void
 */
void mexFunction(int32_T nlhs, mxArray *plhs[], int32_T nrhs, const mxArray
                 *prhs[])
{
  mexAtExit(PGain_from_LinearSlope_atexit);

  /* Initialize the memory manager. */
  /* Module initialization. */
  PGain_from_LinearSlope_initialize();

  /* Dispatch the entry-point. */
  c_PGain_from_LinearSlope_mexFun(nlhs, plhs, nrhs, prhs);
}

/*
 * Arguments    : void
 * Return Type  : emlrtCTX
 */
emlrtCTX mexFunctionCreateRootTLS(void)
{
  emlrtCreateRootTLS(&emlrtRootTLSGlobal, &emlrtContextGlobal, NULL, 1);
  return emlrtRootTLSGlobal;
}

/*
 * File trailer for _coder_PGain_from_LinearSlope_mex.c
 *
 * [EOF]
 */
