/*
 * File: _coder_LinearSlope4Manual_api.h
 *
 * MATLAB Coder version            : 3.1
 * C/C++ source code generated on  : 20-Dec-2017 18:14:33
 */

#ifndef _CODER_LINEARSLOPE4MANUAL_API_H
#define _CODER_LINEARSLOPE4MANUAL_API_H

/* Include Files */
#include "tmwtypes.h"
#include "mex.h"
#include "emlrt.h"
#include <stddef.h>
#include <stdlib.h>
#include "_coder_LinearSlope4Manual_api.h"

/* Variable Declarations */
extern emlrtCTX emlrtRootTLSGlobal;
extern emlrtContext emlrtContextGlobal;

/* Function Declarations */
extern void LinearSlope4Manual(real_T f[1000], real_T r[1000], real_T n, real_T
  MinFreqRange, real_T MaxFreqRange, real_T *SlopeResult, real_T *PlantGain,
  real_T *ErrCode);
extern void LinearSlope4Manual_api(const mxArray *prhs[5], const mxArray *plhs[3]);
extern void LinearSlope4Manual_atexit(void);
extern void LinearSlope4Manual_initialize(void);
extern void LinearSlope4Manual_terminate(void);
extern void LinearSlope4Manual_xil_terminate(void);

#endif

/*
 * File trailer for _coder_LinearSlope4Manual_api.h
 *
 * [EOF]
 */
