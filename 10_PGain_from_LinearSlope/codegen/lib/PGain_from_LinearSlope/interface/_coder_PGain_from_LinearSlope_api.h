/*
 * File: _coder_PGain_from_LinearSlope_api.h
 *
 * MATLAB Coder version            : 3.1
 * C/C++ source code generated on  : 20-Dec-2017 17:38:07
 */

#ifndef _CODER_PGAIN_FROM_LINEARSLOPE_API_H
#define _CODER_PGAIN_FROM_LINEARSLOPE_API_H

/* Include Files */
#include "tmwtypes.h"
#include "mex.h"
#include "emlrt.h"
#include <stddef.h>
#include <stdlib.h>
#include "_coder_PGain_from_LinearSlope_api.h"

/* Variable Declarations */
extern emlrtCTX emlrtRootTLSGlobal;
extern emlrtContext emlrtContextGlobal;

/* Function Declarations */
extern void PGain_from_LinearSlope(real_T inf[1000], real_T r[1000], real_T n,
  real_T MinFreqRange, real_T MaxFreqRange, real_T Slope, real_T Tolerance,
  real_T vPh[1000], real_T *StartFrequency, real_T *EndFrequency, real_T
  *SlopeResult, real_T *NumberOfPoints, real_T *PlantGain, real_T *nErrCode);
extern void PGain_from_LinearSlope_api(const mxArray *prhs[8], const mxArray
  *plhs[6]);
extern void PGain_from_LinearSlope_atexit(void);
extern void PGain_from_LinearSlope_initialize(void);
extern void PGain_from_LinearSlope_terminate(void);
extern void PGain_from_LinearSlope_xil_terminate(void);

#endif

/*
 * File trailer for _coder_PGain_from_LinearSlope_api.h
 *
 * [EOF]
 */
