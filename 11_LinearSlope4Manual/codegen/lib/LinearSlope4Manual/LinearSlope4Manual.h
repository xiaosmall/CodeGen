/*
 * File: LinearSlope4Manual.h
 *
 * MATLAB Coder version            : 3.1
 * C/C++ source code generated on  : 20-Dec-2017 18:14:33
 */

#ifndef LINEARSLOPE4MANUAL_H
#define LINEARSLOPE4MANUAL_H

/* Include Files */
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "rt_nonfinite.h"
#include "rtwtypes.h"
#include "LinearSlope4Manual_types.h"

/* Function Declarations */
extern void LinearSlope4Manual(const double f[1000], const double r[1000],
  double n, double MinFreqRange, double MaxFreqRange, double *SlopeResult,
  double *PlantGain, double *ErrCode);

#endif

/*
 * File trailer for LinearSlope4Manual.h
 *
 * [EOF]
 */
