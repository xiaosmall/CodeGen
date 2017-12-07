/*
 * File: PGain_from_LinearSlope.h
 *
 * MATLAB Coder version            : 3.1
 * C/C++ source code generated on  : 20-Dec-2017 17:38:07
 */

#ifndef PGAIN_FROM_LINEARSLOPE_H
#define PGAIN_FROM_LINEARSLOPE_H

/* Include Files */
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "rt_nonfinite.h"
#include "rtwtypes.h"
#include "PGain_from_LinearSlope_types.h"

/* Function Declarations */
extern void PGain_from_LinearSlope(const double inf[1000], const double r[1000],
  double n, double MinFreqRange, double MaxFreqRange, double Slope, double
  Tolerance, const double vPh[1000], double *StartFrequency, double
  *EndFrequency, double *SlopeResult, double *NumberOfPoints, double *PlantGain,
  double *nErrCode);

#endif

/*
 * File trailer for PGain_from_LinearSlope.h
 *
 * [EOF]
 */
