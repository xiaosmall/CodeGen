/*
 * File: polyfit.h
 *
 * MATLAB Coder version            : 3.1
 * C/C++ source code generated on  : 20-Dec-2017 18:14:33
 */

#ifndef POLYFIT_H
#define POLYFIT_H

/* Include Files */
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "rt_nonfinite.h"
#include "rtwtypes.h"
#include "LinearSlope4Manual_types.h"

/* Function Declarations */
extern void polyfit(const emxArray_real_T *x, const emxArray_real_T *y, double
                    p[2], double S_R_data[], int S_R_size[2], double *S_df,
                    double *S_normr);

#endif

/*
 * File trailer for polyfit.h
 *
 * [EOF]
 */
