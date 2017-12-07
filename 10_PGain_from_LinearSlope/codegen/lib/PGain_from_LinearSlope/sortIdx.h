/*
 * File: sortIdx.h
 *
 * MATLAB Coder version            : 3.1
 * C/C++ source code generated on  : 20-Dec-2017 17:38:07
 */

#ifndef SORTIDX_H
#define SORTIDX_H

/* Include Files */
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "rt_nonfinite.h"
#include "rtwtypes.h"
#include "PGain_from_LinearSlope_types.h"

/* Function Declarations */
extern void b_sortIdx(emxArray_real_T *x, emxArray_int32_T *idx);
extern void sortIdx(const double x_data[], const int x_size[1], int idx_data[],
                    int idx_size[1]);

#endif

/*
 * File trailer for sortIdx.h
 *
 * [EOF]
 */
