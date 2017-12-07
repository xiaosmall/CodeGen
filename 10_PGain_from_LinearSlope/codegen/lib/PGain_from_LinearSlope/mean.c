/*
 * File: mean.c
 *
 * MATLAB Coder version            : 3.1
 * C/C++ source code generated on  : 20-Dec-2017 17:38:07
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "PGain_from_LinearSlope.h"
#include "mean.h"

/* Function Definitions */

/*
 * Arguments    : const emxArray_real_T *x
 * Return Type  : double
 */
double mean(const emxArray_real_T *x)
{
  double y;
  int k;
  y = x->data[0];
  for (k = 2; k <= x->size[0]; k++) {
    y += x->data[k - 1];
  }

  y /= (double)x->size[0];
  return y;
}

/*
 * File trailer for mean.c
 *
 * [EOF]
 */
