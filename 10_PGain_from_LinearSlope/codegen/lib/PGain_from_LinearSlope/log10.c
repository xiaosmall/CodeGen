/*
 * File: log10.c
 *
 * MATLAB Coder version            : 3.1
 * C/C++ source code generated on  : 20-Dec-2017 17:38:07
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "PGain_from_LinearSlope.h"
#include "log10.h"

/* Function Definitions */

/*
 * Arguments    : emxArray_real_T *x
 * Return Type  : void
 */
void b_log10(emxArray_real_T *x)
{
  int nx;
  int k;
  nx = x->size[0];
  for (k = 0; k + 1 <= nx; k++) {
    x->data[k] = log10(x->data[k]);
  }
}

/*
 * File trailer for log10.c
 *
 * [EOF]
 */
