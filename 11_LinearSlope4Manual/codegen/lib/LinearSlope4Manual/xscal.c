/*
 * File: xscal.c
 *
 * MATLAB Coder version            : 3.1
 * C/C++ source code generated on  : 20-Dec-2017 18:14:33
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "LinearSlope4Manual.h"
#include "xscal.h"

/* Function Definitions */

/*
 * Arguments    : int n
 *                double a
 *                emxArray_real_T *x
 *                int ix0
 * Return Type  : void
 */
void xscal(int n, double a, emxArray_real_T *x, int ix0)
{
  int i1;
  int k;
  i1 = (ix0 + n) - 1;
  for (k = ix0; k <= i1; k++) {
    x->data[k - 1] *= a;
  }
}

/*
 * File trailer for xscal.c
 *
 * [EOF]
 */
