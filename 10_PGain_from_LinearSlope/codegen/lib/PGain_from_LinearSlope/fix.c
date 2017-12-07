/*
 * File: fix.c
 *
 * MATLAB Coder version            : 3.1
 * C/C++ source code generated on  : 20-Dec-2017 17:38:07
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "PGain_from_LinearSlope.h"
#include "fix.h"

/* Function Definitions */

/*
 * Arguments    : double *x
 * Return Type  : void
 */
void b_fix(double *x)
{
  if (*x < 0.0) {
    *x = ceil(*x);
  } else {
    *x = floor(*x);
  }
}

/*
 * File trailer for fix.c
 *
 * [EOF]
 */
