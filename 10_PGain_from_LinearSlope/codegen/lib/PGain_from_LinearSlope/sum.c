/*
 * File: sum.c
 *
 * MATLAB Coder version            : 3.1
 * C/C++ source code generated on  : 20-Dec-2017 17:38:07
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "PGain_from_LinearSlope.h"
#include "sum.h"

/* Function Definitions */

/*
 * Arguments    : const double x[10000000]
 * Return Type  : double
 */
double sum(const double x[10000000])
{
  double y;
  int k;
  y = x[0];
  for (k = 0; k < 9999999; k++) {
    y += x[k + 1];
  }

  return y;
}

/*
 * File trailer for sum.c
 *
 * [EOF]
 */
