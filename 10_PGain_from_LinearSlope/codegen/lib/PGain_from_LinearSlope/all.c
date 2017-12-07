/*
 * File: all.c
 *
 * MATLAB Coder version            : 3.1
 * C/C++ source code generated on  : 20-Dec-2017 17:38:07
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "PGain_from_LinearSlope.h"
#include "all.h"

/* Function Definitions */

/*
 * Arguments    : const boolean_T x_data[]
 * Return Type  : boolean_T
 */
boolean_T all(const boolean_T x_data[])
{
  boolean_T y;
  int ix;
  boolean_T exitg1;
  y = true;
  ix = 1;
  exitg1 = false;
  while ((!exitg1) && (ix <= 1)) {
    if (!x_data[0]) {
      y = false;
      exitg1 = true;
    } else {
      ix = 2;
    }
  }

  return y;
}

/*
 * File trailer for all.c
 *
 * [EOF]
 */
