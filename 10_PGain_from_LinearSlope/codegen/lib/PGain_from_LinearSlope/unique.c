/*
 * File: unique.c
 *
 * MATLAB Coder version            : 3.1
 * C/C++ source code generated on  : 20-Dec-2017 17:38:07
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "PGain_from_LinearSlope.h"
#include "unique.h"

/* Function Definitions */

/*
 * Arguments    : const double b_data[]
 *                int nb
 *                int *nMInf
 *                int *nFinite
 *                int *nPInf
 *                int *nNaN
 * Return Type  : void
 */
void count_nonfinites(const double b_data[], int nb, int *nMInf, int *nFinite,
                      int *nPInf, int *nNaN)
{
  int k;
  k = 0;
  while ((k + 1 <= nb) && rtIsInf(b_data[k]) && (b_data[k] < 0.0)) {
    k++;
  }

  *nMInf = k;
  k = nb;
  while ((k >= 1) && rtIsNaN(b_data[k - 1])) {
    k--;
  }

  *nNaN = nb - k;
  while ((k >= 1) && rtIsInf(b_data[k - 1]) && (b_data[k - 1] > 0.0)) {
    k--;
  }

  *nPInf = (nb - k) - *nNaN;
  *nFinite = k - *nMInf;
}

/*
 * File trailer for unique.c
 *
 * [EOF]
 */
