/*
 * File: nullAssignment.c
 *
 * MATLAB Coder version            : 3.1
 * C/C++ source code generated on  : 20-Dec-2017 17:38:07
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "PGain_from_LinearSlope.h"
#include "nullAssignment.h"
#include "PGain_from_LinearSlope_emxutil.h"

/* Function Definitions */

/*
 * Arguments    : emxArray_real_T *x
 *                int idx
 * Return Type  : void
 */
void nullAssignment(emxArray_real_T *x, int idx)
{
  int nxin;
  int nxout;
  int k;
  emxArray_real_T *b_x;
  nxin = x->size[0];
  nxout = x->size[0] - 1;
  for (k = idx; k < nxin; k++) {
    x->data[k - 1] = x->data[k];
  }

  emxInit_real_T1(&b_x, 1);
  k = b_x->size[0];
  b_x->size[0] = nxout;
  emxEnsureCapacity((emxArray__common *)b_x, k, (int)sizeof(double));
  for (k = 0; k < nxout; k++) {
    b_x->data[k] = x->data[k];
  }

  k = x->size[0];
  x->size[0] = b_x->size[0];
  emxEnsureCapacity((emxArray__common *)x, k, (int)sizeof(double));
  nxin = b_x->size[0];
  for (k = 0; k < nxin; k++) {
    x->data[k] = b_x->data[k];
  }

  emxFree_real_T(&b_x);
}

/*
 * File trailer for nullAssignment.c
 *
 * [EOF]
 */
