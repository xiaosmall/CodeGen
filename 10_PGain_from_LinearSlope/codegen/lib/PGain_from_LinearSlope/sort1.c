/*
 * File: sort1.c
 *
 * MATLAB Coder version            : 3.1
 * C/C++ source code generated on  : 20-Dec-2017 17:38:07
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "PGain_from_LinearSlope.h"
#include "sort1.h"
#include "PGain_from_LinearSlope_emxutil.h"
#include "sortIdx.h"

/* Function Definitions */

/*
 * Arguments    : emxArray_real_T *x
 * Return Type  : void
 */
void sort(emxArray_real_T *x)
{
  emxArray_real_T *vwork;
  int i5;
  int k;
  int i6;
  emxArray_int32_T *b_vwork;
  emxInit_real_T1(&vwork, 1);
  i5 = x->size[0];
  k = x->size[0];
  i6 = vwork->size[0];
  vwork->size[0] = k;
  emxEnsureCapacity((emxArray__common *)vwork, i6, (int)sizeof(double));
  emxInit_int32_T(&b_vwork, 1);
  for (k = 0; k + 1 <= i5; k++) {
    vwork->data[k] = x->data[k];
  }

  b_sortIdx(vwork, b_vwork);
  for (k = 0; k + 1 <= i5; k++) {
    x->data[k] = vwork->data[k];
  }

  emxFree_int32_T(&b_vwork);
  emxFree_real_T(&vwork);
}

/*
 * File trailer for sort1.c
 *
 * [EOF]
 */
