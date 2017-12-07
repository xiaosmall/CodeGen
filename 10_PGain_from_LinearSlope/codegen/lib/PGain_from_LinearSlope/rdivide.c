/*
 * File: rdivide.c
 *
 * MATLAB Coder version            : 3.1
 * C/C++ source code generated on  : 20-Dec-2017 17:38:07
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "PGain_from_LinearSlope.h"
#include "rdivide.h"
#include "PGain_from_LinearSlope_emxutil.h"

/* Function Definitions */

/*
 * Arguments    : const emxArray_real_T *x
 *                const emxArray_real_T *y
 *                emxArray_real_T *z
 * Return Type  : void
 */
void rdivide(const emxArray_real_T *x, const emxArray_real_T *y, emxArray_real_T
             *z)
{
  int i2;
  int loop_ub;
  i2 = z->size[0];
  z->size[0] = x->size[0];
  emxEnsureCapacity((emxArray__common *)z, i2, (int)sizeof(double));
  loop_ub = x->size[0];
  for (i2 = 0; i2 < loop_ub; i2++) {
    z->data[i2] = x->data[i2] / y->data[i2];
  }
}

/*
 * File trailer for rdivide.c
 *
 * [EOF]
 */
