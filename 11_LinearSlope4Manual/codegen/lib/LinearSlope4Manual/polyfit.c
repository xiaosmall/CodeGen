/*
 * File: polyfit.c
 *
 * MATLAB Coder version            : 3.1
 * C/C++ source code generated on  : 20-Dec-2017 18:14:33
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "LinearSlope4Manual.h"
#include "polyfit.h"
#include "LinearSlope4Manual_emxutil.h"
#include "xgeqp3.h"

/* Function Definitions */

/*
 * Arguments    : const emxArray_real_T *x
 *                const emxArray_real_T *y
 *                double p[2]
 *                double S_R_data[]
 *                int S_R_size[2]
 *                double *S_df
 *                double *S_normr
 * Return Type  : void
 */
void polyfit(const emxArray_real_T *x, const emxArray_real_T *y, double p[2],
             double S_R_data[], int S_R_size[2], double *S_df, double *S_normr)
{
  emxArray_real_T *V;
  int n;
  int ic;
  int mn;
  int b_mn;
  emxArray_real_T *A;
  double tau_data[2];
  int tau_size[1];
  int jpvt[2];
  double p1[2];
  emxArray_real_T *B;
  int m;
  int ia;
  double scale;
  double R_data[4];
  emxArray_real_T *r;
  unsigned int V_idx_0;
  int br;
  double absxk;
  double t;
  emxInit_real_T(&V, 2);
  n = x->size[1] - 1;
  ic = V->size[0] * V->size[1];
  V->size[0] = x->size[1];
  V->size[1] = 2;
  emxEnsureCapacity((emxArray__common *)V, ic, (int)sizeof(double));
  if (V->size[0] == 0) {
  } else {
    for (mn = 0; mn <= n; mn++) {
      V->data[mn + V->size[0]] = 1.0;
    }

    for (mn = 0; mn <= n; mn++) {
      V->data[mn] = x->data[mn];
    }
  }

  b_mn = V->size[0];
  if (b_mn <= 2) {
  } else {
    b_mn = 2;
  }

  emxInit_real_T(&A, 2);
  ic = A->size[0] * A->size[1];
  A->size[0] = V->size[0];
  A->size[1] = 2;
  emxEnsureCapacity((emxArray__common *)A, ic, (int)sizeof(double));
  mn = V->size[0] * V->size[1];
  for (ic = 0; ic < mn; ic++) {
    A->data[ic] = V->data[ic];
  }

  xgeqp3(A, tau_data, tau_size, jpvt);
  for (n = 0; n < 2; n++) {
    p1[n] = 0.0;
  }

  emxInit_real_T1(&B, 1);
  ic = B->size[0];
  B->size[0] = y->size[1];
  emxEnsureCapacity((emxArray__common *)B, ic, (int)sizeof(double));
  mn = y->size[1];
  for (ic = 0; ic < mn; ic++) {
    B->data[ic] = y->data[ic];
  }

  m = A->size[0];
  mn = A->size[0];
  if (mn <= 2) {
  } else {
    mn = 2;
  }

  for (ia = 0; ia + 1 <= mn; ia++) {
    if (tau_data[ia] != 0.0) {
      scale = B->data[ia];
      for (n = ia + 1; n + 1 <= m; n++) {
        scale += A->data[n + A->size[0] * ia] * B->data[n];
      }

      scale *= tau_data[ia];
      if (scale != 0.0) {
        B->data[ia] -= scale;
        for (n = ia + 1; n + 1 <= m; n++) {
          B->data[n] -= A->data[n + A->size[0] * ia] * scale;
        }
      }
    }
  }

  for (n = 0; n + 1 <= b_mn; n++) {
    p1[jpvt[n] - 1] = B->data[n];
  }

  emxFree_real_T(&B);
  for (ia = b_mn - 1; ia + 1 > 0; ia--) {
    p1[jpvt[ia] - 1] /= A->data[ia + A->size[0] * ia];
    n = 1;
    while (n <= ia) {
      p1[jpvt[0] - 1] -= p1[jpvt[ia] - 1] * A->data[A->size[0] * ia];
      n = 2;
    }
  }

  for (ia = 0; ia < 2; ia++) {
    if (ia + 1 <= b_mn) {
      ic = ia + 1;
    } else {
      ic = b_mn;
    }

    for (n = 0; n + 1 <= ic; n++) {
      R_data[n + (signed char)b_mn * ia] = A->data[n + A->size[0] * ia];
    }

    n = ia + 2;
    while (n <= b_mn) {
      R_data[1 + 2 * ia] = 0.0;
      n = 3;
    }
  }

  emxFree_real_T(&A);
  emxInit_real_T1(&r, 1);
  V_idx_0 = (unsigned int)V->size[0];
  ic = r->size[0];
  r->size[0] = (int)V_idx_0;
  emxEnsureCapacity((emxArray__common *)r, ic, (int)sizeof(double));
  m = V->size[0];
  mn = r->size[0];
  ic = r->size[0];
  r->size[0] = mn;
  emxEnsureCapacity((emxArray__common *)r, ic, (int)sizeof(double));
  for (ic = 0; ic < mn; ic++) {
    r->data[ic] = 0.0;
  }

  if (V->size[0] == 0) {
  } else {
    mn = 0;
    while ((m > 0) && (mn <= 0)) {
      for (ic = 1; ic <= m; ic++) {
        r->data[ic - 1] = 0.0;
      }

      mn = m;
    }

    br = 0;
    mn = 0;
    while ((m > 0) && (mn <= 0)) {
      n = -1;
      for (mn = br; mn + 1 <= br + 2; mn++) {
        if (p1[mn] != 0.0) {
          ia = n;
          for (ic = 0; ic + 1 <= m; ic++) {
            ia++;
            r->data[ic] += p1[mn] * V->data[ia];
          }
        }

        n += m;
      }

      br += 2;
      mn = m;
    }
  }

  emxFree_real_T(&V);
  ic = r->size[0];
  r->size[0] = y->size[1];
  emxEnsureCapacity((emxArray__common *)r, ic, (int)sizeof(double));
  mn = y->size[1];
  for (ic = 0; ic < mn; ic++) {
    r->data[ic] = y->data[ic] - r->data[ic];
  }

  S_R_size[0] = (signed char)b_mn;
  S_R_size[1] = 2;
  mn = (signed char)b_mn << 1;
  for (ic = 0; ic < mn; ic++) {
    S_R_data[ic] = R_data[ic];
  }

  if (0 >= y->size[1] - 2) {
    mn = 0;
  } else {
    mn = y->size[1] - 2;
  }

  *S_df = mn;
  if (r->size[0] == 0) {
    *S_normr = 0.0;
  } else {
    *S_normr = 0.0;
    if (r->size[0] == 1) {
      *S_normr = fabs(r->data[0]);
    } else {
      scale = 2.2250738585072014E-308;
      for (mn = 1; mn <= r->size[0]; mn++) {
        absxk = fabs(r->data[mn - 1]);
        if (absxk > scale) {
          t = scale / absxk;
          *S_normr = 1.0 + *S_normr * t * t;
          scale = absxk;
        } else {
          t = absxk / scale;
          *S_normr += t * t;
        }
      }

      *S_normr = scale * sqrt(*S_normr);
    }
  }

  emxFree_real_T(&r);
  for (ic = 0; ic < 2; ic++) {
    p[ic] = p1[ic];
  }
}

/*
 * File trailer for polyfit.c
 *
 * [EOF]
 */
