/*
 * File: xgeqp3.c
 *
 * MATLAB Coder version            : 3.1
 * C/C++ source code generated on  : 20-Dec-2017 17:38:07
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "PGain_from_LinearSlope.h"
#include "xgeqp3.h"
#include "xscal.h"
#include "xnrm2.h"

/* Function Declarations */
static double rt_hypotd_snf(double u0, double u1);

/* Function Definitions */

/*
 * Arguments    : double u0
 *                double u1
 * Return Type  : double
 */
static double rt_hypotd_snf(double u0, double u1)
{
  double y;
  double a;
  double b;
  a = fabs(u0);
  b = fabs(u1);
  if (a < b) {
    a /= b;
    y = b * sqrt(a * a + 1.0);
  } else if (a > b) {
    b /= a;
    y = a * sqrt(b * b + 1.0);
  } else if (rtIsNaN(b)) {
    y = b;
  } else {
    y = a * 1.4142135623730951;
  }

  return y;
}

/*
 * Arguments    : emxArray_real_T *A
 *                double tau_data[]
 *                int tau_size[1]
 *                int jpvt[2]
 * Return Type  : void
 */
void xgeqp3(emxArray_real_T *A, double tau_data[], int tau_size[1], int jpvt[2])
{
  int m;
  int mn;
  int i3;
  int k;
  double vn1[2];
  double vn2[2];
  int j;
  int i;
  double work[2];
  double xnorm;
  int i_i;
  int mmi;
  double beta1;
  int itemp;
  int pvt;
  double absxk;
  int ix;
  double t;
  int lastv;
  int lastc;
  boolean_T exitg2;
  int exitg1;
  m = A->size[0];
  if (A->size[0] <= 2) {
    mn = A->size[0];
  } else {
    mn = 2;
  }

  tau_size[0] = mn;
  for (i3 = 0; i3 < 2; i3++) {
    jpvt[i3] = 1 + i3;
  }

  if (A->size[0] == 0) {
  } else {
    k = 0;
    for (j = 0; j < 2; j++) {
      work[j] = 0.0;
      xnorm = 0.0;
      if (m < 1) {
      } else if (m == 1) {
        xnorm = fabs(A->data[k]);
      } else {
        beta1 = 2.2250738585072014E-308;
        pvt = k + m;
        for (itemp = k; itemp + 1 <= pvt; itemp++) {
          absxk = fabs(A->data[itemp]);
          if (absxk > beta1) {
            t = beta1 / absxk;
            xnorm = 1.0 + xnorm * t * t;
            beta1 = absxk;
          } else {
            t = absxk / beta1;
            xnorm += t * t;
          }
        }

        xnorm = beta1 * sqrt(xnorm);
      }

      vn1[j] = xnorm;
      vn2[j] = vn1[j];
      k += m;
    }

    for (i = 0; i + 1 <= mn; i++) {
      i_i = i + i * m;
      mmi = (m - i) - 1;
      itemp = 0;
      if ((2 - i > 1) && (vn1[1] > vn1[i])) {
        itemp = 1;
      }

      pvt = i + itemp;
      if (pvt + 1 != i + 1) {
        ix = m * pvt;
        j = m * i;
        for (k = 1; k <= m; k++) {
          xnorm = A->data[ix];
          A->data[ix] = A->data[j];
          A->data[j] = xnorm;
          ix++;
          j++;
        }

        itemp = jpvt[pvt];
        jpvt[pvt] = jpvt[i];
        jpvt[i] = itemp;
        vn1[pvt] = vn1[i];
        vn2[pvt] = vn2[i];
      }

      if (i + 1 < m) {
        t = A->data[i_i];
        absxk = 0.0;
        if (1 + mmi <= 0) {
        } else {
          xnorm = xnrm2(mmi, A, i_i + 2);
          if (xnorm != 0.0) {
            beta1 = rt_hypotd_snf(A->data[i_i], xnorm);
            if (A->data[i_i] >= 0.0) {
              beta1 = -beta1;
            }

            if (fabs(beta1) < 1.0020841800044864E-292) {
              itemp = 0;
              do {
                itemp++;
                xscal(mmi, 9.9792015476736E+291, A, i_i + 2);
                beta1 *= 9.9792015476736E+291;
                t *= 9.9792015476736E+291;
              } while (!(fabs(beta1) >= 1.0020841800044864E-292));

              xnorm = xnrm2(mmi, A, i_i + 2);
              beta1 = rt_hypotd_snf(t, xnorm);
              if (t >= 0.0) {
                beta1 = -beta1;
              }

              absxk = (beta1 - t) / beta1;
              xscal(mmi, 1.0 / (t - beta1), A, i_i + 2);
              for (k = 1; k <= itemp; k++) {
                beta1 *= 1.0020841800044864E-292;
              }

              t = beta1;
            } else {
              absxk = (beta1 - A->data[i_i]) / beta1;
              xscal(mmi, 1.0 / (A->data[i_i] - beta1), A, i_i + 2);
              t = beta1;
            }
          }
        }

        tau_data[i] = absxk;
        A->data[i_i] = t;
      } else {
        tau_data[i] = 0.0;
      }

      if (i + 1 < 2) {
        t = A->data[i_i];
        A->data[i_i] = 1.0;
        if (tau_data[0] != 0.0) {
          lastv = 1 + mmi;
          itemp = i_i + mmi;
          while ((lastv > 0) && (A->data[itemp] == 0.0)) {
            lastv--;
            itemp--;
          }

          lastc = 1;
          exitg2 = false;
          while ((!exitg2) && (lastc > 0)) {
            itemp = m;
            do {
              exitg1 = 0;
              if (itemp + 1 <= m + lastv) {
                if (A->data[itemp] != 0.0) {
                  exitg1 = 1;
                } else {
                  itemp++;
                }
              } else {
                lastc = 0;
                exitg1 = 2;
              }
            } while (exitg1 == 0);

            if (exitg1 == 1) {
              exitg2 = true;
            }
          }
        } else {
          lastv = 0;
          lastc = 0;
        }

        if (lastv > 0) {
          if (lastc == 0) {
          } else {
            work[0] = 0.0;
            j = 0;
            pvt = 1 + m;
            while ((m > 0) && (pvt <= 1 + m)) {
              ix = i_i;
              xnorm = 0.0;
              i3 = (pvt + lastv) - 1;
              for (itemp = pvt; itemp <= i3; itemp++) {
                xnorm += A->data[itemp - 1] * A->data[ix];
                ix++;
              }

              work[j] += xnorm;
              j++;
              pvt += m;
            }
          }

          if (-tau_data[0] == 0.0) {
          } else {
            pvt = m;
            k = 0;
            j = 1;
            while (j <= lastc) {
              if (work[k] != 0.0) {
                xnorm = work[k] * -tau_data[0];
                ix = i_i;
                i3 = lastv + pvt;
                for (itemp = pvt; itemp + 1 <= i3; itemp++) {
                  A->data[itemp] += A->data[ix] * xnorm;
                  ix++;
                }
              }

              k++;
              pvt += m;
              j = 2;
            }
          }
        }

        A->data[i_i] = t;
      }

      j = i + 2;
      while (j < 3) {
        itemp = (i + m) + 1;
        if (vn1[1] != 0.0) {
          xnorm = fabs(A->data[i + A->size[0]]) / vn1[1];
          xnorm = 1.0 - xnorm * xnorm;
          if (xnorm < 0.0) {
            xnorm = 0.0;
          }

          beta1 = vn1[1] / vn2[1];
          beta1 = xnorm * (beta1 * beta1);
          if (beta1 <= 1.4901161193847656E-8) {
            if (i + 1 < m) {
              xnorm = 0.0;
              if (mmi < 1) {
              } else if (mmi == 1) {
                xnorm = fabs(A->data[itemp]);
              } else {
                beta1 = 2.2250738585072014E-308;
                pvt = itemp + mmi;
                while (itemp + 1 <= pvt) {
                  absxk = fabs(A->data[itemp]);
                  if (absxk > beta1) {
                    t = beta1 / absxk;
                    xnorm = 1.0 + xnorm * t * t;
                    beta1 = absxk;
                  } else {
                    t = absxk / beta1;
                    xnorm += t * t;
                  }

                  itemp++;
                }

                xnorm = beta1 * sqrt(xnorm);
              }

              vn1[1] = xnorm;
              vn2[1] = xnorm;
            } else {
              vn1[1] = 0.0;
              vn2[1] = 0.0;
            }
          } else {
            vn1[1] *= sqrt(xnorm);
          }
        }

        j = 3;
      }
    }
  }
}

/*
 * File trailer for xgeqp3.c
 *
 * [EOF]
 */
