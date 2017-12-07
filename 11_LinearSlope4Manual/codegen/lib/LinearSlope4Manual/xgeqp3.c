/*
 * File: xgeqp3.c
 *
 * MATLAB Coder version            : 3.1
 * C/C++ source code generated on  : 20-Dec-2017 18:14:33
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "LinearSlope4Manual.h"
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
  int i0;
  int k;
  double vn1[2];
  double vn2[2];
  int j;
  int i;
  double work[2];
  double smax;
  int i_i;
  int mmi;
  double temp2;
  int pvt;
  int itemp;
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
  for (i0 = 0; i0 < 2; i0++) {
    jpvt[i0] = 1 + i0;
  }

  if (A->size[0] == 0) {
  } else {
    k = 0;
    for (j = 0; j < 2; j++) {
      work[j] = 0.0;
      smax = 0.0;
      if (m < 1) {
      } else if (m == 1) {
        smax = fabs(A->data[k]);
      } else {
        temp2 = 2.2250738585072014E-308;
        pvt = k + m;
        for (itemp = k; itemp + 1 <= pvt; itemp++) {
          absxk = fabs(A->data[itemp]);
          if (absxk > temp2) {
            t = temp2 / absxk;
            smax = 1.0 + smax * t * t;
            temp2 = absxk;
          } else {
            t = absxk / temp2;
            smax += t * t;
          }
        }

        smax = temp2 * sqrt(smax);
      }

      vn1[j] = smax;
      vn2[j] = vn1[j];
      k += m;
    }

    for (i = 0; i + 1 <= mn; i++) {
      i_i = i + i * m;
      mmi = (m - i) - 1;
      if (2 - i < 1) {
        itemp = -1;
      } else {
        itemp = 0;
        if (2 - i > 1) {
          smax = vn1[i];
          k = 2;
          while (k <= 2 - i) {
            if (vn1[1] > smax) {
              itemp = 1;
              smax = vn1[1];
            }

            k = 3;
          }
        }
      }

      pvt = i + itemp;
      if (pvt + 1 != i + 1) {
        ix = m * pvt;
        j = m * i;
        for (k = 1; k <= m; k++) {
          smax = A->data[ix];
          A->data[ix] = A->data[j];
          A->data[j] = smax;
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
        absxk = A->data[i_i];
        temp2 = 0.0;
        if (1 + mmi <= 0) {
        } else {
          smax = xnrm2(mmi, A, i_i + 2);
          if (smax != 0.0) {
            smax = rt_hypotd_snf(A->data[i_i], smax);
            if (A->data[i_i] >= 0.0) {
              smax = -smax;
            }

            if (fabs(smax) < 1.0020841800044864E-292) {
              itemp = 0;
              do {
                itemp++;
                xscal(mmi, 9.9792015476736E+291, A, i_i + 2);
                smax *= 9.9792015476736E+291;
                absxk *= 9.9792015476736E+291;
              } while (!(fabs(smax) >= 1.0020841800044864E-292));

              smax = xnrm2(mmi, A, i_i + 2);
              smax = rt_hypotd_snf(absxk, smax);
              if (absxk >= 0.0) {
                smax = -smax;
              }

              temp2 = (smax - absxk) / smax;
              xscal(mmi, 1.0 / (absxk - smax), A, i_i + 2);
              for (k = 1; k <= itemp; k++) {
                smax *= 1.0020841800044864E-292;
              }

              absxk = smax;
            } else {
              temp2 = (smax - A->data[i_i]) / smax;
              xscal(mmi, 1.0 / (A->data[i_i] - smax), A, i_i + 2);
              absxk = smax;
            }
          }
        }

        tau_data[i] = temp2;
        A->data[i_i] = absxk;
      } else {
        tau_data[i] = 0.0;
      }

      if (i + 1 < 2) {
        absxk = A->data[i_i];
        A->data[i_i] = 1.0;
        if (tau_data[0] != 0.0) {
          lastv = mmi + 1;
          itemp = i_i + mmi;
          while ((lastv > 0) && (A->data[itemp] == 0.0)) {
            lastv--;
            itemp--;
          }

          lastc = 1;
          exitg2 = false;
          while ((!exitg2) && (lastc > 0)) {
            itemp = m + 1;
            do {
              exitg1 = 0;
              if (itemp <= m + lastv) {
                if (A->data[itemp - 1] != 0.0) {
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
            while ((m > 0) && (pvt <= m + 1)) {
              ix = i_i;
              smax = 0.0;
              i0 = (pvt + lastv) - 1;
              for (itemp = pvt; itemp <= i0; itemp++) {
                smax += A->data[itemp - 1] * A->data[ix];
                ix++;
              }

              work[j] += smax;
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
                smax = work[k] * -tau_data[0];
                ix = i_i;
                i0 = lastv + pvt;
                for (itemp = pvt; itemp + 1 <= i0; itemp++) {
                  A->data[itemp] += A->data[ix] * smax;
                  ix++;
                }
              }

              k++;
              pvt += m;
              j = 2;
            }
          }
        }

        A->data[i_i] = absxk;
      }

      j = i + 2;
      while (j < 3) {
        itemp = (i + m) + 1;
        if (vn1[1] != 0.0) {
          smax = fabs(A->data[i + A->size[0]]) / vn1[1];
          smax = 1.0 - smax * smax;
          if (smax < 0.0) {
            smax = 0.0;
          }

          temp2 = vn1[1] / vn2[1];
          temp2 = smax * (temp2 * temp2);
          if (temp2 <= 1.4901161193847656E-8) {
            if (i + 1 < m) {
              smax = 0.0;
              if (mmi < 1) {
              } else if (mmi == 1) {
                smax = fabs(A->data[itemp]);
              } else {
                temp2 = 2.2250738585072014E-308;
                pvt = itemp + mmi;
                while (itemp + 1 <= pvt) {
                  absxk = fabs(A->data[itemp]);
                  if (absxk > temp2) {
                    t = temp2 / absxk;
                    smax = 1.0 + smax * t * t;
                    temp2 = absxk;
                  } else {
                    t = absxk / temp2;
                    smax += t * t;
                  }

                  itemp++;
                }

                smax = temp2 * sqrt(smax);
              }

              vn1[1] = smax;
              vn2[1] = smax;
            } else {
              vn1[1] = 0.0;
              vn2[1] = 0.0;
            }
          } else {
            vn1[1] *= sqrt(smax);
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
