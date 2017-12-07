/*
 * File: interp1.c
 *
 * MATLAB Coder version            : 3.1
 * C/C++ source code generated on  : 20-Dec-2017 17:38:07
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "PGain_from_LinearSlope.h"
#include "interp1.h"
#include "PGain_from_LinearSlope_emxutil.h"

/* Function Definitions */

/*
 * Arguments    : const double varargin_1_data[]
 *                const int varargin_1_size[1]
 *                const double varargin_2_data[]
 *                const int varargin_2_size[1]
 *                const emxArray_real_T *varargin_3
 *                emxArray_real_T *Vq
 * Return Type  : void
 */
void interp1(const double varargin_1_data[], const int varargin_1_size[1], const
             double varargin_2_data[], const int varargin_2_size[1], const
             emxArray_real_T *varargin_3, emxArray_real_T *Vq)
{
  int n;
  int nd2;
  int x_size_idx_0;
  double y_data[1000];
  int nx;
  double x_data[1000];
  unsigned int outsize_idx_0;
  int k;
  int exitg1;
  double r;
  int mid_i;
  n = varargin_2_size[0];
  for (nd2 = 0; nd2 < n; nd2++) {
    y_data[nd2] = varargin_2_data[nd2];
  }

  x_size_idx_0 = varargin_1_size[0];
  n = varargin_1_size[0];
  for (nd2 = 0; nd2 < n; nd2++) {
    x_data[nd2] = varargin_1_data[nd2];
  }

  nx = varargin_1_size[0];
  outsize_idx_0 = (unsigned int)varargin_3->size[0];
  nd2 = Vq->size[0];
  Vq->size[0] = (int)outsize_idx_0;
  emxEnsureCapacity((emxArray__common *)Vq, nd2, (int)sizeof(double));
  n = (int)outsize_idx_0;
  for (nd2 = 0; nd2 < n; nd2++) {
    Vq->data[nd2] = rtNaN;
  }

  if (varargin_3->size[0] == 0) {
  } else {
    k = 1;
    do {
      exitg1 = 0;
      if (k <= nx) {
        if (rtIsNaN(varargin_1_data[k - 1])) {
          exitg1 = 1;
        } else {
          k++;
        }
      } else {
        if (varargin_1_data[1] < varargin_1_data[0]) {
          nd2 = nx >> 1;
          for (n = 1; n <= nd2; n++) {
            r = x_data[n - 1];
            x_data[n - 1] = x_data[nx - n];
            x_data[nx - n] = r;
          }

          if ((!(varargin_2_size[0] == 0)) && (varargin_2_size[0] > 1)) {
            n = varargin_2_size[0];
            nd2 = varargin_2_size[0] >> 1;
            for (k = 1; k <= nd2; k++) {
              r = y_data[k - 1];
              y_data[k - 1] = y_data[n - k];
              y_data[n - k] = r;
            }
          }
        }

        for (k = 0; k + 1 <= varargin_3->size[0]; k++) {
          r = Vq->data[k];
          if (rtIsNaN(varargin_3->data[k])) {
            r = rtNaN;
          } else if ((varargin_3->data[k] > x_data[x_size_idx_0 - 1]) ||
                     (varargin_3->data[k] < x_data[0])) {
          } else {
            n = 1;
            nd2 = 2;
            nx = x_size_idx_0;
            while (nx > nd2) {
              mid_i = (n >> 1) + (nx >> 1);
              if (((n & 1) == 1) && ((nx & 1) == 1)) {
                mid_i++;
              }

              if (varargin_3->data[k] >= x_data[mid_i - 1]) {
                n = mid_i;
                nd2 = mid_i + 1;
              } else {
                nx = mid_i;
              }
            }

            r = (varargin_3->data[k] - x_data[n - 1]) / (x_data[n] - x_data[n -
              1]);
            if (r == 0.0) {
              r = y_data[n - 1];
            } else if (r == 1.0) {
              r = y_data[n];
            } else if (y_data[n - 1] == y_data[n]) {
              r = y_data[n - 1];
            } else {
              r = (1.0 - r) * y_data[n - 1] + r * y_data[n];
            }
          }

          Vq->data[k] = r;
        }

        exitg1 = 1;
      }
    } while (exitg1 == 0);
  }
}

/*
 * File trailer for interp1.c
 *
 * [EOF]
 */
