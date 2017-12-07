/*
 * File: LinearSlope4Manual.c
 *
 * MATLAB Coder version            : 3.1
 * C/C++ source code generated on  : 20-Dec-2017 18:14:33
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "LinearSlope4Manual.h"
#include "LinearSlope4Manual_emxutil.h"
#include "sum.h"
#include "power.h"
#include "polyfit.h"
#include "log10.h"

/* Function Definitions */

/*
 * Input parameters:
 *
 *    f - Frequencies vector. Assumed to have fixed length of 1000 (needs fixed length for codegen)
 *    r - Gains vector, in [dB]. As above.
 *    n - Actual length of the f and r vectors.
 *    MinFreqRange - Start frequency range.
 *    MaxFreqRange - End frequency range.
 *    Output parameters:
 * Output---------
 *    SlopeResult - the calculated slope for this frequency range.
 *    PlantGain - Calculated gain of plant.
 *    ErrCode --- if negative value, Daniel will display the ErrCode;
 *            --- [-3]--- input value out of range;
 *    All output parameters are 0 if the function fails to find a solution.
 *    plantGain calculation is based on each freq point(w) and G/(jw)^2,then make average of the whole interpolated freq points
 * Arguments    : const double f[1000]
 *                const double r[1000]
 *                double n
 *                double MinFreqRange
 *                double MaxFreqRange
 *                double *SlopeResult
 *                double *PlantGain
 *                double *ErrCode
 * Return Type  : void
 */
void LinearSlope4Manual(const double f[1000], const double r[1000], double n,
  double MinFreqRange, double MaxFreqRange, double *SlopeResult, double
  *PlantGain, double *ErrCode)
{
  int low_ip1;
  double anew;
  double apnd;
  double ndbl;
  emxArray_real_T *fnew;
  double cdiff;
  double absa;
  int mid_i;
  double absb;
  int nm1d2;
  int k;
  double varargin_2_data[1000];
  double outsize[2];
  double x_data[1000];
  emxArray_real_T *rdb;
  emxArray_boolean_T *x;
  int exitg2;
  emxArray_int32_T *ii;
  int high_i;
  boolean_T exitg1;
  boolean_T guard1 = false;
  emxArray_real_T *idxG;
  emxArray_real_T *xx;
  emxArray_real_T *b_rdb;
  double expl_temp_data[4];
  int expl_temp_size[2];
  double gainLineRange[1000];
  emxArray_real_T *c_rdb;

  /* init output */
  *SlopeResult = 0.0;
  *PlantGain = 0.0;

  /*  var = f; */
  /*  if isempty(var) || isnan(var)  */
  /*      ErrCode  = -3; */
  /*  end */
  *ErrCode = 0.0;
  if (rtIsNaN(MaxFreqRange)) {
    *ErrCode = -3.0;
  }

  if (*ErrCode < 0.0) {
  } else {
    if ((MinFreqRange <= 0.0) || (MaxFreqRange <= 0.0)) {
      *ErrCode = -4.0;
    }

    if (*ErrCode < 0.0) {
    } else {
      if (rtIsNaN(f[0]) || rtIsNaN(f[(int)n - 1])) {
        low_ip1 = 1;
        anew = rtNaN;
        apnd = f[(int)n - 1];
      } else if (f[(int)n - 1] < f[0]) {
        low_ip1 = 0;
        anew = f[0];
        apnd = f[(int)n - 1];
      } else if (rtIsInf(f[0]) || rtIsInf(f[(int)n - 1])) {
        low_ip1 = 1;
        anew = rtNaN;
        apnd = f[(int)n - 1];
      } else {
        anew = f[0];
        ndbl = floor((f[(int)n - 1] - f[0]) / 0.1 + 0.5);
        apnd = f[0] + ndbl * 0.1;
        cdiff = apnd - f[(int)n - 1];
        absa = fabs(f[0]);
        absb = fabs(f[(int)n - 1]);
        if ((absa >= absb) || rtIsNaN(absb)) {
          absb = absa;
        }

        if (fabs(cdiff) < 4.4408920985006262E-16 * absb) {
          ndbl++;
          apnd = f[(int)n - 1];
        } else if (cdiff > 0.0) {
          apnd = f[0] + (ndbl - 1.0) * 0.1;
        } else {
          ndbl++;
        }

        if (ndbl >= 0.0) {
          low_ip1 = (int)ndbl;
        } else {
          low_ip1 = 0;
        }
      }

      emxInit_real_T(&fnew, 2);
      mid_i = fnew->size[0] * fnew->size[1];
      fnew->size[0] = 1;
      fnew->size[1] = low_ip1;
      emxEnsureCapacity((emxArray__common *)fnew, mid_i, (int)sizeof(double));
      if (low_ip1 > 0) {
        fnew->data[0] = anew;
        if (low_ip1 > 1) {
          fnew->data[low_ip1 - 1] = apnd;
          nm1d2 = (low_ip1 - 1) / 2;
          for (k = 1; k < nm1d2; k++) {
            cdiff = (double)k * 0.1;
            fnew->data[k] = anew + cdiff;
            fnew->data[(low_ip1 - k) - 1] = apnd - cdiff;
          }

          if (nm1d2 << 1 == low_ip1 - 1) {
            fnew->data[nm1d2] = (anew + apnd) / 2.0;
          } else {
            cdiff = (double)nm1d2 * 0.1;
            fnew->data[nm1d2] = anew + cdiff;
            fnew->data[nm1d2 + 1] = apnd - cdiff;
          }
        }
      }

      /* % 1*XX  */
      /* % */
      /*  nLength = length(fnew); */
      /*  % fnew = fnew(:);%% do we need to XX*1 ? temp, will uncomment later; */
      /*  % % // Den = reshape(Den,1,numel(Den));%%1*XX matrix transopse eg */
      nm1d2 = (int)n;
      for (mid_i = 0; mid_i < nm1d2; mid_i++) {
        varargin_2_data[mid_i] = r[mid_i];
      }

      nm1d2 = (int)n;
      for (mid_i = 0; mid_i < nm1d2; mid_i++) {
        x_data[mid_i] = f[mid_i];
      }

      for (mid_i = 0; mid_i < 2; mid_i++) {
        outsize[mid_i] = fnew->size[mid_i];
      }

      emxInit_real_T(&rdb, 2);
      mid_i = rdb->size[0] * rdb->size[1];
      rdb->size[0] = 1;
      rdb->size[1] = (int)outsize[1];
      emxEnsureCapacity((emxArray__common *)rdb, mid_i, (int)sizeof(double));
      nm1d2 = (int)outsize[1];
      for (mid_i = 0; mid_i < nm1d2; mid_i++) {
        rdb->data[mid_i] = rtNaN;
      }

      if (fnew->size[1] == 0) {
      } else {
        k = 1;
        do {
          exitg2 = 0;
          if (k <= (int)n) {
            if (rtIsNaN(f[k - 1])) {
              exitg2 = 1;
            } else {
              k++;
            }
          } else {
            if (f[1] < f[0]) {
              mid_i = (int)n >> 1;
              for (nm1d2 = 1; nm1d2 <= mid_i; nm1d2++) {
                cdiff = x_data[nm1d2 - 1];
                x_data[nm1d2 - 1] = x_data[(int)n - nm1d2];
                x_data[(int)n - nm1d2] = cdiff;
              }

              if ((int)n > 1) {
                nm1d2 = (int)n >> 1;
                for (k = 1; k <= nm1d2; k++) {
                  ndbl = varargin_2_data[k - 1];
                  varargin_2_data[k - 1] = varargin_2_data[(int)n - k];
                  varargin_2_data[(int)n - k] = ndbl;
                }
              }
            }

            for (k = 0; k + 1 <= fnew->size[1]; k++) {
              cdiff = rdb->data[k];
              if (rtIsNaN(fnew->data[k])) {
                cdiff = rtNaN;
              } else if ((fnew->data[k] > x_data[(int)n - 1]) || (fnew->data[k] <
                          x_data[0])) {
              } else {
                nm1d2 = 1;
                low_ip1 = 2;
                high_i = (int)n;
                while (high_i > low_ip1) {
                  mid_i = (nm1d2 >> 1) + (high_i >> 1);
                  if (((nm1d2 & 1) == 1) && ((high_i & 1) == 1)) {
                    mid_i++;
                  }

                  if (fnew->data[k] >= x_data[mid_i - 1]) {
                    nm1d2 = mid_i;
                    low_ip1 = mid_i + 1;
                  } else {
                    high_i = mid_i;
                  }
                }

                cdiff = (fnew->data[k] - x_data[nm1d2 - 1]) / (x_data[nm1d2] -
                  x_data[nm1d2 - 1]);
                if (cdiff == 0.0) {
                  cdiff = varargin_2_data[nm1d2 - 1];
                } else if (cdiff == 1.0) {
                  cdiff = varargin_2_data[nm1d2];
                } else if (varargin_2_data[nm1d2 - 1] == varargin_2_data[nm1d2])
                {
                  cdiff = varargin_2_data[nm1d2 - 1];
                } else {
                  cdiff = (1.0 - cdiff) * varargin_2_data[nm1d2 - 1] + cdiff *
                    varargin_2_data[nm1d2];
                }
              }

              rdb->data[k] = cdiff;
            }

            exitg2 = 1;
          }
        } while (exitg2 == 0);
      }

      emxInit_boolean_T(&x, 2);
      mid_i = x->size[0] * x->size[1];
      x->size[0] = 1;
      x->size[1] = fnew->size[1];
      emxEnsureCapacity((emxArray__common *)x, mid_i, (int)sizeof(boolean_T));
      nm1d2 = fnew->size[0] * fnew->size[1];
      for (mid_i = 0; mid_i < nm1d2; mid_i++) {
        x->data[mid_i] = ((fnew->data[mid_i] > MinFreqRange) && (fnew->
          data[mid_i] < MaxFreqRange));
      }

      emxInit_int32_T(&ii, 2);
      low_ip1 = x->size[1];
      high_i = 0;
      mid_i = ii->size[0] * ii->size[1];
      ii->size[0] = 1;
      ii->size[1] = x->size[1];
      emxEnsureCapacity((emxArray__common *)ii, mid_i, (int)sizeof(int));
      nm1d2 = 1;
      exitg1 = false;
      while ((!exitg1) && (nm1d2 <= low_ip1)) {
        guard1 = false;
        if (x->data[nm1d2 - 1]) {
          high_i++;
          ii->data[high_i - 1] = nm1d2;
          if (high_i >= low_ip1) {
            exitg1 = true;
          } else {
            guard1 = true;
          }
        } else {
          guard1 = true;
        }

        if (guard1) {
          nm1d2++;
        }
      }

      if (x->size[1] == 1) {
        if (high_i == 0) {
          mid_i = ii->size[0] * ii->size[1];
          ii->size[0] = 1;
          ii->size[1] = 0;
          emxEnsureCapacity((emxArray__common *)ii, mid_i, (int)sizeof(int));
        }
      } else {
        mid_i = ii->size[0] * ii->size[1];
        if (1 > high_i) {
          ii->size[1] = 0;
        } else {
          ii->size[1] = high_i;
        }

        emxEnsureCapacity((emxArray__common *)ii, mid_i, (int)sizeof(int));
      }

      emxFree_boolean_T(&x);
      emxInit_real_T(&idxG, 2);
      mid_i = idxG->size[0] * idxG->size[1];
      idxG->size[0] = 1;
      idxG->size[1] = ii->size[1];
      emxEnsureCapacity((emxArray__common *)idxG, mid_i, (int)sizeof(double));
      nm1d2 = ii->size[0] * ii->size[1];
      for (mid_i = 0; mid_i < nm1d2; mid_i++) {
        idxG->data[mid_i] = ii->data[mid_i];
      }

      emxFree_int32_T(&ii);
      if (idxG->size[1] == 0) {
        *ErrCode = -1.0;
      } else if (idxG->size[1] > fnew->size[1]) {
        *ErrCode = -2.0;
      } else {
        if ((int)idxG->data[0] > (int)idxG->data[idxG->size[1] - 1]) {
          mid_i = 0;
          low_ip1 = 0;
        } else {
          mid_i = (int)idxG->data[0] - 1;
          low_ip1 = (int)idxG->data[idxG->size[1] - 1];
        }

        emxInit_real_T(&xx, 2);
        nm1d2 = xx->size[0] * xx->size[1];
        xx->size[0] = 1;
        xx->size[1] = low_ip1 - mid_i;
        emxEnsureCapacity((emxArray__common *)xx, nm1d2, (int)sizeof(double));
        nm1d2 = low_ip1 - mid_i;
        for (low_ip1 = 0; low_ip1 < nm1d2; low_ip1++) {
          xx->data[xx->size[0] * low_ip1] = fnew->data[mid_i + low_ip1];
        }

        b_log10(xx);
        if ((int)idxG->data[0] > (int)idxG->data[idxG->size[1] - 1]) {
          mid_i = 0;
          low_ip1 = 0;
        } else {
          mid_i = (int)idxG->data[0] - 1;
          low_ip1 = (int)idxG->data[idxG->size[1] - 1];
        }

        emxInit_real_T(&b_rdb, 2);
        nm1d2 = b_rdb->size[0] * b_rdb->size[1];
        b_rdb->size[0] = 1;
        b_rdb->size[1] = low_ip1 - mid_i;
        emxEnsureCapacity((emxArray__common *)b_rdb, nm1d2, (int)sizeof(double));
        nm1d2 = low_ip1 - mid_i;
        for (low_ip1 = 0; low_ip1 < nm1d2; low_ip1++) {
          b_rdb->data[b_rdb->size[0] * low_ip1] = rdb->data[mid_i + low_ip1];
        }

        polyfit(xx, b_rdb, outsize, expl_temp_data, expl_temp_size, &cdiff,
                &ndbl);
        *SlopeResult = outsize[0];

        /* PlantGain22 = (2*pi)^2*10^(pCoefficient(2)/20) */
        /* %%for  plant gain calculation; */
        emxFree_real_T(&b_rdb);
        emxFree_real_T(&xx);
        memset(&gainLineRange[0], 0, 1000U * sizeof(double));
        emxInit_real_T(&c_rdb, 2);
        mid_i = c_rdb->size[0] * c_rdb->size[1];
        c_rdb->size[0] = 1;
        c_rdb->size[1] = rdb->size[1];
        emxEnsureCapacity((emxArray__common *)c_rdb, mid_i, (int)sizeof(double));
        nm1d2 = rdb->size[0] * rdb->size[1];
        for (mid_i = 0; mid_i < nm1d2; mid_i++) {
          c_rdb->data[mid_i] = rdb->data[mid_i] / 20.0;
        }

        power(c_rdb, rdb);

        /* %non-dB */
        nm1d2 = 0;
        emxFree_real_T(&c_rdb);
        while (nm1d2 <= idxG->size[1] - 1) {
          gainLineRange[nm1d2] = fabs(rdb->data[(int)idxG->data[nm1d2] - 1]) *
            fnew->data[(int)idxG->data[nm1d2] - 1] * 2.0 * 3.1415926535897931 *
            fnew->data[(int)idxG->data[nm1d2] - 1] * 2.0 * 3.1415926535897931;
          nm1d2++;
        }

        /* %can not use mean of matlab for the zero value  */
        *PlantGain = sum(gainLineRange) / (double)idxG->size[1];
      }

      emxFree_real_T(&idxG);
      emxFree_real_T(&rdb);
      emxFree_real_T(&fnew);
    }
  }
}

/*
 * File trailer for LinearSlope4Manual.c
 *
 * [EOF]
 */
