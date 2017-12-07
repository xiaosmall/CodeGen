/*
 * File: LinearSlope4Manual.c
 *
 * MATLAB Coder version            : 2.6
 * C/C++ source code generated on  : 28-Feb-2016 19:56:09
 */

/* Include files */
#include "rt_nonfinite.h"
#include "LinearSlope4Manual.h"
#include "LinearSlope4Manual_emxutil.h"
#include "sum.h"
#include "polyfit.h"
#include "interp1.h"

/* Function Declarations */
static double rt_powd_snf(double u0, double u1);

/* Function Definitions */

/*
 * Arguments    : double u0
 *                double u1
 * Return Type  : double
 */
static double rt_powd_snf(double u0, double u1)
{
  double y;
  double d0;
  double d1;
  if (rtIsNaN(u0) || rtIsNaN(u1)) {
    y = rtNaN;
  } else {
    d0 = fabs(u0);
    d1 = fabs(u1);
    if (rtIsInf(u1)) {
      if (d0 == 1.0) {
        y = rtNaN;
      } else if (d0 > 1.0) {
        if (u1 > 0.0) {
          y = rtInf;
        } else {
          y = 0.0;
        }
      } else if (u1 > 0.0) {
        y = 0.0;
      } else {
        y = rtInf;
      }
    } else if (d1 == 0.0) {
      y = 1.0;
    } else if (d1 == 1.0) {
      if (u1 > 0.0) {
        y = u0;
      } else {
        y = 1.0 / u0;
      }
    } else if (u1 == 2.0) {
      y = u0 * u0;
    } else if ((u1 == 0.5) && (u0 >= 0.0)) {
      y = sqrt(u0);
    } else if ((u0 < 0.0) && (u1 > floor(u1))) {
      y = rtNaN;
    } else {
      y = pow(u0, u1);
    }
  }

  return y;
}

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
  int idx;
  double anew;
  double apnd;
  double ndbl;
  double cdiff;
  double absa;
  double absb;
  emxArray_real_T *fnew;
  int i0;
  int nm1d2;
  int k;
  emxArray_real_T *rdb;
  double f_data[1000];
  int f_size[1];
  double r_data[1000];
  int r_size[1];
  emxArray_real_T *r0;
  emxArray_boolean_T *x;
  emxArray_int32_T *ii;
  boolean_T exitg1;
  boolean_T guard1 = false;
  emxArray_int32_T *b_ii;
  emxArray_real_T *idxG;
  emxArray_real_T *xx;
  emxArray_real_T *b_xx;
  emxArray_real_T *b_rdb;
  int expl_temp_size[2];
  double expl_temp_data[4];
  double pCoefficient[2];
  static double rnondb[1000000];
  static double gainLineRange[1000000];

  /* init output */
  *SlopeResult = 0.0;
  *PlantGain = 0.0;
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
        idx = 0;
        anew = rtNaN;
        apnd = f[(int)n - 1];
      } else if (f[(int)n - 1] < f[0]) {
        idx = -1;
        anew = f[0];
        apnd = f[(int)n - 1];
      } else if (rtIsInf(f[0]) || rtIsInf(f[(int)n - 1])) {
        idx = 0;
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
          idx = (int)ndbl - 1;
        } else {
          idx = -1;
        }
      }

      emxInit_real_T(&fnew, 2);
      i0 = fnew->size[0] * fnew->size[1];
      fnew->size[0] = 1;
      fnew->size[1] = idx + 1;
      emxEnsureCapacity((emxArray__common *)fnew, i0, (int)sizeof(double));
      if (idx + 1 > 0) {
        fnew->data[0] = anew;
        if (idx + 1 > 1) {
          fnew->data[idx] = apnd;
          nm1d2 = idx / 2;
          for (k = 1; k < nm1d2; k++) {
            cdiff = (double)k * 0.1;
            fnew->data[k] = anew + cdiff;
            fnew->data[idx - k] = apnd - cdiff;
          }

          if (nm1d2 << 1 == idx) {
            fnew->data[nm1d2] = (anew + apnd) / 2.0;
          } else {
            cdiff = (double)nm1d2 * 0.1;
            fnew->data[nm1d2] = anew + cdiff;
            fnew->data[nm1d2 + 1] = apnd - cdiff;
          }
        }
      }

      emxInit_real_T(&rdb, 2);

      /* % 1*XX  */
      /* % */
      /*  % fnew = fnew(:);%% do we need to XX*1 ? temp, will uncomment later; */
      /*  % % // Den = reshape(Den,1,numel(Den));%%1*XX matrix transopse eg */
      i0 = rdb->size[0] * rdb->size[1];
      rdb->size[0] = 1;
      emxEnsureCapacity((emxArray__common *)rdb, i0, (int)sizeof(double));
      nm1d2 = fnew->size[1];
      i0 = rdb->size[0] * rdb->size[1];
      rdb->size[1] = nm1d2;
      emxEnsureCapacity((emxArray__common *)rdb, i0, (int)sizeof(double));
      nm1d2 = fnew->size[1];
      for (i0 = 0; i0 < nm1d2; i0++) {
        rdb->data[i0] = 0.0;
      }

      /* init for codegen */
      f_size[0] = (int)n;
      nm1d2 = (int)n;
      for (i0 = 0; i0 < nm1d2; i0++) {
        f_data[i0] = f[i0];
      }

      r_size[0] = (int)n;
      nm1d2 = (int)n;
      for (i0 = 0; i0 < nm1d2; i0++) {
        r_data[i0] = r[i0];
      }

      emxInit_real_T(&r0, 2);
      interp1(f_data, f_size, r_data, r_size, fnew, r0);
      nm1d2 = r0->size[1];
      for (i0 = 0; i0 < nm1d2; i0++) {
        rdb->data[i0] = r0->data[r0->size[0] * i0];
      }

      emxInit_boolean_T(&x, 2);
      i0 = x->size[0] * x->size[1];
      x->size[0] = 1;
      x->size[1] = fnew->size[1];
      emxEnsureCapacity((emxArray__common *)x, i0, (int)sizeof(boolean_T));
      nm1d2 = fnew->size[0] * fnew->size[1];
      for (i0 = 0; i0 < nm1d2; i0++) {
        x->data[i0] = ((fnew->data[i0] > MinFreqRange) && (fnew->data[i0] <
          MaxFreqRange));
      }

      emxInit_int32_T(&ii, 2);
      idx = 0;
      i0 = ii->size[0] * ii->size[1];
      ii->size[0] = 1;
      ii->size[1] = x->size[1];
      emxEnsureCapacity((emxArray__common *)ii, i0, (int)sizeof(int));
      nm1d2 = 1;
      exitg1 = false;
      while ((!exitg1) && (nm1d2 <= x->size[1])) {
        guard1 = false;
        if (x->data[nm1d2 - 1]) {
          idx++;
          ii->data[idx - 1] = nm1d2;
          if (idx >= x->size[1]) {
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
        if (idx == 0) {
          i0 = ii->size[0] * ii->size[1];
          ii->size[0] = 1;
          ii->size[1] = 0;
          emxEnsureCapacity((emxArray__common *)ii, i0, (int)sizeof(int));
        }
      } else {
        if (1 > idx) {
          nm1d2 = 0;
        } else {
          nm1d2 = idx;
        }

        emxInit_int32_T(&b_ii, 2);
        i0 = b_ii->size[0] * b_ii->size[1];
        b_ii->size[0] = 1;
        b_ii->size[1] = nm1d2;
        emxEnsureCapacity((emxArray__common *)b_ii, i0, (int)sizeof(int));
        for (i0 = 0; i0 < nm1d2; i0++) {
          b_ii->data[b_ii->size[0] * i0] = ii->data[i0];
        }

        i0 = ii->size[0] * ii->size[1];
        ii->size[0] = 1;
        ii->size[1] = b_ii->size[1];
        emxEnsureCapacity((emxArray__common *)ii, i0, (int)sizeof(int));
        nm1d2 = b_ii->size[1];
        for (i0 = 0; i0 < nm1d2; i0++) {
          ii->data[ii->size[0] * i0] = b_ii->data[b_ii->size[0] * i0];
        }

        emxFree_int32_T(&b_ii);
      }

      emxFree_boolean_T(&x);
      emxInit_real_T(&idxG, 2);
      i0 = idxG->size[0] * idxG->size[1];
      idxG->size[0] = 1;
      idxG->size[1] = ii->size[1];
      emxEnsureCapacity((emxArray__common *)idxG, i0, (int)sizeof(double));
      nm1d2 = ii->size[0] * ii->size[1];
      for (i0 = 0; i0 < nm1d2; i0++) {
        idxG->data[i0] = ii->data[i0];
      }

      emxFree_int32_T(&ii);
      if (idxG->size[1] == 0) {
        *ErrCode = -1.0;
      } else if (idxG->size[1] > fnew->size[1]) {
        *ErrCode = -2.0;
      } else {
        emxInit_real_T(&xx, 2);
        i0 = xx->size[0] * xx->size[1];
        xx->size[0] = 1;
        emxEnsureCapacity((emxArray__common *)xx, i0, (int)sizeof(double));
        nm1d2 = idxG->size[1];
        i0 = xx->size[0] * xx->size[1];
        xx->size[1] = nm1d2;
        emxEnsureCapacity((emxArray__common *)xx, i0, (int)sizeof(double));
        nm1d2 = idxG->size[1];
        for (i0 = 0; i0 < nm1d2; i0++) {
          xx->data[i0] = 0.0;
        }

        /* %init for codegen; */
        if (idxG->data[0] > idxG->data[idxG->size[1] - 1]) {
          i0 = 1;
          idx = 0;
        } else {
          i0 = (int)idxG->data[0];
          idx = (int)idxG->data[idxG->size[1] - 1];
        }

        k = r0->size[0] * r0->size[1];
        r0->size[0] = 1;
        r0->size[1] = (idx - i0) + 1;
        emxEnsureCapacity((emxArray__common *)r0, k, (int)sizeof(double));
        nm1d2 = idx - i0;
        for (k = 0; k <= nm1d2; k++) {
          r0->data[r0->size[0] * k] = fnew->data[(i0 + k) - 1];
        }

        i0 = idx - i0;
        for (k = 0; k <= i0; k++) {
          r0->data[k] = log10(r0->data[k]);
        }

        nm1d2 = r0->size[1];
        for (i0 = 0; i0 < nm1d2; i0++) {
          xx->data[i0] = r0->data[r0->size[0] * i0];
        }

        nm1d2 = idxG->size[1];
        if (idxG->data[0] > idxG->data[idxG->size[1] - 1]) {
          i0 = 0;
          idx = 0;
        } else {
          i0 = (int)idxG->data[0] - 1;
          idx = (int)idxG->data[idxG->size[1] - 1];
        }

        emxInit_real_T(&b_xx, 2);
        k = b_xx->size[0] * b_xx->size[1];
        b_xx->size[0] = 1;
        b_xx->size[1] = nm1d2;
        emxEnsureCapacity((emxArray__common *)b_xx, k, (int)sizeof(double));
        for (k = 0; k < nm1d2; k++) {
          b_xx->data[b_xx->size[0] * k] = xx->data[k];
        }

        emxInit_real_T(&b_rdb, 2);
        k = b_rdb->size[0] * b_rdb->size[1];
        b_rdb->size[0] = 1;
        b_rdb->size[1] = idx - i0;
        emxEnsureCapacity((emxArray__common *)b_rdb, k, (int)sizeof(double));
        nm1d2 = idx - i0;
        for (idx = 0; idx < nm1d2; idx++) {
          b_rdb->data[b_rdb->size[0] * idx] = rdb->data[i0 + idx];
        }

        polyfit(b_xx, b_rdb, pCoefficient, expl_temp_data, expl_temp_size,
                &cdiff, &ndbl);
        *SlopeResult = pCoefficient[0];

        /* PlantGain22 = (2*pi)^2*10^(pCoefficient(2)/20) */
        /* %%for  plant gain calculation; */
        emxFree_real_T(&b_rdb);
        emxFree_real_T(&b_xx);
        memset(&rnondb[0], 0, 1000000U * sizeof(double));
        nm1d2 = fnew->size[1];
        i0 = xx->size[0] * xx->size[1];
        xx->size[0] = 1;
        xx->size[1] = nm1d2;
        emxEnsureCapacity((emxArray__common *)xx, i0, (int)sizeof(double));
        for (i0 = 0; i0 < nm1d2; i0++) {
          xx->data[xx->size[0] * i0] = rdb->data[i0] / 20.0;
        }

        for (i0 = 0; i0 < 2; i0++) {
          pCoefficient[i0] = xx->size[i0];
        }

        i0 = r0->size[0] * r0->size[1];
        r0->size[0] = 1;
        r0->size[1] = (int)pCoefficient[1];
        emxEnsureCapacity((emxArray__common *)r0, i0, (int)sizeof(double));
        for (k = 0; k < (int)pCoefficient[1]; k++) {
          r0->data[k] = rt_powd_snf(10.0, xx->data[k]);
        }

        emxFree_real_T(&xx);
        nm1d2 = r0->size[1];
        for (i0 = 0; i0 < nm1d2; i0++) {
          rnondb[i0] = r0->data[r0->size[0] * i0];
        }

        /* %non-dB */
        memset(&gainLineRange[0], 0, 1000000U * sizeof(double));
        for (nm1d2 = 0; nm1d2 < idxG->size[1]; nm1d2++) {
          gainLineRange[nm1d2] = rnondb[(int)idxG->data[nm1d2] - 1] * fnew->
            data[(int)idxG->data[nm1d2] - 1] * 2.0 * 3.1415926535897931 *
            fnew->data[(int)idxG->data[nm1d2] - 1] * 2.0 * 3.1415926535897931;
        }

        /* %can not use mean of matlab for the zero value  */
        *PlantGain = sum(gainLineRange) / (double)idxG->size[1];
      }

      emxFree_real_T(&r0);
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
