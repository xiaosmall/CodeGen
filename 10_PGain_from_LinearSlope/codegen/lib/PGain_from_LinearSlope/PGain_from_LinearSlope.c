/*
 * File: PGain_from_LinearSlope.c
 *
 * MATLAB Coder version            : 3.1
 * C/C++ source code generated on  : 20-Dec-2017 17:38:07
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "PGain_from_LinearSlope.h"
#include "PGain_from_LinearSlope_emxutil.h"
#include "mldivide.h"
#include "sum.h"
#include "mean.h"
#include "all.h"
#include "abs.h"
#include "fix.h"
#include "rdivide.h"
#include "nullAssignment.h"
#include "sort1.h"
#include "diff.h"
#include "log10.h"
#include "interp1.h"
#include "unique.h"
#include "sortIdx.h"

/* Function Declarations */
static double rt_powd_snf(double u0, double u1);
static double rt_roundd_snf(double u);

/* Function Definitions */

/*
 * Arguments    : double u0
 *                double u1
 * Return Type  : double
 */
static double rt_powd_snf(double u0, double u1)
{
  double y;
  double d1;
  double d2;
  if (rtIsNaN(u0) || rtIsNaN(u1)) {
    y = rtNaN;
  } else {
    d1 = fabs(u0);
    d2 = fabs(u1);
    if (rtIsInf(u1)) {
      if (d1 == 1.0) {
        y = rtNaN;
      } else if (d1 > 1.0) {
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
    } else if (d2 == 0.0) {
      y = 1.0;
    } else if (d2 == 1.0) {
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
 * Arguments    : double u
 * Return Type  : double
 */
static double rt_roundd_snf(double u)
{
  double y;
  if (fabs(u) < 4.503599627370496E+15) {
    if (u >= 0.5) {
      y = floor(u + 0.5);
    } else if (u > -0.5) {
      y = u * 0.0;
    } else {
      y = ceil(u - 0.5);
    }
  } else {
    y = u;
  }

  return y;
}

/*
 * [StartFrequency, EndFrequency, SlopeResult, NumberOfPoints,PlantGain] = PGain_from_LinearSlope(f, r, n, MinFreqRange, MaxFreqRange, Slope, Tolerance)
 * Ver02 - add nErrCode, other bug is fixed at the LinearSlopeRange(Ver06)
 * Ver03 -  df chg to 0.05 from 0.01,otherwise size too large;
 *          MaxMemorySize is chged from 1000 to 1e7;
 * Arguments    : const double inf[1000]
 *                const double r[1000]
 *                double n
 *                double MinFreqRange
 *                double MaxFreqRange
 *                double Slope
 *                double Tolerance
 *                const double vPh[1000]
 *                double *StartFrequency
 *                double *EndFrequency
 *                double *SlopeResult
 *                double *NumberOfPoints
 *                double *PlantGain
 *                double *nErrCode
 * Return Type  : void
 */
void PGain_from_LinearSlope(const double inf[1000], const double r[1000], double
  n, double MinFreqRange, double MaxFreqRange, double Slope, double Tolerance,
  const double vPh[1000], double *StartFrequency, double *EndFrequency, double
  *SlopeResult, double *NumberOfPoints, double *PlantGain, double *nErrCode)
{
  static double gainLineRange[10000000];
  int loop_ub;
  int b_n;
  double anew;
  double apnd;
  double ndbl;
  emxArray_real_T *fnew;
  double cdiff;
  double absa;
  int i0;
  double absb;
  double inf_data[1000];
  int nm1d2;
  int ii_data[1];
  int k;
  int idx_data[1000];
  int idx_size[1];
  double kd;
  int khi;
  int idx;
  int nb;
  int exitg11;
  int pos_data[1000];
  int exponent;
  boolean_T p;
  int index_f_afterUnique_data[1000];
  double b_data[1000];
  int b_size[1];
  double r_data[1000];
  int r_size[1];
  emxArray_real_T *rdb;
  int b_fnew[1];
  emxArray_real_T c_fnew;
  int b_b_size[1];
  int vPh_size[1];
  emxArray_real_T *phdeg;
  emxArray_boolean_T *x;
  int d_fnew[1];
  emxArray_int32_T *ii;
  boolean_T exitg10;
  boolean_T guard5 = false;
  emxArray_real_T *idxPhLarger180;
  emxArray_real_T *b_log10f;
  double tol;
  boolean_T exitg9;
  double MinIndexRange;
  emxArray_real_T *y;
  boolean_T exitg8;
  boolean_T guard4 = false;
  double d0;
  emxArray_real_T *c_log10f;
  emxArray_real_T *b_rdb;
  emxArray_real_T *b_y;
  emxArray_real_T *varargin_2;
  double b_x[2];
  boolean_T exitg7;
  boolean_T guard3 = false;
  emxArray_real_T *kk;
  boolean_T exitg6;
  boolean_T guard2 = false;
  emxArray_real_T *c_y;
  int exitg5;
  emxArray_real_T *b_kk;
  emxArray_real_T *c_kk;
  double i2mean;
  boolean_T exitg4;
  double l[3];
  boolean_T exitg3;
  emxArray_real_T *IndexRange;
  emxArray_real_T *d_log10f;
  emxArray_real_T *c_rdb;
  emxArray_real_T *b_phdeg;
  boolean_T exitg2;
  unsigned int k_data[1];
  boolean_T flag_data[1];
  double StartFrequency_data[1];
  double k1_data[1];
  double EndFrequency_data[1];
  double SlopeResult_data[1];
  double k2_data[1];
  boolean_T exitg1;
  boolean_T guard1 = false;

  /* %nErrCode */
  /*  -1 : input parameter wrong: MinFreqRange or MaxFreqRange outof range; */
  /*  -2 : input parameter (MinFreqRange> MaxFreqRange) conflict; */
  /*  -3: not find any period within the tolorance; */
  /*  -5; calculation fail; */
  /*  -10: calculate failed; */
  /*  -11: freq range is out of middle variable size; */
  /* to use phase info to limit the slope searching range, any freq whose phase */
  /* smaller than -240deg, will be not the linear slope range ,will not be the */
  /* search range: */
  /*  idxPhLarger180 = find( phdeg <=( -180-60 )); */
  /*  EndIndex = idxPhLarger180(1); */
  /* %% */
  /*  8 input: */
  /*  inf, r, n, MinFreqRange, MaxFreqRange, Slope, Tolerance,vPh */
  /*  6 output */
  /*  % StartFrequency, EndFrequency, SlopeResult, NumberOfPoints,PlantGain, nErrCode */
  /* %%%%%%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
  *StartFrequency = 0.0;

  /*  to be ready if we return due to some error */
  *EndFrequency = 0.0;
  *SlopeResult = 0.0;
  *NumberOfPoints = 0.0;
  *PlantGain = 0.0;
  *nErrCode = 0.0;
  memset(&gainLineRange[0], 0, 10000000U * sizeof(double));

  /* %1000 is for the c codegen,since c code will assign other n number of point to other non-zero value, this will affect to the following sum calculation */
  /*      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%input value protection; */
  /* -------------------------------------------------------------------------% */
  if (1.0 > n) {
    loop_ub = 0;
  } else {
    loop_ub = (int)n;
  }

  if (MinFreqRange > MaxFreqRange) {
    *nErrCode = -2.0;
  }

  if (*nErrCode < 0.0) {
  } else {
    if ((MinFreqRange < inf[0]) || (MaxFreqRange > inf[loop_ub - 1])) {
      *nErrCode = -1.0;
    }

    if (*nErrCode < 0.0) {
    } else {
      /* -------------------------------------------------------------------------% */
      /*  Generate fix resolution response vector (0.1Hz) */
      /*  */
      if (rtIsNaN(inf[0]) || rtIsNaN(inf[(int)n - 1])) {
        b_n = 1;
        anew = rtNaN;
        apnd = inf[(int)n - 1];
      } else if (inf[(int)n - 1] < inf[0]) {
        b_n = 0;
        anew = inf[0];
        apnd = inf[(int)n - 1];
      } else if (rtIsInf(inf[0]) || rtIsInf(inf[(int)n - 1])) {
        b_n = 1;
        anew = rtNaN;
        apnd = inf[(int)n - 1];
      } else {
        anew = inf[0];
        ndbl = floor((inf[(int)n - 1] - inf[0]) / 0.05 + 0.5);
        apnd = inf[0] + ndbl * 0.05;
        cdiff = apnd - inf[(int)n - 1];
        absa = fabs(inf[0]);
        absb = fabs(inf[(int)n - 1]);
        if ((absa >= absb) || rtIsNaN(absb)) {
          absb = absa;
        }

        if (fabs(cdiff) < 4.4408920985006262E-16 * absb) {
          ndbl++;
          apnd = inf[(int)n - 1];
        } else if (cdiff > 0.0) {
          apnd = inf[0] + (ndbl - 1.0) * 0.05;
        } else {
          ndbl++;
        }

        if (ndbl >= 0.0) {
          b_n = (int)ndbl;
        } else {
          b_n = 0;
        }
      }

      emxInit_real_T(&fnew, 2);
      i0 = fnew->size[0] * fnew->size[1];
      fnew->size[0] = 1;
      fnew->size[1] = b_n;
      emxEnsureCapacity((emxArray__common *)fnew, i0, (int)sizeof(double));
      if (b_n > 0) {
        fnew->data[0] = anew;
        if (b_n > 1) {
          fnew->data[b_n - 1] = apnd;
          nm1d2 = (b_n - 1) / 2;
          for (k = 1; k < nm1d2; k++) {
            kd = (double)k * 0.05;
            fnew->data[k] = anew + kd;
            fnew->data[(b_n - k) - 1] = apnd - kd;
          }

          if (nm1d2 << 1 == b_n - 1) {
            fnew->data[nm1d2] = (anew + apnd) / 2.0;
          } else {
            kd = (double)nm1d2 * 0.05;
            fnew->data[nm1d2] = anew + kd;
            fnew->data[nm1d2 + 1] = apnd - kd;
          }
        }
      }

      /*  rdb = interp1(f(1:n),r(1:n),fnew); */
      /* to solve f is non-unique problem; */
      ii_data[0] = loop_ub;
      for (i0 = 0; i0 < loop_ub; i0++) {
        inf_data[i0] = inf[i0];
      }

      sortIdx(inf_data, ii_data, idx_data, idx_size);
      for (k = 0; k + 1 <= loop_ub; k++) {
        inf_data[k] = inf[idx_data[k] - 1];
      }

      count_nonfinites(inf_data, loop_ub, &k, &khi, &idx, &b_n);
      nb = -1;
      if (k > 0) {
        nb = 0;
      }

      khi += k;
      while (k + 1 <= khi) {
        ndbl = inf_data[k];
        nm1d2 = k;
        do {
          exitg11 = 0;
          k++;
          if (k + 1 > khi) {
            exitg11 = 1;
          } else {
            kd = fabs(ndbl / 2.0);
            if ((!rtIsInf(kd)) && (!rtIsNaN(kd))) {
              if (kd <= 2.2250738585072014E-308) {
                kd = 4.94065645841247E-324;
              } else {
                frexp(kd, &exponent);
                kd = ldexp(1.0, exponent - 53);
              }
            } else {
              kd = rtNaN;
            }

            if ((fabs(ndbl - inf_data[k]) < kd) || (rtIsInf(inf_data[k]) &&
                 rtIsInf(ndbl) && ((inf_data[k] > 0.0) == (ndbl > 0.0)))) {
              p = true;
            } else {
              p = false;
            }

            if (!p) {
              exitg11 = 1;
            }
          }
        } while (exitg11 == 0);

        nb++;
        inf_data[nb] = ndbl;
        idx_data[nb] = idx_data[nm1d2];
      }

      if (idx > 0) {
        nb++;
        inf_data[nb] = inf_data[khi];
        idx_data[nb] = idx_data[khi];
      }

      k = (khi + idx) - 1;
      for (khi = 1; khi <= b_n; khi++) {
        nb++;
        inf_data[nb] = inf_data[k + khi];
        idx_data[nb] = idx_data[k + khi];
      }

      for (k = 0; k + 1 <= nb + 1; k++) {
        pos_data[k] = idx_data[k];
      }

      loop_ub = nb + 1;
      for (i0 = 0; i0 < loop_ub; i0++) {
        index_f_afterUnique_data[i0] = pos_data[i0];
      }

      /* sort and remove duplicates */
      if (1 > nb + 1) {
        loop_ub = 0;
      } else {
        loop_ub = nb + 1;
      }

      if (1 > nb + 1) {
        khi = 0;
      } else {
        khi = nb + 1;
      }

      b_size[0] = loop_ub;
      for (i0 = 0; i0 < loop_ub; i0++) {
        b_data[i0] = inf_data[i0];
      }

      r_size[0] = khi;
      for (i0 = 0; i0 < khi; i0++) {
        r_data[i0] = r[index_f_afterUnique_data[i0] - 1];
      }

      emxInit_real_T1(&rdb, 1);
      b_fnew[0] = fnew->size[1];
      c_fnew = *fnew;
      c_fnew.size = (int *)&b_fnew;
      c_fnew.numDimensions = 1;
      interp1(b_data, b_size, r_data, r_size, &c_fnew, rdb);

      /* EndIndex = length(fnew); */
      /* improve @20160728 */
      if (1 > nb + 1) {
        loop_ub = 0;
      } else {
        loop_ub = nb + 1;
      }

      if (1 > nb + 1) {
        khi = 0;
      } else {
        khi = nb + 1;
      }

      b_b_size[0] = loop_ub;
      for (i0 = 0; i0 < loop_ub; i0++) {
        b_data[i0] = inf_data[i0];
      }

      vPh_size[0] = khi;
      for (i0 = 0; i0 < khi; i0++) {
        inf_data[i0] = vPh[i0];
      }

      emxInit_real_T1(&phdeg, 1);
      emxInit_boolean_T(&x, 1);
      d_fnew[0] = fnew->size[1];
      c_fnew = *fnew;
      c_fnew.size = (int *)&d_fnew;
      c_fnew.numDimensions = 1;
      interp1(b_data, b_b_size, inf_data, vPh_size, &c_fnew, phdeg);

      /* %improve program,add vPh info */
      /* %%not use the ending freq;but use the vPh=-180deg to set the EndIndex; */
      i0 = x->size[0];
      x->size[0] = phdeg->size[0];
      emxEnsureCapacity((emxArray__common *)x, i0, (int)sizeof(boolean_T));
      loop_ub = phdeg->size[0];
      for (i0 = 0; i0 < loop_ub; i0++) {
        x->data[i0] = (phdeg->data[i0] <= -240.0);
      }

      emxInit_int32_T(&ii, 1);
      khi = x->size[0];
      idx = 0;
      i0 = ii->size[0];
      ii->size[0] = x->size[0];
      emxEnsureCapacity((emxArray__common *)ii, i0, (int)sizeof(int));
      b_n = 1;
      exitg10 = false;
      while ((!exitg10) && (b_n <= khi)) {
        guard5 = false;
        if (x->data[b_n - 1]) {
          idx++;
          ii->data[idx - 1] = b_n;
          if (idx >= khi) {
            exitg10 = true;
          } else {
            guard5 = true;
          }
        } else {
          guard5 = true;
        }

        if (guard5) {
          b_n++;
        }
      }

      if (x->size[0] == 1) {
        if (idx == 0) {
          i0 = ii->size[0];
          ii->size[0] = 0;
          emxEnsureCapacity((emxArray__common *)ii, i0, (int)sizeof(int));
        }
      } else {
        i0 = ii->size[0];
        if (1 > idx) {
          ii->size[0] = 0;
        } else {
          ii->size[0] = idx;
        }

        emxEnsureCapacity((emxArray__common *)ii, i0, (int)sizeof(int));
      }

      emxInit_real_T1(&idxPhLarger180, 1);
      i0 = idxPhLarger180->size[0];
      idxPhLarger180->size[0] = ii->size[0];
      emxEnsureCapacity((emxArray__common *)idxPhLarger180, i0, (int)sizeof
                        (double));
      loop_ub = ii->size[0];
      for (i0 = 0; i0 < loop_ub; i0++) {
        idxPhLarger180->data[i0] = ii->data[i0];
      }

      if (idxPhLarger180->size[0] == 0) {
        nm1d2 = 1;
        b_n = phdeg->size[0];
        ndbl = phdeg->data[0];
        if (phdeg->size[0] > 1) {
          if (rtIsNaN(phdeg->data[0])) {
            idx = 2;
            exitg9 = false;
            while ((!exitg9) && (idx <= b_n)) {
              nm1d2 = idx;
              if (!rtIsNaN(phdeg->data[idx - 1])) {
                ndbl = phdeg->data[idx - 1];
                exitg9 = true;
              } else {
                idx++;
              }
            }
          }

          if (nm1d2 < phdeg->size[0]) {
            while (nm1d2 + 1 <= b_n) {
              if (phdeg->data[nm1d2] < ndbl) {
                ndbl = phdeg->data[nm1d2];
              }

              nm1d2++;
            }
          }
        }

        i0 = x->size[0];
        x->size[0] = phdeg->size[0];
        emxEnsureCapacity((emxArray__common *)x, i0, (int)sizeof(boolean_T));
        loop_ub = phdeg->size[0];
        for (i0 = 0; i0 < loop_ub; i0++) {
          x->data[i0] = (phdeg->data[i0] == ndbl);
        }

        khi = x->size[0];
        idx = 0;
        i0 = ii->size[0];
        ii->size[0] = x->size[0];
        emxEnsureCapacity((emxArray__common *)ii, i0, (int)sizeof(int));
        b_n = 1;
        exitg8 = false;
        while ((!exitg8) && (b_n <= khi)) {
          guard4 = false;
          if (x->data[b_n - 1]) {
            idx++;
            ii->data[idx - 1] = b_n;
            if (idx >= khi) {
              exitg8 = true;
            } else {
              guard4 = true;
            }
          } else {
            guard4 = true;
          }

          if (guard4) {
            b_n++;
          }
        }

        if (x->size[0] == 1) {
          if (idx == 0) {
            i0 = ii->size[0];
            ii->size[0] = 0;
            emxEnsureCapacity((emxArray__common *)ii, i0, (int)sizeof(int));
          }
        } else {
          i0 = ii->size[0];
          if (1 > idx) {
            ii->size[0] = 0;
          } else {
            ii->size[0] = idx;
          }

          emxEnsureCapacity((emxArray__common *)ii, i0, (int)sizeof(int));
        }

        i0 = idxPhLarger180->size[0];
        idxPhLarger180->size[0] = ii->size[0];
        emxEnsureCapacity((emxArray__common *)idxPhLarger180, i0, (int)sizeof
                          (double));
        loop_ub = ii->size[0];
        for (i0 = 0; i0 < loop_ub; i0++) {
          idxPhLarger180->data[i0] = ii->data[i0];
        }
      }

      emxInit_real_T1(&b_log10f, 1);

      /*  */
      /*  */
      tol = fabs(Tolerance * Slope);
      i0 = b_log10f->size[0];
      b_log10f->size[0] = fnew->size[1];
      emxEnsureCapacity((emxArray__common *)b_log10f, i0, (int)sizeof(double));
      loop_ub = fnew->size[1];
      for (i0 = 0; i0 < loop_ub; i0++) {
        b_log10f->data[i0] = fnew->data[i0];
      }

      b_log10(b_log10f);
      MinIndexRange = MinFreqRange / 0.05;
      b_fix(&MinIndexRange);

      /*  */
      /*  Must be at least the minimal number of elements. */
      /*  */
      if (MinFreqRange > inf[nb] - inf[0]) {
        *nErrCode = -5.0;

        /* %??? */
      } else {
        emxInit_real_T1(&y, 1);

        /*  */
        i0 = y->size[0];
        y->size[0] = (int)idxPhLarger180->data[0];
        emxEnsureCapacity((emxArray__common *)y, i0, (int)sizeof(double));
        loop_ub = (int)idxPhLarger180->data[0];
        for (i0 = 0; i0 < loop_ub; i0++) {
          y->data[i0] = 0.0;
        }

        /*  */
        /*  Claculate the slopes of all the gain responses in MinIndexRange window */
        /*  and put the results in y vector. */
        /*  */
        d0 = idxPhLarger180->data[0] - MinIndexRange;
        k = 0;
        emxInit_real_T(&c_log10f, 2);
        emxInit_real_T1(&b_rdb, 1);
        while (k <= (int)d0 - 1) {
          kd = ((1.0 + (double)k) + MinIndexRange) - 1.0;
          if (1.0 + (double)k > kd) {
            i0 = 1;
            nm1d2 = 0;
          } else {
            i0 = k + 1;
            nm1d2 = (int)kd;
          }

          kd = ((1.0 + (double)k) + MinIndexRange) - 1.0;
          if (1.0 + (double)k > kd) {
            idx = 1;
            b_n = 0;
          } else {
            idx = k + 1;
            b_n = (int)kd;
          }

          khi = c_log10f->size[0] * c_log10f->size[1];
          c_log10f->size[0] = (nm1d2 - i0) + 1;
          c_log10f->size[1] = 2;
          emxEnsureCapacity((emxArray__common *)c_log10f, khi, (int)sizeof
                            (double));
          loop_ub = (nm1d2 - i0) + 1;
          for (nm1d2 = 0; nm1d2 < loop_ub; nm1d2++) {
            c_log10f->data[nm1d2] = b_log10f->data[(i0 + nm1d2) - 1];
          }

          loop_ub = (int)MinIndexRange;
          for (i0 = 0; i0 < loop_ub; i0++) {
            c_log10f->data[i0 + c_log10f->size[0]] = 1.0;
          }

          i0 = b_rdb->size[0];
          b_rdb->size[0] = (b_n - idx) + 1;
          emxEnsureCapacity((emxArray__common *)b_rdb, i0, (int)sizeof(double));
          loop_ub = (b_n - idx) + 1;
          for (i0 = 0; i0 < loop_ub; i0++) {
            b_rdb->data[i0] = rdb->data[(idx + i0) - 1];
          }

          mldivide(c_log10f, b_rdb, b_x);
          y->data[k] = b_x[0];
          k++;
        }

        emxFree_real_T(&b_rdb);
        emxFree_real_T(&c_log10f);
        emxInit_real_T1(&b_y, 1);

        /*  */
        /*  Calculate the longest series of successive */
        /*  indexes where the slope meets the required */
        /*  accuracy. */
        /*  i2mean is index into the mean value of these */
        /*  series */
        /*  */
        /*  kk - contains the indexes of all the slopes that meet the required */
        /*  accuracy */
        /*  */
        i0 = b_y->size[0];
        b_y->size[0] = y->size[0];
        emxEnsureCapacity((emxArray__common *)b_y, i0, (int)sizeof(double));
        loop_ub = y->size[0];
        for (i0 = 0; i0 < loop_ub; i0++) {
          b_y->data[i0] = y->data[i0] - Slope;
        }

        emxInit_real_T1(&varargin_2, 1);
        b_abs(b_y, varargin_2);
        i0 = x->size[0];
        x->size[0] = varargin_2->size[0];
        emxEnsureCapacity((emxArray__common *)x, i0, (int)sizeof(boolean_T));
        loop_ub = varargin_2->size[0];
        emxFree_real_T(&b_y);
        for (i0 = 0; i0 < loop_ub; i0++) {
          x->data[i0] = (varargin_2->data[i0] < tol);
        }

        khi = x->size[0];
        idx = 0;
        i0 = ii->size[0];
        ii->size[0] = x->size[0];
        emxEnsureCapacity((emxArray__common *)ii, i0, (int)sizeof(int));
        b_n = 1;
        exitg7 = false;
        while ((!exitg7) && (b_n <= khi)) {
          guard3 = false;
          if (x->data[b_n - 1]) {
            idx++;
            ii->data[idx - 1] = b_n;
            if (idx >= khi) {
              exitg7 = true;
            } else {
              guard3 = true;
            }
          } else {
            guard3 = true;
          }

          if (guard3) {
            b_n++;
          }
        }

        if (x->size[0] == 1) {
          if (idx == 0) {
            i0 = ii->size[0];
            ii->size[0] = 0;
            emxEnsureCapacity((emxArray__common *)ii, i0, (int)sizeof(int));
          }
        } else {
          i0 = ii->size[0];
          if (1 > idx) {
            ii->size[0] = 0;
          } else {
            ii->size[0] = idx;
          }

          emxEnsureCapacity((emxArray__common *)ii, i0, (int)sizeof(int));
        }

        emxInit_real_T1(&kk, 1);
        i0 = kk->size[0];
        kk->size[0] = ii->size[0];
        emxEnsureCapacity((emxArray__common *)kk, i0, (int)sizeof(double));
        loop_ub = ii->size[0];
        for (i0 = 0; i0 < loop_ub; i0++) {
          kk->data[i0] = ii->data[i0];
        }

        diff(kk, varargin_2);
        i0 = x->size[0];
        x->size[0] = varargin_2->size[0];
        emxEnsureCapacity((emxArray__common *)x, i0, (int)sizeof(boolean_T));
        loop_ub = varargin_2->size[0];
        for (i0 = 0; i0 < loop_ub; i0++) {
          x->data[i0] = (varargin_2->data[i0] == 1.0);
        }

        khi = x->size[0] - 1;
        b_n = 0;
        for (nm1d2 = 0; nm1d2 <= khi; nm1d2++) {
          if (x->data[nm1d2]) {
            b_n++;
          }
        }

        idx = 0;
        for (nm1d2 = 0; nm1d2 <= khi; nm1d2++) {
          if (x->data[nm1d2]) {
            kk->data[idx] = kk->data[nm1d2];
            idx++;
          }
        }

        i0 = kk->size[0];
        kk->size[0] = b_n;
        emxEnsureCapacity((emxArray__common *)kk, i0, (int)sizeof(double));
        if (kk->size[0] == 0) {
          *nErrCode = -3.0;

          /* %tolerance too tighten,pls increase tolerance */
        } else {
          /*  */
          /*  jj - is set of indexes of kk. kk(jj(1)), kk(jj(2), ..., kk(jj(N) */
          /*  are the indexes that point to the starting and the ending of sequence of kk. */
          diff(kk, varargin_2);
          i0 = x->size[0];
          x->size[0] = varargin_2->size[0];
          emxEnsureCapacity((emxArray__common *)x, i0, (int)sizeof(boolean_T));
          loop_ub = varargin_2->size[0];
          for (i0 = 0; i0 < loop_ub; i0++) {
            x->data[i0] = (varargin_2->data[i0] != 1.0);
          }

          khi = x->size[0];
          idx = 0;
          i0 = ii->size[0];
          ii->size[0] = x->size[0];
          emxEnsureCapacity((emxArray__common *)ii, i0, (int)sizeof(int));
          b_n = 1;
          exitg6 = false;
          while ((!exitg6) && (b_n <= khi)) {
            guard2 = false;
            if (x->data[b_n - 1]) {
              idx++;
              ii->data[idx - 1] = b_n;
              if (idx >= khi) {
                exitg6 = true;
              } else {
                guard2 = true;
              }
            } else {
              guard2 = true;
            }

            if (guard2) {
              b_n++;
            }
          }

          if (x->size[0] == 1) {
            if (idx == 0) {
              i0 = ii->size[0];
              ii->size[0] = 0;
              emxEnsureCapacity((emxArray__common *)ii, i0, (int)sizeof(int));
            }
          } else {
            i0 = ii->size[0];
            if (1 > idx) {
              ii->size[0] = 0;
            } else {
              ii->size[0] = idx;
            }

            emxEnsureCapacity((emxArray__common *)ii, i0, (int)sizeof(int));
          }

          i0 = phdeg->size[0];
          phdeg->size[0] = ii->size[0];
          emxEnsureCapacity((emxArray__common *)phdeg, i0, (int)sizeof(double));
          loop_ub = ii->size[0];
          for (i0 = 0; i0 < loop_ub; i0++) {
            phdeg->data[i0] = ii->data[i0];
          }

          /* %modify @20160728 */
          if (phdeg->size[0] == 0) {
            emxInit_real_T1(&c_y, 1);

            /* means only one stage inside of the tolorance */
            *StartFrequency = fnew->data[(int)kk->data[0] - 1];
            *EndFrequency = fnew->data[(int)kk->data[kk->size[0] - 1] - 1];
            loop_ub = kk->size[0];
            i0 = c_y->size[0];
            c_y->size[0] = loop_ub;
            emxEnsureCapacity((emxArray__common *)c_y, i0, (int)sizeof(double));
            for (i0 = 0; i0 < loop_ub; i0++) {
              c_y->data[i0] = y->data[(int)kk->data[i0] - 1];
            }

            *SlopeResult = mean(c_y);
            *NumberOfPoints = kk->size[0];

            /* means success; */
            /* %if return here,calculate the plantGain;!!!!! */
            b_n = 0;
            emxFree_real_T(&c_y);
            do {
              exitg5 = 0;
              i0 = kk->size[0];
              if (b_n <= i0 - 1) {
                gainLineRange[b_n] = rt_powd_snf(10.0, rdb->data[(int)kk->
                  data[b_n] - 1] / 20.0) * fnew->data[(int)kk->data[b_n] - 1] *
                  2.0 * 3.1415926535897931 * fnew->data[(int)kk->data[b_n] - 1] *
                  2.0 * 3.1415926535897931;
                b_n++;
              } else {
                exitg5 = 1;
              }
            } while (exitg5 == 0);

            *PlantGain = sum(gainLineRange) / (double)kk->size[0];

            /* %can not use mean of matlab for the zero value */
            /* %%finish cal plant gain */
          } else {
            emxInit_real_T1(&c_y, 1);
            i0 = c_y->size[0];
            c_y->size[0] = (phdeg->size[0] + phdeg->size[0]) + 1;
            emxEnsureCapacity((emxArray__common *)c_y, i0, (int)sizeof(double));
            c_y->data[0] = 1.0;
            loop_ub = phdeg->size[0];
            for (i0 = 0; i0 < loop_ub; i0++) {
              c_y->data[i0 + 1] = phdeg->data[i0];
            }

            loop_ub = phdeg->size[0];
            for (i0 = 0; i0 < loop_ub; i0++) {
              c_y->data[(i0 + phdeg->size[0]) + 1] = phdeg->data[i0] + 1.0;
            }

            i0 = phdeg->size[0];
            phdeg->size[0] = c_y->size[0];
            emxEnsureCapacity((emxArray__common *)phdeg, i0, (int)sizeof(double));
            loop_ub = c_y->size[0];
            for (i0 = 0; i0 < loop_ub; i0++) {
              phdeg->data[i0] = c_y->data[i0];
            }

            emxFree_real_T(&c_y);
            emxInit_real_T1(&b_kk, 1);
            sort(phdeg);
            i0 = phdeg->size[0];
            nullAssignment(phdeg, i0);

            /*  */
            /*  ii - is an index of jj. kk(jj(ii)) is the index to starting point of the */
            /*  longest sequence. */
            /*  */
            /*  Note that we are looking for the longest sequence from the point of view */
            /*  of a user looking at the Bode graph. This means, when the X axis is */
            /*  logarithmic. */
            /*  As a result, we use the "/" (division) of end and start frequencies of */
            /*  each range (and not "-"). */
            /*  */
            i0 = phdeg->size[0] - 2;
            nm1d2 = b_kk->size[0];
            b_kk->size[0] = (i0 >> 1) + 1;
            emxEnsureCapacity((emxArray__common *)b_kk, nm1d2, (int)sizeof
                              (double));
            loop_ub = i0 >> 1;
            for (i0 = 0; i0 <= loop_ub; i0++) {
              b_kk->data[i0] = kk->data[(int)phdeg->data[1 + (i0 << 1)] - 1];
            }

            emxInit_real_T1(&c_kk, 1);
            i0 = (int)((double)phdeg->size[0] - 1.0) - 1;
            nm1d2 = c_kk->size[0];
            c_kk->size[0] = (i0 >> 1) + 1;
            emxEnsureCapacity((emxArray__common *)c_kk, nm1d2, (int)sizeof
                              (double));
            loop_ub = i0 >> 1;
            for (i0 = 0; i0 <= loop_ub; i0++) {
              c_kk->data[i0] = kk->data[(int)phdeg->data[i0 << 1] - 1];
            }

            rdivide(b_kk, c_kk, varargin_2);
            nm1d2 = 1;
            b_n = varargin_2->size[0];
            ndbl = varargin_2->data[0];
            khi = 0;
            emxFree_real_T(&c_kk);
            emxFree_real_T(&b_kk);
            if (varargin_2->size[0] > 1) {
              if (rtIsNaN(varargin_2->data[0])) {
                idx = 2;
                exitg4 = false;
                while ((!exitg4) && (idx <= b_n)) {
                  nm1d2 = idx;
                  if (!rtIsNaN(varargin_2->data[idx - 1])) {
                    ndbl = varargin_2->data[idx - 1];
                    khi = idx - 1;
                    exitg4 = true;
                  } else {
                    idx++;
                  }
                }
              }

              if (nm1d2 < varargin_2->size[0]) {
                while (nm1d2 + 1 <= b_n) {
                  if (varargin_2->data[nm1d2] > ndbl) {
                    ndbl = varargin_2->data[nm1d2];
                    khi = nm1d2;
                  }

                  nm1d2++;
                }
              }
            }

            /* %bug:if the tolorance too small,the ii maybe empty */
            i2mean = rt_roundd_snf(sqrt(kk->data[(int)phdeg->data[khi + 1] - 1] *
              kk->data[(int)phdeg->data[khi] - 1]));

            /*  Geometrical mean */
            /*  */
            /*  Calculate the permmited range that the min range can be extended to */
            /*  */
            l[0] = 2.0 * (i2mean - 1.0) + MinIndexRange;

            /*  limitation by the first index of the vectors. */
            l[1] = 2.0 * (idxPhLarger180->data[0] - ((i2mean - 1.0) +
              MinIndexRange)) + MinIndexRange;

            /*  limitation by the last index of the vectors. */
            l[2] = MaxFreqRange / 0.05;
            b_fix(&l[2]);

            /*  limitation as set by the user. */
            nm1d2 = 1;
            ndbl = l[0];
            if (rtIsNaN(l[0])) {
              idx = 2;
              exitg3 = false;
              while ((!exitg3) && (idx < 4)) {
                nm1d2 = idx;
                if (!rtIsNaN(l[idx - 1])) {
                  ndbl = l[idx - 1];
                  exitg3 = true;
                } else {
                  idx++;
                }
              }
            }

            if (nm1d2 < 3) {
              while (nm1d2 + 1 < 4) {
                if (l[nm1d2] < ndbl) {
                  ndbl = l[nm1d2];
                }

                nm1d2++;
              }
            }

            /*  */
            /*  Claculate the slopes of all gain responses around the optimal point */
            /*  and put the results in z vector. */
            /*  */
            kd = (ndbl - MinIndexRange) / 2.0;
            i0 = phdeg->size[0];
            phdeg->size[0] = (int)kd;
            emxEnsureCapacity((emxArray__common *)phdeg, i0, (int)sizeof(double));
            loop_ub = (int)kd;
            for (i0 = 0; i0 < loop_ub; i0++) {
              phdeg->data[i0] = 0.0;
            }

            d0 = (ndbl - MinIndexRange) / 2.0;
            k = 0;
            emxInit_real_T(&IndexRange, 2);
            emxInit_real_T(&d_log10f, 2);
            emxInit_real_T1(&c_rdb, 1);
            while (k <= (int)d0 - 1) {
              anew = i2mean - (1.0 + (double)k);
              kd = ((i2mean + MinIndexRange) + (1.0 + (double)k)) - 1.0;
              if (rtIsNaN(anew) || rtIsNaN(kd)) {
                b_n = 1;
                anew = rtNaN;
                apnd = kd;
              } else if (kd < anew) {
                b_n = 0;
                apnd = kd;
              } else if (rtIsInf(anew) || rtIsInf(kd)) {
                b_n = 1;
                anew = rtNaN;
                apnd = kd;
              } else {
                ndbl = floor((kd - anew) + 0.5);
                apnd = anew + ndbl;
                cdiff = apnd - kd;
                absa = fabs(anew);
                absb = fabs(kd);
                if ((absa >= absb) || rtIsNaN(absb)) {
                  absb = absa;
                }

                if (fabs(cdiff) < 4.4408920985006262E-16 * absb) {
                  ndbl++;
                  apnd = kd;
                } else if (cdiff > 0.0) {
                  apnd = anew + (ndbl - 1.0);
                } else {
                  ndbl++;
                }

                if (ndbl >= 0.0) {
                  b_n = (int)ndbl;
                } else {
                  b_n = 0;
                }
              }

              i0 = IndexRange->size[0] * IndexRange->size[1];
              IndexRange->size[0] = 1;
              IndexRange->size[1] = b_n;
              emxEnsureCapacity((emxArray__common *)IndexRange, i0, (int)sizeof
                                (double));
              if (b_n > 0) {
                IndexRange->data[0] = anew;
                if (b_n > 1) {
                  IndexRange->data[b_n - 1] = apnd;
                  nm1d2 = (b_n - 1) / 2;
                  for (khi = 1; khi < nm1d2; khi++) {
                    IndexRange->data[khi] = anew + (double)khi;
                    IndexRange->data[(b_n - khi) - 1] = apnd - (double)khi;
                  }

                  if (nm1d2 << 1 == b_n - 1) {
                    IndexRange->data[nm1d2] = (anew + apnd) / 2.0;
                  } else {
                    IndexRange->data[nm1d2] = anew + (double)nm1d2;
                    IndexRange->data[nm1d2 + 1] = apnd - (double)nm1d2;
                  }
                }
              }

              khi = IndexRange->size[1];
              i0 = d_log10f->size[0] * d_log10f->size[1];
              d_log10f->size[0] = IndexRange->size[1];
              d_log10f->size[1] = 2;
              emxEnsureCapacity((emxArray__common *)d_log10f, i0, (int)sizeof
                                (double));
              loop_ub = IndexRange->size[1];
              for (i0 = 0; i0 < loop_ub; i0++) {
                d_log10f->data[i0] = b_log10f->data[(int)IndexRange->
                  data[IndexRange->size[0] * i0] - 1];
              }

              for (i0 = 0; i0 < khi; i0++) {
                d_log10f->data[i0 + d_log10f->size[0]] = 1.0;
              }

              i0 = c_rdb->size[0];
              c_rdb->size[0] = IndexRange->size[1];
              emxEnsureCapacity((emxArray__common *)c_rdb, i0, (int)sizeof
                                (double));
              loop_ub = IndexRange->size[1];
              for (i0 = 0; i0 < loop_ub; i0++) {
                c_rdb->data[i0] = rdb->data[(int)IndexRange->data
                  [IndexRange->size[0] * i0] - 1];
              }

              mldivide(d_log10f, c_rdb, b_x);
              phdeg->data[k] = b_x[0];
              k++;
            }

            emxFree_real_T(&c_rdb);
            emxFree_real_T(&d_log10f);
            emxFree_real_T(&IndexRange);
            emxInit_real_T1(&b_phdeg, 1);

            /*  */
            /*  */
            /*  Takes the last one, meaning the wider window that still sarisfies the slope demand */
            /*  */
            i0 = b_phdeg->size[0];
            b_phdeg->size[0] = phdeg->size[0];
            emxEnsureCapacity((emxArray__common *)b_phdeg, i0, (int)sizeof
                              (double));
            loop_ub = phdeg->size[0];
            for (i0 = 0; i0 < loop_ub; i0++) {
              b_phdeg->data[i0] = phdeg->data[i0] - Slope;
            }

            b_abs(b_phdeg, varargin_2);
            i0 = x->size[0];
            x->size[0] = varargin_2->size[0];
            emxEnsureCapacity((emxArray__common *)x, i0, (int)sizeof(boolean_T));
            loop_ub = varargin_2->size[0];
            emxFree_real_T(&b_phdeg);
            for (i0 = 0; i0 < loop_ub; i0++) {
              x->data[i0] = (varargin_2->data[i0] < tol);
            }

            k = (1 <= x->size[0]);
            idx = 0;
            b_n = x->size[0];
            exitg2 = false;
            while ((!exitg2) && (b_n > 0)) {
              if (x->data[b_n - 1]) {
                idx = 1;
                ii_data[0] = b_n;
                exitg2 = true;
              } else {
                b_n--;
              }
            }

            if (k == 1) {
              if (idx == 0) {
                k = 0;
              }
            } else {
              k = !(1 > idx);
            }

            for (i0 = 0; i0 < k; i0++) {
              k_data[i0] = (unsigned int)ii_data[i0];
            }

            for (i0 = 0; i0 < k; i0++) {
              flag_data[i0] = true;
            }

            /*  */
            /*  Take the values only if a solution was found */
            /*  */
            if ((!(k == 0)) && all(flag_data)) {
              /*  all() is required by codegen, checking isempty is required as all([]) = 1, surprisingly */
              for (i0 = 0; i0 < 1; i0++) {
                k1_data[0] = i2mean - (double)k_data[0];
              }

              i2mean += MinIndexRange;
              for (i0 = 0; i0 < 1; i0++) {
                k2_data[0] = (i2mean + (double)k_data[0]) - 1.0;
              }

              for (i0 = 0; i0 < 1; i0++) {
                StartFrequency_data[0] = fnew->data[(int)k1_data[0] - 1];
              }

              for (i0 = 0; i0 < 1; i0++) {
                EndFrequency_data[0] = fnew->data[(int)k2_data[0] - 1];
              }

              for (i0 = 0; i0 < 1; i0++) {
                SlopeResult_data[0] = phdeg->data[(int)k_data[0] - 1];
              }

              for (i0 = 0; i0 < 1; i0++) {
                k2_data[0] = (k2_data[0] - k1_data[0]) + 1.0;
              }
            } else {
              /*  */
              /*  We didn't suceed to extend the initial solution, so we return the initial solution */
              /*  */
              StartFrequency_data[0] = fnew->data[(int)i2mean - 1];
              EndFrequency_data[0] = fnew->data[(int)((i2mean + MinIndexRange) -
                1.0) - 1];
              SlopeResult_data[0] = y->data[(int)i2mean - 1];
              k2_data[0] = MinIndexRange;
            }

            /*  */
            /*  Required for codegen to create scalar output parameters */
            /*  */
            *StartFrequency = StartFrequency_data[0];
            *EndFrequency = EndFrequency_data[0];
            *SlopeResult = SlopeResult_data[0];
            *NumberOfPoints = k2_data[0];

            /* %%for  plant gain calculation; */
            i0 = x->size[0];
            x->size[0] = fnew->size[1];
            emxEnsureCapacity((emxArray__common *)x, i0, (int)sizeof(boolean_T));
            loop_ub = fnew->size[1];
            for (i0 = 0; i0 < loop_ub; i0++) {
              x->data[i0] = ((fnew->data[i0] > StartFrequency_data[0]) &&
                             (fnew->data[i0] < EndFrequency_data[0]));
            }

            khi = x->size[0];
            idx = 0;
            i0 = ii->size[0];
            ii->size[0] = x->size[0];
            emxEnsureCapacity((emxArray__common *)ii, i0, (int)sizeof(int));
            b_n = 1;
            exitg1 = false;
            while ((!exitg1) && (b_n <= khi)) {
              guard1 = false;
              if (x->data[b_n - 1]) {
                idx++;
                ii->data[idx - 1] = b_n;
                if (idx >= khi) {
                  exitg1 = true;
                } else {
                  guard1 = true;
                }
              } else {
                guard1 = true;
              }

              if (guard1) {
                b_n++;
              }
            }

            if (x->size[0] == 1) {
              if (idx == 0) {
                i0 = ii->size[0];
                ii->size[0] = 0;
                emxEnsureCapacity((emxArray__common *)ii, i0, (int)sizeof(int));
              }
            } else {
              i0 = ii->size[0];
              if (1 > idx) {
                ii->size[0] = 0;
              } else {
                ii->size[0] = idx;
              }

              emxEnsureCapacity((emxArray__common *)ii, i0, (int)sizeof(int));
            }

            i0 = phdeg->size[0];
            phdeg->size[0] = ii->size[0];
            emxEnsureCapacity((emxArray__common *)phdeg, i0, (int)sizeof(double));
            loop_ub = ii->size[0];
            for (i0 = 0; i0 < loop_ub; i0++) {
              phdeg->data[i0] = ii->data[i0];
            }

            if (phdeg->size[0] == 0) {
              *nErrCode = -10.0;

              /* program inside wrong;calculate failed; */
            } else if (phdeg->size[0] > 10000000) {
              *nErrCode = -11.0;

              /* freq range is out of middle variable size; */
            } else {
              for (b_n = 0; b_n < phdeg->size[0]; b_n++) {
                gainLineRange[b_n] = rt_powd_snf(10.0, rdb->data[(int)
                  phdeg->data[b_n] - 1] / 20.0) * fnew->data[(int)phdeg->
                  data[b_n] - 1] * 2.0 * 3.1415926535897931 * fnew->data[(int)
                  phdeg->data[b_n] - 1] * 2.0 * 3.1415926535897931;
              }

              *PlantGain = sum(gainLineRange) / (double)phdeg->size[0];

              /* %can not use mean of matlab for the zero value */
            }
          }
        }

        emxFree_real_T(&varargin_2);
        emxFree_real_T(&kk);
        emxFree_real_T(&y);
      }

      emxFree_int32_T(&ii);
      emxFree_boolean_T(&x);
      emxFree_real_T(&b_log10f);
      emxFree_real_T(&idxPhLarger180);
      emxFree_real_T(&phdeg);
      emxFree_real_T(&rdb);
      emxFree_real_T(&fnew);
    }
  }
}

/*
 * File trailer for PGain_from_LinearSlope.c
 *
 * [EOF]
 */
