/*
 * File: PGain_from_LinearSlope.c
 *
 * MATLAB Coder version            : 2.6
 * C/C++ source code generated on  : 20-Oct-2016 19:08:04
 */

/* Include files */
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
#include "eml_sort.h"
#include "diff.h"
#include "log10.h"
#include "interp1.h"
#include "unique.h"

/* Function Declarations */
static void eml_li_find(const emxArray_boolean_T *x, emxArray_int32_T *y);
static void eml_null_assignment(emxArray_real_T *x, double idx);
static double rt_powd_snf(double u0, double u1);
static double rt_roundd_snf(double u);

/* Function Definitions */

/*
 * Arguments    : const emxArray_boolean_T *x
 *                emxArray_int32_T *y
 * Return Type  : void
 */
static void eml_li_find(const emxArray_boolean_T *x, emxArray_int32_T *y)
{
  int k;
  int i;
  int j;
  k = 0;
  for (i = 1; i <= x->size[0]; i++) {
    if (x->data[i - 1]) {
      k++;
    }
  }

  j = y->size[0];
  y->size[0] = k;
  emxEnsureCapacity((emxArray__common *)y, j, (int)sizeof(int));
  j = 0;
  for (i = 1; i <= x->size[0]; i++) {
    if (x->data[i - 1]) {
      y->data[j] = i;
      j++;
    }
  }
}

/*
 * Arguments    : emxArray_real_T *x
 *                double idx
 * Return Type  : void
 */
static void eml_null_assignment(emxArray_real_T *x, double idx)
{
  int nxin;
  int nxout;
  int k;
  emxArray_real_T *b_x;
  nxin = x->size[0];
  nxout = x->size[0] - 1;
  for (k = (int)idx; k < nxin; k++) {
    x->data[k - 1] = x->data[k];
  }

  b_emxInit_real_T(&b_x, 1);
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
  double x;
  double apnd;
  double ndbl;
  double cdiff;
  double absa;
  double absb;
  emxArray_real_T *fnew;
  int i0;
  int nm1d2;
  int k;
  double inf_data[1000];
  int inf_size[1];
  emxArray_int32_T *ii;
  emxArray_real_T b_inf_data;
  int khi;
  int idx_data[1000];
  double b_data[1000];
  int nNaN;
  int nb;
  int32_T exitg10;
  int i2mean;
  boolean_T p;
  double b_b_data[1000];
  int pos_data[1000];
  int index_f_afterUnique_data[1000];
  int ii_data[1];
  double r_data[1000];
  int r_size[1];
  emxArray_real_T *rdb;
  int b_fnew[1];
  int b_size[1];
  int vPh_size[1];
  emxArray_real_T *phdeg;
  emxArray_boolean_T *b_x;
  int c_fnew[1];
  boolean_T exitg9;
  boolean_T guard5 = false;
  emxArray_int32_T *b_ii;
  emxArray_real_T *idxPhLarger180;
  boolean_T exitg8;
  boolean_T exitg7;
  boolean_T guard4 = false;
  emxArray_int32_T *c_ii;
  emxArray_real_T *b_log10f;
  double tol;
  double MinIndexRange;
  emxArray_real_T *y;
  double d0;
  emxArray_real_T *c_x;
  emxArray_real_T *c_log10f;
  emxArray_real_T *b_rdb;
  emxArray_real_T *d_log10f;
  emxArray_real_T *r0;
  emxArray_real_T *b_y;
  boolean_T exitg6;
  boolean_T guard3 = false;
  emxArray_int32_T *d_ii;
  emxArray_real_T *kk;
  emxArray_boolean_T *d_x;
  emxArray_real_T *b_kk;
  boolean_T exitg5;
  boolean_T guard2 = false;
  emxArray_int32_T *e_ii;
  emxArray_real_T *c_y;
  emxArray_real_T *b_phdeg;
  emxArray_real_T *c_kk;
  emxArray_real_T *d_kk;
  boolean_T exitg4;
  double l[3];
  boolean_T exitg3;
  emxArray_real_T *IndexRange;
  emxArray_real_T *e_log10f;
  emxArray_real_T *c_rdb;
  emxArray_real_T *f_log10f;
  emxArray_real_T *r1;
  unsigned int b_absa;
  double u0;
  emxArray_real_T *c_phdeg;
  boolean_T exitg2;
  unsigned int k_data[1];
  double k1_data[1];
  double k2_data[1];
  double StartFrequency_data[1];
  double EndFrequency_data[1];
  double SlopeResult_data[1];
  boolean_T exitg1;
  boolean_T guard1 = false;
  emxArray_int32_T *f_ii;

  /* %nErrCode */
  /*  -1 : input parameter wrong: MinFreqRange or MaxFreqRange outof range; */
  /*  -2 : input parameter (MinFreqRange> MaxFreqRange) conflict; */
  /*  -3: not find any period within the tolorance; */
  /*  -5; calculation fail; */
  /*  -10: calculate failed; */
  /*  -11: freq range is out of middle variable size; */
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
        b_n = 0;
        x = rtNaN;
        apnd = inf[(int)n - 1];
      } else if (inf[(int)n - 1] < inf[0]) {
        b_n = -1;
        x = inf[0];
        apnd = inf[(int)n - 1];
      } else if (rtIsInf(inf[0]) || rtIsInf(inf[(int)n - 1])) {
        b_n = 0;
        x = rtNaN;
        apnd = inf[(int)n - 1];
      } else {
        x = inf[0];
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
          b_n = (int)ndbl - 1;
        } else {
          b_n = -1;
        }
      }

      emxInit_real_T(&fnew, 2);
      i0 = fnew->size[0] * fnew->size[1];
      fnew->size[0] = 1;
      fnew->size[1] = b_n + 1;
      emxEnsureCapacity((emxArray__common *)fnew, i0, (int)sizeof(double));
      if (b_n + 1 > 0) {
        fnew->data[0] = x;
        if (b_n + 1 > 1) {
          fnew->data[b_n] = apnd;
          nm1d2 = b_n / 2;
          for (k = 1; k < nm1d2; k++) {
            absa = (double)k * 0.05;
            fnew->data[k] = x + absa;
            fnew->data[b_n - k] = apnd - absa;
          }

          if (nm1d2 << 1 == b_n) {
            fnew->data[nm1d2] = (x + apnd) / 2.0;
          } else {
            absa = (double)nm1d2 * 0.05;
            fnew->data[nm1d2] = x + absa;
            fnew->data[nm1d2 + 1] = apnd - absa;
          }
        }
      }

      /*  rdb = interp1(f(1:n),r(1:n),fnew); */
      /* to solve f is non-unique problem; */
      inf_size[0] = loop_ub;
      for (i0 = 0; i0 < loop_ub; i0++) {
        inf_data[i0] = inf[i0];
      }

      emxInit_int32_T(&ii, 1);
      b_inf_data.data = (double *)&inf_data;
      b_inf_data.size = (int *)&inf_size;
      b_inf_data.allocatedSize = 1000;
      b_inf_data.numDimensions = 1;
      b_inf_data.canFreeData = false;
      eml_sort_idx(&b_inf_data, ii);
      khi = ii->size[0];
      for (i0 = 0; i0 < khi; i0++) {
        idx_data[i0] = ii->data[i0];
      }

      for (k = 0; k + 1 <= loop_ub; k++) {
        b_data[k] = inf[idx_data[k] - 1];
      }

      count_nonfinites(b_data, loop_ub, &k, &khi, &b_n, &nNaN);
      nb = -1;
      if (k > 0) {
        nb = 0;
      }

      khi += k;
      while (k + 1 <= khi) {
        x = b_data[k];
        nm1d2 = k;
        do {
          exitg10 = 0;
          k++;
          if (k + 1 > khi) {
            exitg10 = 1;
          } else {
            absa = fabs(x / 2.0);
            if ((!rtIsInf(absa)) && (!rtIsNaN(absa))) {
              if (absa <= 2.2250738585072014E-308) {
                absa = 4.94065645841247E-324;
              } else {
                frexp(absa, &i2mean);
                absa = ldexp(1.0, i2mean - 53);
              }
            } else {
              absa = rtNaN;
            }

            if ((fabs(x - b_data[k]) < absa) || (rtIsInf(b_data[k]) && rtIsInf(x)
                 && ((b_data[k] > 0.0) == (x > 0.0)))) {
              p = true;
            } else {
              p = false;
            }

            if (!p) {
              exitg10 = 1;
            }
          }
        } while (exitg10 == 0);

        nb++;
        b_data[nb] = x;
        idx_data[nb] = idx_data[nm1d2];
      }

      if (b_n > 0) {
        nb++;
        b_data[nb] = b_data[khi];
        idx_data[nb] = idx_data[khi];
      }

      k = khi + b_n;
      for (b_n = 0; b_n + 1 <= nNaN; b_n++) {
        nb++;
        b_data[nb] = b_data[k + b_n];
        idx_data[nb] = idx_data[k + b_n];
      }

      if (1 > nb + 1) {
        loop_ub = -1;
      } else {
        loop_ub = nb;
      }

      for (i0 = 0; i0 <= loop_ub; i0++) {
        b_b_data[i0] = b_data[i0];
      }

      loop_ub++;
      for (i0 = 0; i0 < loop_ub; i0++) {
        b_data[i0] = b_b_data[i0];
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

      ii_data[0] = loop_ub;
      for (i0 = 0; i0 < loop_ub; i0++) {
        b_b_data[i0] = b_data[i0];
      }

      r_size[0] = khi;
      for (i0 = 0; i0 < khi; i0++) {
        r_data[i0] = r[index_f_afterUnique_data[i0] - 1];
      }

      b_emxInit_real_T(&rdb, 1);
      b_fnew[0] = fnew->size[1];
      b_inf_data = *fnew;
      b_inf_data.size = (int *)&b_fnew;
      b_inf_data.numDimensions = 1;
      interp1(b_b_data, ii_data, r_data, r_size, &b_inf_data, rdb);

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

      b_size[0] = loop_ub;
      for (i0 = 0; i0 < loop_ub; i0++) {
        b_b_data[i0] = b_data[i0];
      }

      vPh_size[0] = khi;
      for (i0 = 0; i0 < khi; i0++) {
        b_data[i0] = vPh[i0];
      }

      b_emxInit_real_T(&phdeg, 1);
      emxInit_boolean_T(&b_x, 1);
      c_fnew[0] = fnew->size[1];
      b_inf_data = *fnew;
      b_inf_data.size = (int *)&c_fnew;
      b_inf_data.numDimensions = 1;
      interp1(b_b_data, b_size, b_data, vPh_size, &b_inf_data, phdeg);

      /* %improve program,add vPh info */
      /* %%not use the ending freq;but use the vPh=-180deg to set the EndIndex; */
      i0 = b_x->size[0];
      b_x->size[0] = phdeg->size[0];
      emxEnsureCapacity((emxArray__common *)b_x, i0, (int)sizeof(boolean_T));
      loop_ub = phdeg->size[0];
      for (i0 = 0; i0 < loop_ub; i0++) {
        b_x->data[i0] = (phdeg->data[i0] <= -240.0);
      }

      khi = 0;
      i0 = ii->size[0];
      ii->size[0] = b_x->size[0];
      emxEnsureCapacity((emxArray__common *)ii, i0, (int)sizeof(int));
      b_n = 1;
      exitg9 = false;
      while ((!exitg9) && (b_n <= b_x->size[0])) {
        guard5 = false;
        if (b_x->data[b_n - 1]) {
          khi++;
          ii->data[khi - 1] = b_n;
          if (khi >= b_x->size[0]) {
            exitg9 = true;
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

      if (b_x->size[0] == 1) {
        if (khi == 0) {
          i0 = ii->size[0];
          ii->size[0] = 0;
          emxEnsureCapacity((emxArray__common *)ii, i0, (int)sizeof(int));
        }
      } else {
        if (1 > khi) {
          loop_ub = 0;
        } else {
          loop_ub = khi;
        }

        emxInit_int32_T(&b_ii, 1);
        i0 = b_ii->size[0];
        b_ii->size[0] = loop_ub;
        emxEnsureCapacity((emxArray__common *)b_ii, i0, (int)sizeof(int));
        for (i0 = 0; i0 < loop_ub; i0++) {
          b_ii->data[i0] = ii->data[i0];
        }

        i0 = ii->size[0];
        ii->size[0] = b_ii->size[0];
        emxEnsureCapacity((emxArray__common *)ii, i0, (int)sizeof(int));
        loop_ub = b_ii->size[0];
        for (i0 = 0; i0 < loop_ub; i0++) {
          ii->data[i0] = b_ii->data[i0];
        }

        emxFree_int32_T(&b_ii);
      }

      b_emxInit_real_T(&idxPhLarger180, 1);
      i0 = idxPhLarger180->size[0];
      idxPhLarger180->size[0] = ii->size[0];
      emxEnsureCapacity((emxArray__common *)idxPhLarger180, i0, (int)sizeof
                        (double));
      loop_ub = ii->size[0];
      for (i0 = 0; i0 < loop_ub; i0++) {
        idxPhLarger180->data[i0] = ii->data[i0];
      }

      if (idxPhLarger180->size[0] == 0) {
        khi = 1;
        x = phdeg->data[0];
        if (phdeg->size[0] > 1) {
          if (rtIsNaN(phdeg->data[0])) {
            b_n = 2;
            exitg8 = false;
            while ((!exitg8) && (b_n <= phdeg->size[0])) {
              khi = b_n;
              if (!rtIsNaN(phdeg->data[b_n - 1])) {
                x = phdeg->data[b_n - 1];
                exitg8 = true;
              } else {
                b_n++;
              }
            }
          }

          if (khi < phdeg->size[0]) {
            while (khi + 1 <= phdeg->size[0]) {
              if (phdeg->data[khi] < x) {
                x = phdeg->data[khi];
              }

              khi++;
            }
          }
        }

        i0 = b_x->size[0];
        b_x->size[0] = phdeg->size[0];
        emxEnsureCapacity((emxArray__common *)b_x, i0, (int)sizeof(boolean_T));
        loop_ub = phdeg->size[0];
        for (i0 = 0; i0 < loop_ub; i0++) {
          b_x->data[i0] = (phdeg->data[i0] == x);
        }

        khi = 0;
        i0 = ii->size[0];
        ii->size[0] = b_x->size[0];
        emxEnsureCapacity((emxArray__common *)ii, i0, (int)sizeof(int));
        b_n = 1;
        exitg7 = false;
        while ((!exitg7) && (b_n <= b_x->size[0])) {
          guard4 = false;
          if (b_x->data[b_n - 1]) {
            khi++;
            ii->data[khi - 1] = b_n;
            if (khi >= b_x->size[0]) {
              exitg7 = true;
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

        if (b_x->size[0] == 1) {
          if (khi == 0) {
            i0 = ii->size[0];
            ii->size[0] = 0;
            emxEnsureCapacity((emxArray__common *)ii, i0, (int)sizeof(int));
          }
        } else {
          if (1 > khi) {
            loop_ub = 0;
          } else {
            loop_ub = khi;
          }

          emxInit_int32_T(&c_ii, 1);
          i0 = c_ii->size[0];
          c_ii->size[0] = loop_ub;
          emxEnsureCapacity((emxArray__common *)c_ii, i0, (int)sizeof(int));
          for (i0 = 0; i0 < loop_ub; i0++) {
            c_ii->data[i0] = ii->data[i0];
          }

          i0 = ii->size[0];
          ii->size[0] = c_ii->size[0];
          emxEnsureCapacity((emxArray__common *)ii, i0, (int)sizeof(int));
          loop_ub = c_ii->size[0];
          for (i0 = 0; i0 < loop_ub; i0++) {
            ii->data[i0] = c_ii->data[i0];
          }

          emxFree_int32_T(&c_ii);
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

      b_emxInit_real_T(&b_log10f, 1);

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
        b_emxInit_real_T(&y, 1);

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
        b_emxInit_real_T(&c_x, 1);
        emxInit_real_T(&c_log10f, 2);
        b_emxInit_real_T(&b_rdb, 1);
        b_emxInit_real_T(&d_log10f, 1);
        b_emxInit_real_T(&r0, 1);
        while (k <= (int)d0 - 1) {
          absa = ((1.0 + (double)k) + MinIndexRange) - 1.0;
          if (1.0 + (double)k > absa) {
            i0 = 0;
            khi = 0;
          } else {
            i0 = k;
            khi = (int)absa;
          }

          absa = ((1.0 + (double)k) + MinIndexRange) - 1.0;
          if (1.0 + (double)k > absa) {
            b_n = 0;
            nNaN = 0;
          } else {
            b_n = k;
            nNaN = (int)absa;
          }

          i2mean = khi - i0;
          nm1d2 = d_log10f->size[0];
          d_log10f->size[0] = khi - i0;
          emxEnsureCapacity((emxArray__common *)d_log10f, nm1d2, (int)sizeof
                            (double));
          loop_ub = khi - i0;
          for (khi = 0; khi < loop_ub; khi++) {
            d_log10f->data[khi] = b_log10f->data[i0 + khi];
          }

          i0 = r0->size[0];
          r0->size[0] = (int)MinIndexRange;
          emxEnsureCapacity((emxArray__common *)r0, i0, (int)sizeof(double));
          loop_ub = (int)MinIndexRange;
          for (i0 = 0; i0 < loop_ub; i0++) {
            r0->data[i0] = 1.0;
          }

          i0 = c_log10f->size[0] * c_log10f->size[1];
          c_log10f->size[0] = i2mean;
          c_log10f->size[1] = 2;
          emxEnsureCapacity((emxArray__common *)c_log10f, i0, (int)sizeof(double));
          for (i0 = 0; i0 < i2mean; i0++) {
            c_log10f->data[i0] = d_log10f->data[i0];
          }

          loop_ub = (int)MinIndexRange;
          for (i0 = 0; i0 < loop_ub; i0++) {
            c_log10f->data[i0 + c_log10f->size[0]] = r0->data[i0];
          }

          i0 = b_rdb->size[0];
          b_rdb->size[0] = nNaN - b_n;
          emxEnsureCapacity((emxArray__common *)b_rdb, i0, (int)sizeof(double));
          loop_ub = nNaN - b_n;
          for (i0 = 0; i0 < loop_ub; i0++) {
            b_rdb->data[i0] = rdb->data[b_n + i0];
          }

          mldivide(c_log10f, b_rdb, c_x);
          y->data[k] = c_x->data[0];
          k++;
        }

        emxFree_real_T(&r0);
        emxFree_real_T(&d_log10f);
        emxFree_real_T(&b_rdb);
        emxFree_real_T(&c_log10f);
        b_emxInit_real_T(&b_y, 1);

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

        b_abs(b_y, c_x);
        i0 = b_x->size[0];
        b_x->size[0] = c_x->size[0];
        emxEnsureCapacity((emxArray__common *)b_x, i0, (int)sizeof(boolean_T));
        loop_ub = c_x->size[0];
        emxFree_real_T(&b_y);
        for (i0 = 0; i0 < loop_ub; i0++) {
          b_x->data[i0] = (c_x->data[i0] < tol);
        }

        khi = 0;
        i0 = ii->size[0];
        ii->size[0] = b_x->size[0];
        emxEnsureCapacity((emxArray__common *)ii, i0, (int)sizeof(int));
        b_n = 1;
        exitg6 = false;
        while ((!exitg6) && (b_n <= b_x->size[0])) {
          guard3 = false;
          if (b_x->data[b_n - 1]) {
            khi++;
            ii->data[khi - 1] = b_n;
            if (khi >= b_x->size[0]) {
              exitg6 = true;
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

        if (b_x->size[0] == 1) {
          if (khi == 0) {
            i0 = ii->size[0];
            ii->size[0] = 0;
            emxEnsureCapacity((emxArray__common *)ii, i0, (int)sizeof(int));
          }
        } else {
          if (1 > khi) {
            loop_ub = 0;
          } else {
            loop_ub = khi;
          }

          emxInit_int32_T(&d_ii, 1);
          i0 = d_ii->size[0];
          d_ii->size[0] = loop_ub;
          emxEnsureCapacity((emxArray__common *)d_ii, i0, (int)sizeof(int));
          for (i0 = 0; i0 < loop_ub; i0++) {
            d_ii->data[i0] = ii->data[i0];
          }

          i0 = ii->size[0];
          ii->size[0] = d_ii->size[0];
          emxEnsureCapacity((emxArray__common *)ii, i0, (int)sizeof(int));
          loop_ub = d_ii->size[0];
          for (i0 = 0; i0 < loop_ub; i0++) {
            ii->data[i0] = d_ii->data[i0];
          }

          emxFree_int32_T(&d_ii);
        }

        b_emxInit_real_T(&kk, 1);
        i0 = kk->size[0];
        kk->size[0] = ii->size[0];
        emxEnsureCapacity((emxArray__common *)kk, i0, (int)sizeof(double));
        loop_ub = ii->size[0];
        for (i0 = 0; i0 < loop_ub; i0++) {
          kk->data[i0] = ii->data[i0];
        }

        emxInit_boolean_T(&d_x, 1);
        diff(kk, c_x);
        i0 = d_x->size[0];
        d_x->size[0] = c_x->size[0];
        emxEnsureCapacity((emxArray__common *)d_x, i0, (int)sizeof(boolean_T));
        loop_ub = c_x->size[0];
        for (i0 = 0; i0 < loop_ub; i0++) {
          d_x->data[i0] = (c_x->data[i0] == 1.0);
        }

        b_emxInit_real_T(&b_kk, 1);
        eml_li_find(d_x, ii);
        i0 = b_kk->size[0];
        b_kk->size[0] = ii->size[0];
        emxEnsureCapacity((emxArray__common *)b_kk, i0, (int)sizeof(double));
        loop_ub = ii->size[0];
        emxFree_boolean_T(&d_x);
        for (i0 = 0; i0 < loop_ub; i0++) {
          b_kk->data[i0] = kk->data[ii->data[i0] - 1];
        }

        i0 = kk->size[0];
        kk->size[0] = b_kk->size[0];
        emxEnsureCapacity((emxArray__common *)kk, i0, (int)sizeof(double));
        loop_ub = b_kk->size[0];
        for (i0 = 0; i0 < loop_ub; i0++) {
          kk->data[i0] = b_kk->data[i0];
        }

        emxFree_real_T(&b_kk);
        if (kk->size[0] == 0) {
          *nErrCode = -3.0;

          /* %tolerance too tighten,pls increase tolerance */
        } else {
          /*  */
          /*  jj - is set of indexes of kk. kk(jj(1)), kk(jj(2), ..., kk(jj(N) */
          /*  are the indexes that point to the starting and the ending of sequence of kk. */
          diff(kk, c_x);
          i0 = b_x->size[0];
          b_x->size[0] = c_x->size[0];
          emxEnsureCapacity((emxArray__common *)b_x, i0, (int)sizeof(boolean_T));
          loop_ub = c_x->size[0];
          for (i0 = 0; i0 < loop_ub; i0++) {
            b_x->data[i0] = (c_x->data[i0] != 1.0);
          }

          khi = 0;
          i0 = ii->size[0];
          ii->size[0] = b_x->size[0];
          emxEnsureCapacity((emxArray__common *)ii, i0, (int)sizeof(int));
          b_n = 1;
          exitg5 = false;
          while ((!exitg5) && (b_n <= b_x->size[0])) {
            guard2 = false;
            if (b_x->data[b_n - 1]) {
              khi++;
              ii->data[khi - 1] = b_n;
              if (khi >= b_x->size[0]) {
                exitg5 = true;
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

          if (b_x->size[0] == 1) {
            if (khi == 0) {
              i0 = ii->size[0];
              ii->size[0] = 0;
              emxEnsureCapacity((emxArray__common *)ii, i0, (int)sizeof(int));
            }
          } else {
            if (1 > khi) {
              loop_ub = 0;
            } else {
              loop_ub = khi;
            }

            emxInit_int32_T(&e_ii, 1);
            i0 = e_ii->size[0];
            e_ii->size[0] = loop_ub;
            emxEnsureCapacity((emxArray__common *)e_ii, i0, (int)sizeof(int));
            for (i0 = 0; i0 < loop_ub; i0++) {
              e_ii->data[i0] = ii->data[i0];
            }

            i0 = ii->size[0];
            ii->size[0] = e_ii->size[0];
            emxEnsureCapacity((emxArray__common *)ii, i0, (int)sizeof(int));
            loop_ub = e_ii->size[0];
            for (i0 = 0; i0 < loop_ub; i0++) {
              ii->data[i0] = e_ii->data[i0];
            }

            emxFree_int32_T(&e_ii);
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
            b_emxInit_real_T(&c_y, 1);

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
            i0 = kk->size[0];
            b_n = 0;
            emxFree_real_T(&c_y);
            while (b_n <= i0 - 1) {
              gainLineRange[b_n] = rt_powd_snf(10.0, rdb->data[(int)kk->data[b_n]
                - 1] / 20.0) * fnew->data[(int)kk->data[b_n] - 1] * 2.0 *
                3.1415926535897931 * fnew->data[(int)kk->data[b_n] - 1] * 2.0 *
                3.1415926535897931;
              b_n++;
            }

            *PlantGain = sum(gainLineRange) / (double)kk->size[0];

            /* %can not use mean of matlab for the zero value */
            /* %%finish cal plant gain */
          } else {
            b_emxInit_real_T(&c_y, 1);
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
            b_emxInit_real_T(&b_phdeg, 1);
            i0 = b_phdeg->size[0];
            b_phdeg->size[0] = phdeg->size[0];
            emxEnsureCapacity((emxArray__common *)b_phdeg, i0, (int)sizeof
                              (double));
            loop_ub = phdeg->size[0];
            for (i0 = 0; i0 < loop_ub; i0++) {
              b_phdeg->data[i0] = phdeg->data[i0];
            }

            b_emxInit_real_T(&c_kk, 1);
            eml_sort(b_phdeg, phdeg);
            i0 = phdeg->size[0];
            eml_null_assignment(phdeg, i0);

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
            khi = c_kk->size[0];
            c_kk->size[0] = (i0 >> 1) + 1;
            emxEnsureCapacity((emxArray__common *)c_kk, khi, (int)sizeof(double));
            loop_ub = i0 >> 1;
            emxFree_real_T(&b_phdeg);
            for (i0 = 0; i0 <= loop_ub; i0++) {
              c_kk->data[i0] = kk->data[(int)phdeg->data[1 + (i0 << 1)] - 1];
            }

            b_emxInit_real_T(&d_kk, 1);
            i0 = (int)((double)phdeg->size[0] - 1.0) - 1;
            khi = d_kk->size[0];
            d_kk->size[0] = (i0 >> 1) + 1;
            emxEnsureCapacity((emxArray__common *)d_kk, khi, (int)sizeof(double));
            loop_ub = i0 >> 1;
            for (i0 = 0; i0 <= loop_ub; i0++) {
              d_kk->data[i0] = kk->data[(int)phdeg->data[i0 << 1] - 1];
            }

            rdivide(c_kk, d_kk, c_x);
            khi = 1;
            x = c_x->data[0];
            nm1d2 = 0;
            emxFree_real_T(&d_kk);
            emxFree_real_T(&c_kk);
            if (c_x->size[0] > 1) {
              if (rtIsNaN(c_x->data[0])) {
                b_n = 2;
                exitg4 = false;
                while ((!exitg4) && (b_n <= c_x->size[0])) {
                  khi = b_n;
                  if (!rtIsNaN(c_x->data[b_n - 1])) {
                    x = c_x->data[b_n - 1];
                    nm1d2 = b_n - 1;
                    exitg4 = true;
                  } else {
                    b_n++;
                  }
                }
              }

              if (khi < c_x->size[0]) {
                while (khi + 1 <= c_x->size[0]) {
                  if (c_x->data[khi] > x) {
                    x = c_x->data[khi];
                    nm1d2 = khi;
                  }

                  khi++;
                }
              }
            }

            /* %bug:if the tolorance too small,the ii maybe empty */
            i2mean = (int)rt_roundd_snf(sqrt(kk->data[(int)phdeg->data[nm1d2 + 1]
              - 1] * kk->data[(int)phdeg->data[nm1d2] - 1]));

            /*  Geometrical mean */
            /*  */
            /*  Calculate the permmited range that the min range can be extended to */
            /*  */
            l[0] = 2.0 * ((double)i2mean - 1.0) + MinIndexRange;

            /*  limitation by the first index of the vectors. */
            l[1] = 2.0 * (idxPhLarger180->data[0] - (((double)i2mean - 1.0) +
              MinIndexRange)) + MinIndexRange;

            /*  limitation by the last index of the vectors. */
            l[2] = MaxFreqRange / 0.05;
            b_fix(&l[2]);

            /*  limitation as set by the user. */
            khi = 1;
            x = l[0];
            if (rtIsNaN(l[0])) {
              b_n = 2;
              exitg3 = false;
              while ((!exitg3) && (b_n < 4)) {
                khi = b_n;
                if (!rtIsNaN(l[b_n - 1])) {
                  x = l[b_n - 1];
                  exitg3 = true;
                } else {
                  b_n++;
                }
              }
            }

            if (khi < 3) {
              while (khi + 1 < 4) {
                if (l[khi] < x) {
                  x = l[khi];
                }

                khi++;
              }
            }

            /*  */
            /*  Claculate the slopes of all gain responses around the optimal point */
            /*  and put the results in z vector. */
            /*  */
            absa = (x - MinIndexRange) / 2.0;
            i0 = phdeg->size[0];
            phdeg->size[0] = (int)absa;
            emxEnsureCapacity((emxArray__common *)phdeg, i0, (int)sizeof(double));
            loop_ub = (int)absa;
            for (i0 = 0; i0 < loop_ub; i0++) {
              phdeg->data[i0] = 0.0;
            }

            d0 = (x - MinIndexRange) / 2.0;
            k = 0;
            emxInit_real_T(&IndexRange, 2);
            emxInit_real_T(&e_log10f, 2);
            b_emxInit_real_T(&c_rdb, 1);
            b_emxInit_real_T(&f_log10f, 1);
            b_emxInit_real_T(&r1, 1);
            while (k <= (int)d0 - 1) {
              x = (double)i2mean - (1.0 + (double)k);
              absa = (((double)i2mean + MinIndexRange) + (1.0 + (double)k)) -
                1.0;
              if (rtIsNaN(absa)) {
                b_n = 0;
                x = rtNaN;
                apnd = absa;
              } else if (absa < x) {
                b_n = -1;
                apnd = absa;
              } else if (rtIsInf(absa)) {
                b_n = 0;
                x = rtNaN;
                apnd = absa;
              } else {
                ndbl = floor((absa - x) + 0.5);
                apnd = x + ndbl;
                cdiff = apnd - absa;
                b_absa = (unsigned int)fabs(x);
                absb = fabs(absa);
                u0 = b_absa;
                if ((u0 >= absb) || rtIsNaN(absb)) {
                  absb = u0;
                }

                if (fabs(cdiff) < 4.4408920985006262E-16 * absb) {
                  ndbl++;
                  apnd = absa;
                } else if (cdiff > 0.0) {
                  apnd = x + (ndbl - 1.0);
                } else {
                  ndbl++;
                }

                if (ndbl >= 0.0) {
                  b_n = (int)ndbl - 1;
                } else {
                  b_n = -1;
                }
              }

              i0 = IndexRange->size[0] * IndexRange->size[1];
              IndexRange->size[0] = 1;
              IndexRange->size[1] = b_n + 1;
              emxEnsureCapacity((emxArray__common *)IndexRange, i0, (int)sizeof
                                (double));
              if (b_n + 1 > 0) {
                IndexRange->data[0] = x;
                if (b_n + 1 > 1) {
                  IndexRange->data[b_n] = apnd;
                  nm1d2 = b_n / 2;
                  for (khi = 1; khi < nm1d2; khi++) {
                    IndexRange->data[khi] = x + (double)khi;
                    IndexRange->data[b_n - khi] = apnd - (double)khi;
                  }

                  if (nm1d2 << 1 == b_n) {
                    IndexRange->data[nm1d2] = (x + apnd) / 2.0;
                  } else {
                    IndexRange->data[nm1d2] = x + (double)nm1d2;
                    IndexRange->data[nm1d2 + 1] = apnd - (double)nm1d2;
                  }
                }
              }

              i0 = f_log10f->size[0];
              f_log10f->size[0] = IndexRange->size[1];
              emxEnsureCapacity((emxArray__common *)f_log10f, i0, (int)sizeof
                                (double));
              loop_ub = IndexRange->size[1];
              for (i0 = 0; i0 < loop_ub; i0++) {
                f_log10f->data[i0] = b_log10f->data[(int)IndexRange->
                  data[IndexRange->size[0] * i0] - 1];
              }

              nm1d2 = IndexRange->size[1];
              khi = IndexRange->size[1];
              b_n = IndexRange->size[1];
              i0 = r1->size[0];
              r1->size[0] = khi;
              emxEnsureCapacity((emxArray__common *)r1, i0, (int)sizeof(double));
              for (i0 = 0; i0 < khi; i0++) {
                r1->data[i0] = 1.0;
              }

              i0 = e_log10f->size[0] * e_log10f->size[1];
              e_log10f->size[0] = nm1d2;
              e_log10f->size[1] = 2;
              emxEnsureCapacity((emxArray__common *)e_log10f, i0, (int)sizeof
                                (double));
              for (i0 = 0; i0 < nm1d2; i0++) {
                e_log10f->data[i0] = f_log10f->data[i0];
              }

              for (i0 = 0; i0 < b_n; i0++) {
                e_log10f->data[i0 + e_log10f->size[0]] = r1->data[i0];
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

              mldivide(e_log10f, c_rdb, c_x);
              phdeg->data[k] = c_x->data[0];
              k++;
            }

            emxFree_real_T(&r1);
            emxFree_real_T(&f_log10f);
            emxFree_real_T(&c_rdb);
            emxFree_real_T(&e_log10f);
            emxFree_real_T(&IndexRange);
            b_emxInit_real_T(&c_phdeg, 1);

            /*  */
            /*  */
            /*  Takes the last one, meaning the wider window that still sarisfies the slope demand */
            /*  */
            i0 = c_phdeg->size[0];
            c_phdeg->size[0] = phdeg->size[0];
            emxEnsureCapacity((emxArray__common *)c_phdeg, i0, (int)sizeof
                              (double));
            loop_ub = phdeg->size[0];
            for (i0 = 0; i0 < loop_ub; i0++) {
              c_phdeg->data[i0] = phdeg->data[i0] - Slope;
            }

            b_abs(c_phdeg, c_x);
            i0 = b_x->size[0];
            b_x->size[0] = c_x->size[0];
            emxEnsureCapacity((emxArray__common *)b_x, i0, (int)sizeof(boolean_T));
            loop_ub = c_x->size[0];
            emxFree_real_T(&c_phdeg);
            for (i0 = 0; i0 < loop_ub; i0++) {
              b_x->data[i0] = (c_x->data[i0] < tol);
            }

            if (1 <= b_x->size[0]) {
              k = 1;
            } else {
              k = 0;
            }

            khi = 0;
            b_n = b_x->size[0];
            exitg2 = false;
            while ((!exitg2) && (b_n > 0)) {
              if (b_x->data[b_n - 1]) {
                khi = 1;
                ii_data[0] = b_n;
                exitg2 = true;
              } else {
                b_n--;
              }
            }

            if (k == 1) {
              if (khi == 0) {
                k = 0;
              }
            } else {
              if (1 > khi) {
                loop_ub = -1;
              } else {
                loop_ub = 0;
              }

              i0 = 0;
              while (i0 <= loop_ub) {
                r_size[0] = ii_data[0];
                i0 = 1;
              }

              k = loop_ub + 1;
              khi = loop_ub + 1;
              i0 = 0;
              while (i0 <= khi - 1) {
                ii_data[0] = r_size[0];
                i0 = 1;
              }

              i0 = loop_ub + 1;
              khi = i0 / 2;
              b_n = 1;
              while (b_n <= khi) {
                nm1d2 = ii_data[0];
                ii_data[loop_ub] = nm1d2;
                b_n = 2;
              }
            }

            i0 = 0;
            while (i0 <= k - 1) {
              k_data[0] = (unsigned int)ii_data[0];
              i0 = 1;
            }

            /*  */
            /*  Take the values only if a solution was found */
            /*  */
            if ((!(k == 0)) && all()) {
              /*  all() is required by codegen, checking isempty is required as all([]) = 1, surprisingly */
              i0 = 0;
              while (i0 <= 0) {
                k1_data[0] = (double)i2mean - (double)k_data[0];
                i0 = 1;
              }

              absa = (double)i2mean + MinIndexRange;
              i0 = 0;
              while (i0 <= 0) {
                k2_data[0] = (absa + (double)k_data[0]) - 1.0;
                i0 = 1;
              }

              i0 = 0;
              while (i0 <= 0) {
                StartFrequency_data[0] = fnew->data[(int)k1_data[0] - 1];
                i0 = 1;
              }

              i0 = 0;
              while (i0 <= 0) {
                EndFrequency_data[0] = fnew->data[(int)k2_data[0] - 1];
                i0 = 1;
              }

              i0 = 0;
              while (i0 <= 0) {
                SlopeResult_data[0] = phdeg->data[(int)k_data[0] - 1];
                i0 = 1;
              }

              i0 = 0;
              while (i0 <= 0) {
                k2_data[0] = (k2_data[0] - k1_data[0]) + 1.0;
                i0 = 1;
              }
            } else {
              /*  */
              /*  We didn't suceed to extend the initial solution, so we return the initial solution */
              /*  */
              StartFrequency_data[0] = fnew->data[i2mean - 1];
              EndFrequency_data[0] = fnew->data[(int)(((double)i2mean +
                MinIndexRange) - 1.0) - 1];
              SlopeResult_data[0] = y->data[i2mean - 1];
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
            i0 = b_x->size[0];
            b_x->size[0] = fnew->size[1];
            emxEnsureCapacity((emxArray__common *)b_x, i0, (int)sizeof(boolean_T));
            loop_ub = fnew->size[1];
            for (i0 = 0; i0 < loop_ub; i0++) {
              b_x->data[i0] = ((fnew->data[i0] > StartFrequency_data[0]) &&
                               (fnew->data[i0] < EndFrequency_data[0]));
            }

            khi = 0;
            i0 = ii->size[0];
            ii->size[0] = b_x->size[0];
            emxEnsureCapacity((emxArray__common *)ii, i0, (int)sizeof(int));
            b_n = 1;
            exitg1 = false;
            while ((!exitg1) && (b_n <= b_x->size[0])) {
              guard1 = false;
              if (b_x->data[b_n - 1]) {
                khi++;
                ii->data[khi - 1] = b_n;
                if (khi >= b_x->size[0]) {
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

            if (b_x->size[0] == 1) {
              if (khi == 0) {
                i0 = ii->size[0];
                ii->size[0] = 0;
                emxEnsureCapacity((emxArray__common *)ii, i0, (int)sizeof(int));
              }
            } else {
              if (1 > khi) {
                loop_ub = 0;
              } else {
                loop_ub = khi;
              }

              emxInit_int32_T(&f_ii, 1);
              i0 = f_ii->size[0];
              f_ii->size[0] = loop_ub;
              emxEnsureCapacity((emxArray__common *)f_ii, i0, (int)sizeof(int));
              for (i0 = 0; i0 < loop_ub; i0++) {
                f_ii->data[i0] = ii->data[i0];
              }

              i0 = ii->size[0];
              ii->size[0] = f_ii->size[0];
              emxEnsureCapacity((emxArray__common *)ii, i0, (int)sizeof(int));
              loop_ub = f_ii->size[0];
              for (i0 = 0; i0 < loop_ub; i0++) {
                ii->data[i0] = f_ii->data[i0];
              }

              emxFree_int32_T(&f_ii);
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

        emxFree_real_T(&kk);
        emxFree_real_T(&c_x);
        emxFree_real_T(&y);
      }

      emxFree_int32_T(&ii);
      emxFree_boolean_T(&b_x);
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
