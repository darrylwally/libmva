/*
* Filename: mvprprocess.h
* Description: Contains structures and functions for preprocessing a matrix.
* Author: Darryl Wallace <wallacdj@gmail.com
* Copyright (c) 2014 - Darryl Wallace
*
* License:
* libmva is subject to the terms of the Mozilla Public
* License, v. 2.0. If a copy of the MPL was not distributed with this
* file, You can obtain one at http://mozilla.org/MPL/2.0/.
*
* libmva is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
*/

#ifndef MVPREPROCESS_H
#define MVPREPROCESS_H

#include "mvmatrix.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef enum MVPreprocessTransformation_e {
    MV_TRANSFORM_NONE = 0,          /*! None */
    MV_TRANSFORM_LINEAR,            /*! Ax + B      */
    MV_TRANSFORM_LOGARITHMIC,       /*! log(Ax + B) */
    MV_TRANSFORM_EXPONENTIAL,       /*! e^(Ax+ B)   */
    MV_TRANSFORM_POWER              /*! (Ax + B)^C  */

} MVPreprocessTransformation;

typedef struct MVPreprocessCoefficients_s {
    double A;
    double B;
    double C;
}MVPreprocessCoeffs;

typedef enum MVPreprocessCentering_e {
    MV_CENTERING_NONE = 0,          /*! No centering */
    MV_CENTERING_MEAN               /*! Mean centering */
} MVPreprocessCentering;

typedef enum MVPreprocessScaling_e {
    MV_SCALING_NONE = 0,            /*! No scaling */
    MV_SCALING_UV,                  /*! Unit Variance (1 / std-dev)*/
    MV_SCALING_PARETO               /*! Pareto scaling (1 / sqrt(std-dev))*/

} MVPreprocessScaling;

typedef struct MVPreprocessColumn_s {
    MVPreprocessTransformation t;
    MVPreprocessCoeffs t_coeffs;
    MVPreprocessCentering c;
    MVPreprocessScaling s;
    double multiplier;

    /* These items are "private" and are set automatically when calling
      mvpreprocess_set_column() */
    MVMAT_FUNC_PTR t_func;
    MVMAT_FUNC_PTR inv_t_func;

} MVPreprocessColumnInfo;

typedef struct MVPreprocessContext_s {
    MVPreprocessColumnInfo * preprocess_info; /*! Array of size (1 x M) */
    MVMat * centering;                    /*! Array of size (1 X M) */
    MVMat * scaling;                      /*! Array of size (1 x M) */

} MVPreprocessContext;

MVPreprocessContext * mvpreprocess_alloc();
MVPreprocessContext * mvpreprocess_alloc_init_mat(int ncolumns);
MVPreprocessContext * mvpreprocess_alloc_init_mcuv(MVMat *matrix);

int mvpreprocess_init_mat(MVPreprocessContext *ctx, int ncolumns);
int mvpreprocess_init_mcuv(MVPreprocessContext *ctx, MVMat *matrix);

int mvpreprocess_free(MVPreprocessContext **ctx);

int mvpreprocess_set_column(MVPreprocessContext *ctx, int column,
                            const MVPreprocessColumnInfo *col_info);
int mvpreprocess_get_column(MVPreprocessContext *ctx, MVPreprocessColumnInfo *col_info_out,
                            int column);

/*! Prepares the preprocess context
    The centering and scaling row vectors are computed based on the columns
    of ``matrix``.

    The centering and scaling are computed after any transformations are set.

    This step must be performed before calling mvpreprocess_do or
    mvprepreocess_undo.

    \arg ctx The MVPreprocessContext
    \arg matrix The matrix that the context is prepared with
    \return MVReturnCode

    \sa mvpreprocess_do, mvpreprocess_undo
  */
int mvpreprocess_prep(MVPreprocessContext *ctx, MVMat *matrix);

/*! Performs the preprocessing
    Preprocessing is done in the following order:

        - Functional Transformations
        - Subtract Centering
        - Divide by Scaling
        - Multiplier
  */
int mvpreprocess_do(MVPreprocessContext *ctx, MVMat *preprocessed_out, MVMat *raw_in);

/*! Unpreprocesses a matrix

    Unpreprocessing is done in the following order:

        - Divide by multipler
        - Multiply by scaling
        - Add centering
        - Functional inverse transformation
  */
int mvpreprocess_undo(MVPreprocessContext *ctx, MVMat *raw_out, MVMat *preprocessed_in);


/* Function pointers for transformations */

/*! Linear transformation function

  Performs y = Ax + B

  \arg x the input value
  \arg opaque pointer to an MVPreprocessCoeffs
  \return preprocessed value or NaN if transformation fails
  \sa MVPreprocessCoeffs
  */
double mvpreprocess_linear(double x, void *opaque);


/*! Inverse Linear transformation function

  Performs x = (y - B) / A
  \arg y the input value
  \arg opaque pointer to an MVPreprocessCoeffs
  \return preprocessed value or NaN if transformation fails
  \sa MVPreprocessCoeffs
  */
double mvpreprocess_invlinear(double y, void *opaque);

/*! Log transformation function

  Performs y = log(Ax + B)
  \arg y the input value
  \arg opaque pointer to an MVPreprocessCoeffs
  \return preprocessed value or NaN if transformation fails
  \sa MVPreprocessCoeffs
  */
double mvpreprocess_log(double x, void *opaque);

/*! Inverse Log transformation function

  Performs x = (10**y - B)/A
  \arg y the input value
  \arg opaque pointer to an MVPreprocessCoeffs
  \return preprocessed value or NaN if transformation fails
  \sa MVPreprocessCoeffs
  */
double mvpreprocess_invlog(double y, void *opaque);


/*! Exponential function

  Performes y = e^(Ax + B)
  \arg x the input value
  \arg opaque pointer to an MVPreprocessCoeffs
  \return preprocessed value or NaN if transformation fails
  \sa MVPreprocessCoeffs
  */
double mvpreprocess_exp(double x, void *opaque);


/*! Inverse Exponential function

  Performes x = (ln(y) - B) / A
  \arg y the input value
  \arg opaque pointer to an MVPreprocessCoeffs
  \return preprocessed value or NaN if transformation fails
  \sa MVPreprocessCoeffs
  */
double mvpreprocess_invexp(double y, void *opaque);


/*! Power function

  Performes y = (Ax + B)^C
  \arg x the input value
  \arg opaque pointer to an MVPreprocessCoeffs
  \return preprocessed value or NaN if transformation fails
  \sa MVPreprocessCoeffs
  */
double mvpreprocess_power(double x, void *opaque);


/*! Inverse power function

  Performes x = (y^(1/C) - B) / A
  \arg y the input value
  \arg opaque pointer to an MVPreprocessCoeffs
  \return preprocessed value or NaN if transformation fails
  \sa MVPreprocessCoeffs
  */
double mvpreprocess_invpower(double y, void *opaque);


#endif // MVPREPROCESS_H
