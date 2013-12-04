/*! \file mvpreprocess.h
  \brief Various preprocess functions.
  \author Darryl Wallace - wallacdj@gmail.com, darryl@wallynet.ca
  \date 2012

  THIS CODE IS COPYRIGHT DARRYL WALLACE (c) 2012.
  UNAUTHORIZED DISTRIBUTION, USE, OR VIEWING THE CONTENTS OF THIS FILE
  IS STRICTLY PROHIBITED.
  IF YOU HAVE RECEIVED THIS CODE WITHOUT AUTHORIZATION IT MUST BE DESTROYED
  IMMEDIATELY
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
    MV_SCALING_UV,                  /*! Unit Variance */
    MV_SCALING_PARETO               /*! Pareto scaling */

} MVPreprocessScaling;

typedef struct MVPreprocessColumn_s {
    MVPreprocessTransformation t;
    MVPreprocessCoeffs t_coeffs;
    MVPreprocessCentering c;
    MVPreprocessScaling s;
    double multiplier;

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

int mvpreprocess_prep_cs(MVPreprocessContext *ctx, MVMat *matrix);

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

#endif // MVPREPROCESS_H
