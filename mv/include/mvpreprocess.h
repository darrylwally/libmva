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
    double custom_scaling;

} MVPreprocessColumn;

typedef struct MVPreprocessContext_s {
    MVMat * matrix;                 /*! Input matrix of size (N x M) */

    /* Internal */
    MVPreprocessColumn *preprocess_info; /*Array of size (1 x M) */

} MVPreprocessContext;

MVPreprocessContext * mvpreprocess_alloc();
MVPreprocessContext * mvpreprocess_alloc_mat(MVMat *matrix);
MVPreprocessContext * mvpreprocess_alloc_mcuv(MVMat *matrix);

int mvpreprocess_init_mat(MVPreprocessContext *ctx, MVMat *matrix);
int mvpreprocess_init_mcuv(MVPreprocessContext *ctx, MVMat *matrix);

int mvpreprocess_free(MVPreprocessContext **ctx);

int mvpreprocess_set_column(MVPreprocessContext *ctx, int column,
                            MVPreprocessColumn *col_info);
int mvpreprocess_get_column(MVPreprocessContext *ctx, MVPreprocessColumn *col_info_out,
                            int column);
int mvpreprocess_do(MVPreprocessContext *ctx, MVMat *preprocessed_out, MVMat *raw_in);
int mvpreprocess_undo(MVPreprocessContext *ctx, MVMat *raw_out, MVMat *preprocessed_in);

#endif // MVPREPROCESS_H
