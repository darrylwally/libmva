/*! \file mvconstants.h
  \brief Contains various constants used in mvlib.
  \author Darryl Wallace - wallacdj@gmail.com, darryl@wallynet.ca
  \date 2012

  THIS CODE IS COPYRIGHT DARRYL WALLACE (c) 2012.
  UNAUTHORIZED DISTRIBUTION, USE, OR VIEWING THE CONTENTS OF THIS FILE
  IS STRICTLY PROHIBITED.
  IF YOU HAVE RECEIVED THIS CODE WITHOUT AUTHORIZATION IT MUST BE DESTROYED
  IMMEDIATELY
  */

#ifndef MVCONSTANTS_H
#define MVCONSTANTS_H

#ifdef __cplusplus
extern "C"{
#endif

#include <math.h>

#define DATA_MISSING 0
#define DATA_PRESENT 1

typedef enum MVReturnCode_enum {
    MV_ERROR = -1,
    MV_INCORRECT_DIMENSIONS = -2,
    MV_INDEX_OUT_OF_BOUNDS = -3,
    MV_ATTEMPT_TO_FREE_NULL_MATRIX = -4,
    MV_REFERENCE_MAT_REQUIRED = -5,
    MV_WRONG_MODEL_TYPE = -6,
    MV_UNKNOWN_MODEL_TYPE = -7,
    MV_UNKOWN_CROSSVALIDATION_TYPE = -8,
    MV_SUCCESS = 0,
    MV_CROSSVAL_RULE1 = 1,
    MV_CROSSVAL_RULE2 = 2,
    MV_CROSSVAL_RULE3 = 3,
    MV_CROSSVAL_RULE4 = 4
} MVReturnCode;

#ifdef WIN32
#include <float.h>
#define MVISNAN_FUNC _isnan
#else
#define MVISNAN_FUNC isnan
#endif

// EPS
extern const double MV_EPS;
extern const double MV_SQRT_EPS;
extern const float MV_EPS_F;
extern const float MV_SQRT_EPS_F;


// NANS
double mv_NaN();
float mv_NaNf();

// Infinite
double mv_inf();
double mv_neg_inf();
float mv_inff();
float mv_neg_inff();

#ifdef __cplusplus
}
#endif

#endif // MVCONSTANTS_H
