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

#define DATA_MISSING 0
#define DATA_PRESENT 1

typedef enum MVErrorCode_enum {
    INCORRECT_DIMENSIONS = -1,
    INDEX_OUT_OF_BOUNDS = -2,
    ATTEMPT_TO_FREE_NULL_MATRIX = -2,
    REFERENCE_MAT_REQUIRED = -3,
    WRONG_MODEL_TYPE = -4,
    UNKNOWN_MODEL_TYPE = -5,
    UNKOWN_CROSSVALIDATION_TYPE = -6,
    SUCCESS = 0
} MVErrorCode;

// EPS
#define MV_DBL_EPS 2.2204460492503131e-016
#define MV_DBL_SQRT_EPS 1.4901161193847656e-08
#define MV_FLOAT_EPS 1.192092896e-07
#define MV_FLOAT_SQRT_EPS 3.452669831e-04

// NANS
double mvNaN();
float mvNaNf();

// Infinite
double mvInf();
double mvNegInf();
float mvInff();
float mvNegInff();

#ifdef __cplusplus
}
#endif

#endif // MVCONSTANTS_H
