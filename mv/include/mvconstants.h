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

typedef enum MVReturnCode_enum {
    INCORRECT_DIMENSIONS = -1,
    INDEX_OUT_OF_BOUNDS = -2,
    ATTEMPT_TO_FREE_NULL_MATRIX = -2,
    REFERENCE_MAT_REQUIRED = -3,
    WRONG_MODEL_TYPE = -4,
    UNKNOWN_MODEL_TYPE = -5,
    UNKOWN_CROSSVALIDATION_TYPE = -6,
    SUCCESS = 0,
    CROSSVAL_RULE1 = 1,
    CROSSVAL_RULE2 = 2,
    CROSSVAL_RULE3 = 3,
    CROSSVAL_RULE4 = 4
} MVReturnCode;

// EPS
extern const double MV_EPS;
extern const double MV_SQRT_EPS;
extern const float MV_EPS_F;
extern const float MV_SQRT_EPS_F;


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
